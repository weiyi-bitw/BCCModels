require(rms)
require(MASS)
require(survival)
require(predictiveModeling)
require(gbm)
require(caret)
require(randomSurvivalForest)
require(devtools)

#  download.file("http://dl.dropbox.com/u/11986954/Oslo/DreamBox7_0.214.tar.gz", destfile="./DreamBox7_0.214.tar.gz")
#  install.packages("./DreamBox7_0.214.tar.gz", repos=NULL)

install_github(repo="DreamBox7", username="weiyi-bitw", ref="master")
library(DreamBox7)

setRefClass(Class = "PredictiveModel")

#' GoldiModel
#'
#' Modified from DemoClinicalOnlyModel from BCC challenge.
#'
#' @author Wei-Yi Cheng
#' @export

GoldiModel <- setRefClass(Class  = "GoldiModel",
                          contains = "PredictiveModel",
                          fields   = c("model", "attractome", "annot", "predictions", "chosenProbes"),
                          methods  = list(
                            
                            initialize = function(...){
                              	data(attractome.minimalist)
				.self$attractome = attractome.minimalist
				data(map)
                              	.self$annot = map
                              return(.self)
                            },
                            
                            customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...)
                            {
                              if(class(clinicalSurvData) != "Surv"){
                                stop("Expecting 'responseData' object of type 'Surv'")
                              }

                              # no need for cnv data, remove it
                              rm(copyData)
				clnc = lazyImputeDFClncOslo(clinicalFeaturesData)
                              	clinical <- expandClncOslo(clnc)
				clinical$tr.CT = NULL; clinical$tr.HT = NULL

			cat("Create metagene space...");flush.console()
                              o = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)
				meta = o$metaSpace
				.self$chosenProbes = o$pbs
				#meta = .self$meta.old
				#meta = t( sapply(.self$attractome, function(pbs){ apply(exprs(exprData)[pbs,] , 2, mean) }  ) )
				ls = meta["ls",]
				meta = t(apply(meta, 1, function(x){ x - median(x)}))
				#meta = meta - matrix(.self$mdns, nrow=nrow(meta), ncol=ncol(meta))
			#===== Conditioning mesenchymal transition metagene by lymph node status
				idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
				mes.lymphneg = meta["mt",] * idx
				mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
				meta = rbind(meta, mes.lymphneg)

				idx = (clinical[,"lymph_nodes_positive"] > 3)
				ls.lymphpos = ls * idx
				ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
				meta = rbind(meta, ls.lymphpos)

				idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
				ls.erneg = ls * idx
				ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
				meta = rbind(meta, ls.erneg)

				lrrc48 = exprs(exprData)["ILMN_2353862",]
				lrrc48 = lrrc48 - median(lrrc48)

				lym.N = factor(clnc$lymph_nodes_positive < 1)
				lymph = clinical[,"lymph_nodes_positive"]
				lymph = sapply(lymph, function(x){min(x, 7)})
				lsxlymph = ls * (7 - lymph)

				#meta = rbind(meta, gs)
				gpr4 = exprs(exprData)["ILMN_2074477",]

			rm(exprData)

                        cat("done!\n");flush.console()
			
#===== 1. AIC models clinical only =====
			cat("1.  Training AIC model using only clinical features...");flush.console()
				X = clinical
                                upper = terms(clinicalSurvData~(.), data = X)
				cm = step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                        cat("done!\n");flush.console()
#===== 2. GBM model =====
			cat("2.  Training gbm model...");flush.console()
				X = X[, attr(cm$term, "term.labels")]
                                cgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=8, cv.folds=5, verbose=F, seed=53) 
			cat("done!\n");flush.console()
#===== 3. AIC models metagenes only =======
                        cat("3.  Training model with metagenes ...");flush.console()
				X = data.frame(t(meta))
                              upper = terms(clinicalSurvData~(.), data = X)
                              coxmodel = step(coxph(clinicalSurvData~., data=X), scope=upper, direction="both", k=2, trace=FALSE)
			cat("done!\n");flush.console()
#===== 4. GBM model =========
			cat("4.  Training gbm model...");flush.console()
				X = data.frame(t(meta))
				gbmmodel = gbm.cvrun(clinicalSurvData~.,data=X, distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=6, cv.folds=5, verbose=F, seed=913) # my bday ;D
			cat("done!\n");flush.console()
#===== 5. KNN model ======
			cat("5.  Creating KNN database ...");flush.console()
			t = clinicalSurvData[,1]
			defSurvSamples = which(clinicalSurvData[,2]==1 | clinicalSurvData[,1] > 365 * 10)
			ccdi = getAllCCDIWz(meta, clinicalSurvData)
			idx = ccdi[c("er", "mitotic", "puf60", "erbb2", "chr7p11.2", "ls")]
			
			knnmodel = list()
			knnmodel$x.train = list(meta=meta[names(idx), defSurvSamples], time=t[defSurvSamples], concordance=idx)
			
			#clnc = lazyImputeDFClnc(clinicalFeaturesData)
			knnmodel$c.train = preproClncKNN(clnc, clinicalSurvData, ccdi.upper=0.6, ccdi.lower=0.4)

			cat("done!\n");flush.console()

#===== 6. Training cox model using minimal features ======= 
#
#  The minimalist features were selected by brute force to find the best 8 features that maximize the concordance index in the training set
#
#=====
                        cat("6.  Model A: Protective Minimalist features...");flush.console()
				X = data.frame( cbind(meta["mitotic",], meta["ls.erneg",], clinical$lymph_nodes_positive, meta["mes.lymphneg",], meta["susd3",], clinical$age_at_diagnosis) )
				colnames(X) = c("CIN", "LYM_ERNeg", "lymNum", "MES_lymNumNeg", "SUSD3", "age")
				
                              cox.a = coxph(clinicalSurvData~., data=X)
                        cat("done!\n");

#===== 7. Training GBM model using minimal features =========
			cat("7.  Training gbm model...");flush.console()
                        	#gbmmodel = gbm.fit(X,clinicalSurvData,distribution="coxph", shrinkage=0.001, n.trees=2000, interaction.depth=7, bag.fraction=1, train.fraction=1, verbose=F)
                                gbm.a <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=6, cv.folds=5, verbose=F, seed=913) #seed = my birthday :D !
			cat("done!\n");flush.console()

#===== 8. Training cox model using inducing minimal features ======= 
#
#  The minimalist features were selected by brute force to find the best features that maximize the concordance index in the training set
# based on the 484 overlapping samples between old and new training set using the old survival definition. 
#
#=====

                        cat("8.  Model B: Disease-specific Minimalist features...\n");flush.console()
				X = data.frame( cbind(meta["mitotic",], meta["mes.lymphneg",], lsxlymph, clinical$size, clinical$h.IDCnMED, gpr4) )
				colnames(X) = c("CIN", "MES_lymNumNeg", "LYMxlymNum", "size", "MED", "GPR4_g" )
				
                              cox.b = coxph(clinicalSurvData~., data=X)
			cat("done!\n");flush.console()

                              .self$model <- list(
						cm=cm, 
						cgbm=cgbm, 
						coxmodel=coxmodel,
						gbmmodel=gbmmodel, 
						knnmodel=knnmodel, 
						cox.a=cox.a,
						gbm.a=gbm.a,
						cox.b=cox.b
						)
                            },
#================================
#
# customPredict
#
#====================================

                            customPredict = function(exprData, copyData, clinicalFeaturesData)
                            {
                              rm(copyData)

				clnc = lazyImputeDFClncOslo(clinicalFeaturesData)
                              	clinical <- expandClncOslo(clnc)
				clinical$tr.CT = NULL; clinical$tr.HT = NULL
				
                        cat("Create metagene space...");flush.console()
                              meta = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)$metaSpace
                              #meta = CreateMetageneSpace(exprs(exprData), chosenProbes = .self$chosenProbes)
				#meta = t(sapply(.self$chosenProbes, function(pbs){apply(exprs(exprData)[pbs,], 2, mean)}))
				#meta = t( sapply(.self$attractome, function(pbs){ apply(exprs(exprData)[pbs,] , 2, function(a){mean(a, na.rm=TRUE)}) }  ) )
				ls = meta["ls",]
				meta = t(apply(meta, 1, function(x){ x - median(x)}))
			#===== Conditioning mesenchymal transition metagene by lymph node status
				idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
				mes.lymphneg = meta["mt",] * idx
				mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
				meta = rbind(meta, mes.lymphneg)
				
				idx = (clinical[,"lymph_nodes_positive"] > 3)
				ls.lymphpos = ls * idx
				ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
				meta = rbind(meta, ls.lymphpos)
				
				idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
				ls.erneg = ls * idx
				ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
				meta = rbind(meta, ls.erneg)

				lrrc48 = exprs(exprData)["ILMN_2353862",]
				lrrc48 = lrrc48 - median(lrrc48)

				lym.N = factor(clnc$lymph_nodes_positive < 1)
				lymph = clinical[,"lymph_nodes_positive"]
				lymph = sapply(lymph, function(x){min(x, 7)})
				lsxlymph = ls * (7 - lymph)

				gpr4 = exprs(exprData)["ILMN_2074477",]
			rm(exprData)

                        cat("done!\n");flush.console()

			p = matrix(NA,nrow=length(.self$model), ncol=ncol(meta))
#===== 1. Predict using clinical only model =====
			cat("1.  Predicting using clinical AIC model...");flush.console()
				X = clinical
                              p[1,] = predict(.self$model$cm, X)
			cat("done!\n");flush.console()
#===== 2. Predict using gbm model =====
			cat("2.  Predicting using clinical GBM model...");flush.console()
				X = X[, attr(.self$model$cm$term, "term.labels")]
                              best.iter=gbm.perf(.self$model$cgbm, method="cv", plot.it=FALSE)
			      cat("Best iter: ", best.iter, "\n", sep="");flush.console()
			      p[2,] = predict.gbm(.self$model$cgbm, X, best.iter)
			cat("done!\n");flush.console()
#===== 3. Predict using metagenes only model =====
			cat("3.  Predicting using metagenes AIC model...");flush.console()
                              X = data.frame(t(meta))
                              p[3,] = predict(.self$model$coxmodel, X)
                        cat("done!\n");flush.console()
#===== 4. Predict using GBM model =====
			cat("4.  Predicting using gbm model...");flush.console()
                              X = data.frame(t(meta))
                              best.iter=gbm.perf(.self$model$gbmmodel, method="cv", plot.it=FALSE)
			      cat("Best iter: ", best.iter, "\n", sep="");flush.console()
			      p[4,] = predict.gbm(.self$model$gbmmodel, X, best.iter)
			      #p3 = predict.gbm(.self$model[[2]], X, 1500)
			cat("done!\n");flush.console()
#===== 5. Predict using KNN model =====
			cat("5.  Predicting using KNN model ...");flush.console()
				knnmodel = .self$model$knnmodel
				qX = meta[names(knnmodel$x.train$concordance),]
				qC = clnc[,names(knnmodel$c.train$distWeight)]
				qC = t(preproClncKNN(qC, isFactorIn=knnmodel$c.train$isFactor, dwIn=knnmodel$c.train$distWeight)$clinical)
				wvec = c(abs(knnmodel$x.train$concordance-0.5), abs(knnmodel$c.train$concordance-0.5))
				qAll = rbind(qX, qC)
				trainDB = rbind(knnmodel$x.train$meta, t(knnmodel$c.train$clinical))
				trainTime = knnmodel$x.train$time
				out=ewknn.predict(trainDB, trainTime, qAll, wvec, k=floor(0.1*ncol(trainDB)))
				p[5,] = 365/out

			cat("done!\n");flush.console()
#===== 6. Predict using protective minimalist model =====
			cat("6.  Predicting using protective minimalist model ...");flush.console()
				X = data.frame( cbind(meta["mitotic",], meta["ls.erneg",], clinical$lymph_nodes_positive, meta["mes.lymphneg",], meta["susd3",], clinical$age_at_diagnosis) )
				colnames(X) = c("CIN", "LYM_ERNeg", "lymNum", "MES_lymNumNeg", "SUSD3", "age")
                              p[6,] = predict(.self$model$cox.a, X)
                        cat("done!\n");flush.console()
#===== 7. Predict using gbm model=====
			cat("7.  Predicting using gbm model...");flush.console()
                              best.iter=gbm.perf(.self$model$gbm.a, method="cv", plot.it=FALSE)
			      cat("Best iter: ", best.iter, "\n", sep="");flush.console()
			      p[7,] = predict.gbm(.self$model$gbm.a, X, best.iter)
			      #p2 = predict.gbm(.self$model[[2]], X, 2000)
			cat("done!\n");flush.console()
#===== 8. Predict using DS minimalist features=====
			cat("8. Predicting using disease specific Minimalist model ...");flush.console()
				X = data.frame( cbind(meta["mitotic",], meta["mes.lymphneg",], lsxlymph, clinical$size, clinical$h.IDCnMED, gpr4) )
				colnames(X) = c("CIN", "MES_lymNumNeg", "LYMxlymNum", "size", "MED", "GPR4_g" )
                              p[8,] = predict(.self$model[[8]], X)
			cat("done!\n");flush.console()

#===== Combining predictions =====
                              .self$predictions = p
				pout = matrix(NA, nrow=2, ncol=ncol(p))
				pout[1,] = apply(p[c(1:5, 8),], 2, mean)
				pout[2,] = apply(p[6:7,], 2, mean)
				pz = apply(pout, 1, scale)
                              p = apply(pz, 1, mean)
				#p = apply(pz, 1, mean)
                              return (p)
                            }
                            
                            )
                          )
