require(rms)
require(MASS)
require(survival)
require(predictiveModeling)
require(gbm)
require(caret)
require(devtools)

#  download.file("http://dl.dropbox.com/u/11986954/Oslo/DreamBox7_0.214.tar.gz", destfile="./DreamBox7_0.214.tar.gz")
#  install.packages("./DreamBox7_0.214.tar.gz", repos=NULL)

  install_github(repo="DreamBox7", username="weiyi-bitw")
  library(DreamBox7)

setRefClass(Class = "PredictiveModel")

#' GoldiModel
#'
#' Modified from DemoClinicalOnlyModel from BCC challenge.
#'
#' @author Wei-Yi Cheng 
#' @Revise Tai-Hsien Ou Yang
#' @export

GoldiloxModel <- setRefClass(Class  = "GoldiloxModel",
                          contains = "PredictiveModel",
                          fields   = c("model", "attractome", "annot", "predictions","mdns", "chosenProbes_g","chosenProbes", "dssurv", "w"),
                          methods  = list(
                            
                            initialize = function(...){
  				download.file("http://dl.dropbox.com/u/11986954/Oslo/metabric.surv.ds.rdata", destfile="./metabric.surv.ds.rdata")
				load("metabric.surv.ds.rdata")
				.self$dssurv = surv.ds
                              return(.self)
                            },
                            
                            customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,clinicalSurvData_DS,...)
                            {
                              if(class(clinicalSurvData) != "Surv"){
                                stop("Expecting 'responseData' object of type 'Surv'")
                              }
				source("http://dl.dropbox.com/u/11986954/Oslo/select/syn1417992_fix.R")
							
   				gdModel.os <- GoldiModel$new()
				gdModel.os$customTrain(exprData, copyData, clinicalFeaturesData,clinicalSurvData) 
                                p1 <- gdModel.os$customPredict(exprData, copyData, clinicalFeaturesData)

				gdModel.ds <- GoldiModel$new()
				gdModel.ds$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData_DS)
				p2 = gdModel.ds$customPredict(exprData, copyData, clinicalFeaturesData)

                        rm(copyData)
			rm(exprData)

				pp = rbind(p1, p2)
				weights = BFFW(pp, clinicalSurvData, w = rep(1, 2))			
				.self$w = weights

                              .self$model <- list(
						osmodel = gdModel.os,
						dsmodel = gdModel.ds
						)
                            },
#================================
#
# customPredict
#
#====================================

                            customPredict = function(exprData, copyData, clinicalFeaturesData)
                            {
				#meta.cnv = CreateMetageneSpace(exprs(copyData), .self$cnvome, .self$annot.cnv)$metaSpace
				

				p1 = .self$model$osmodel$customPredict(exprData, copyData, clinicalFeaturesData)
				p2 = .self$model$dsmodel$customPredict(exprData, copyData, clinicalFeaturesData)

				rm(exprData)
                                rm(copyData)
				
                              .self$predictions = rbind(p1, p2)
                              p = p1 + .self$w[2] * p2
				#p = apply(pz, 1, mean)
                              return (p)
                            }
                            
                            )
                          )





#source("retrain.R")
#demoPredictiveModel <-GoldiModel$new()
#demoPredictiveModel$customTrain(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
#trainPredictions <- demoPredictiveModel$customPredict(trainingData$exprData, trainingData$copyData, trainingData$clinicalFeaturesData)
#tr.ci=getCCDIdx(trainPredictions,trainingData$clinicalSurvData)
#print(tr.ci)

#submitCompetitionModel(modelName="Attractor Metagenes Model 101507 OSDS",trainedModel=demoPredictiveModel,rFiles=c("retrain.R","trainingData.clinicalSurvData.ds.rdata"),isPracticeModel=FALSE) 







