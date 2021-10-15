### Title:    Estimation of CCPCA and SPCA
### Author:   Rosember Guerra-Urzola
### Created:  16-07-2021
### Modified: 

rm(list = ls(all.names = TRUE))

library(doParallel)

# Setting the working directory to the source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# Registering the cluster
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# library(elasticnet)
# library(combinat)
# library(MASS)
load('../data/SimDesig.RData')

estimation = foreach(i=1:SimDesig$Ndatasets,.packages = c("elasticnet","combinat","MASS"),
                     .combine=rbind)%dopar%{
                       source('ccpca_functions.R',local = TRUE)
                       source('tucker.R',local = TRUE )
                       
                       # setwd(file.path(DataDir))
                       # i=2
                       load(paste0("../data/SimData",i,".RData")) # generated data
                       
                       ## initial values ##
                       Xsvd = svd(out$data)
                       alpha = Xsvd$d[1]^2
                       rho = colSums(out$P[,1:out$K] !=0)
                       TotalVariance = sum(Xsvd$d^2)
                       W0= Xsvd$v[,1:out$K]
                
                       ## estimation ##
                       ccpca_estimation = CCPCA(X = out$data,K=out$K,alpha = alpha,rho = rho,
                                                W=W0)
                       spca_estimation = spca(x=out$data ,K = out$K, sparse = "varnum", 
                                              para =rho ,type = "predictor")
                       
                       
                       
                       ## measures ##
                       tucker_ccpca = tuckerCongruence(ccpca_estimation$W,out$P[,1:out$K])
                       tucker_spca = tuckerCongruence(spca_estimation$loadings,out$P[,1:out$K])
                       
                       # support_P = out$P !=0
                       # dim(support_P)
                       corr_class_ccpca=  corClass(out$P[,1:out$K], ccpca_estimation$W, comdisspar = out$supp)
                       corr_class_spca=  corClass(out$P[,1:out$K], spca_estimation$loadings, comdisspar = out$supp)
                       # 
                       recovery_ccpca = sum(rho*corr_class_ccpca$nonZeroCoefFound)/sum(rho)
                       recovery_spca = sum(rho*corr_class_spca$nonZeroCoefFound)/sum(rho)
                       # 
                       zero_recovery_ccpca= sum((out$J-rho)*corr_class_ccpca$zeroCoefFound)/sum((out$J-rho))
                       zero_recovery_spca= sum((out$J-rho)*corr_class_spca$zeroCoefFound)/sum((out$J-rho))
                       # 
                       # total_recovery_ccpca = (sum(rho)*recovery_ccpca+sum(out$J-rho)*zero_recovery_ccpca)/(out$J*out$K)
                       # total_recovery_spca = (sum(rho)*recovery_spca+sum(out$J-rho)*zero_recovery_spca)/(out$J*out$K)
                       total_recovery_ccpca = TotalRecovery(W=ccpca_estimation$W,supp = out$supp)
                       total_recovery_spca = TotalRecovery(W=spca_estimation$loadings, supp = out$supp)
                       
                       pev_ccpca= sum(pev(W= ccpca_estimation$W,X=out$data,totalvariance = TotalVariance))
                       pev_spca= sum(pev(W= spca_estimation$loadings,X=out$data,totalvariance = TotalVariance))
                       
                       
                       measure_summary = matrix(c('CCPCA',out$I,out$J,out$K,out$PS,out$ERROR,tucker_ccpca,recovery_ccpca,zero_recovery_ccpca,total_recovery_ccpca,
                                                  pev_ccpca,
                                                  'SPCA',out$I,out$J,out$K,out$PS,out$ERROR,tucker_spca,recovery_spca,zero_recovery_spca,total_recovery_spca,
                                                  pev_spca),
                                                nrow = 2,byrow = TRUE)
                       measure_summary = as.data.frame(measure_summary)
                       colnames(measure_summary ) = c('Method','I','J','K','PS','Noise','Tucker','Rec','Zero_Rec','Total_Rec','PEV')
                       # View(measure_summary)
                       output = list(W = out$P[,1:out$K],W_ccpca=ccpca_estimation$W,
                                     W_spca = spca_estimation$loadings, measures = measure_summary)
                       
                       save(output, file = paste0("../output/estimation",i,".RData"))
                     }
                       
# stop Cluster
stopCluster(c1)
