### Title:    Data generation
### Author:   Rosember Guerra-Urzola
### Created:  16-07-2021
### Modified: 

# install.packages("doParallel")    # Install doParallel package
# install.packages("MASS")          # Install MASS package
# install.packages("mvtnorm")       # Install mvtnorm package

rm(list = ls(all.names = TRUE))

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

# setting the number of cores 
library(doParallel)

# sparse PCA data simulation #

set.seed(2019)

# sparse PCA data simulation #
error = c(.05, .2,.8)           # Proportion of explained variance
ps = c(.2,.8)         # Proportion of sparsity
k = 3        # Number of components
i = c(100)           # Sample size
j = c(50,100,500)  # Number of variables
replications = c(100)         # Number of repetitions


design_matrix <- expand.grid( J=j,I=i,PS=ps,K=k,ERROR =error)
rep_desig = rep(1:nrow(design_matrix), times = replications)
design_matrix_replication <- design_matrix[rep_desig, ]
Index_design = matrix(rep(0,replications*nrow(design_matrix)),ncol = replications)
for (r in 1:nrow(design_matrix)) {
  Index_design[r,]= which(rep_desig == r)
}
SimDesig = list(index = Index_design, Matrix = design_matrix,Ndatasets = nrow(design_matrix_replication))
save(SimDesig, file=paste0("../data/SimDesig.RData"))
rm(list=setdiff(ls(), "design_matrix_replication"))

variance_comp = c(31, 30, 29)

library(doParallel)
no_cores <- detectCores() 
c1 <- makePSOCKcluster(floor( no_cores*.3))
registerDoParallel(c1)


results_sim1_data1 <- foreach(i=1:nrow(design_matrix_replication),
                              .options.RNG = 2018,
                              .packages = c("MASS","gtools"),
                              .combine=rbind)%dopar%{
                                source('makeData.R',local = TRUE)
                                #i=2
                                J = design_matrix_replication$J[i]
                                s_sample =design_matrix_replication$I[i]
                                K = design_matrix_replication$K[i]
                                PS = design_matrix_replication$PS[i]
                                ERROR = design_matrix_replication$ERROR[i]
                                
                                support <- matrix(1, J, K)
                                support <- sparsify(support, PS)
                                variances <- makeVariance(variance_comp, p = J, error = ERROR)
                                datObject <- makeDat(n = s_sample, p = J, ncomp = K, support, 
                                                     variances = variances)
                                # Saving data
                                out = list(data = datObject$X, P = datObject$P, Sigma = datObject$Sigma,
                                           I=s_sample,J=J, K = K, PS = PS, ERROR=ERROR,supp=support, teller = i)
                                save(out, file=paste0("../data/SimData",i, ".RData"))
                                
                              }

# stop Cluster
stopCluster(c1)
