# Download and prepocessing data set: Authism 
# Rosember Guerra-Urzola 
# created: 21-06-2021
# edited: 15-07-2021

# Note: Before running this script please be aware that this application was conducted
#       Using a a super computer with the following specificaitons.
#       OS version: Windows 2012 R2
#       CPU: 24 Core 6.60 GHz Intel Xeon(R) Gold 6126 (Hyper-Threaded)
#       Memory: 523779 MB


# Using this file, we can download the data set from the
# url: 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz'
# name of data set: GSE7329_series_matrix.txt.gz
# The data is downloaded in as .gz format
# The data 'GSE7329_series_matrix.txt.gz' is save in the same directory as this file.

rm(list = ls(all.names = TRUE))

## pre-processing data set ##

# Setting the working directory to the source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

###  Download data set ### Do this only one time !!!
# Data_loacation= 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz'
# destfile= paste0(getwd(),"/GSE7329_series_matrix.txt.gz")
# download.file(Data_loacation,destfile)

full_data=gzfile('GSE7329_series_matrix.txt.gz','rt')  
full_data=read.csv(full_data,header=F,sep = '\t')

### Column names from the data ###
col_names= c()
for (i in 736:751) {
  # col_names=c(col_names, as.character(full_data[i,1]),as.character(dat[i,2]))
  col_names=c(col_names,sapply(full_data[i,], as.character))
} 
#View(col_names)

### Data neede for the application ###
matrix1= full_data[752:703647,]
matrix2= matrix(rep(0, prod(dim(matrix1))), ncol = dim(matrix1)[2])
for(i in 1:nrow(matrix1) ){
  matrix2[i,]= as.numeric(sapply(matrix1[i,],as.character))
}

Data_Autism=matrix(data = t(matrix2),byrow = TRUE, ncol = length(col_names))
colnames(Data_Autism)=col_names # adding names
IndexCol= Data_Autism[,1]
Data_Autism=Data_Autism[,-c(1,dim(Data_Autism)[2])] # removing first and last column 
Data_Autism=Data_Autism[,-c(12,13,27)] # Removing mislabel individuals 
Data_Autism= t(Data_Autism)



Data_Autism = Data_Autism[,-which(colSums(Data_Autism)==0)] # This is needed to remove the NaN
dim(Data_Autism)
Data_Autism=scale(Data_Autism, center = TRUE, scale = TRUE)
sum(is.na(Data_Autism))
Data_Autism=as.data.frame(Data_Autism) # Before with transpose
#write.csv(Data_Autism, file = 'Data_Autism.csv',sep = ',')

### --- estimation of the index of sparseness ---###

library(elasticnet)
library(MASS)
library(ggplot2)
library(doParallel)

## PCA estimation using full data set ##

# estimation for pca values#

X = as.matrix(Data_Autism)

rm(list=setdiff(ls(), "X"))
K= 3 # number of components

Xsvd= svd(X, nu = K, nv = K)
totalvariance = sum(Xsvd$d^2)
pev_pca = sum(Xsvd$d[1:K]^2)/totalvariance
alpha = Xsvd$d[1]^2
CorrM = t(X)%*%X
## plot pca ##
scores_pca = Xsvd$u%*%diag(Xsvd$d[1:K])
scores_pca = as.data.frame(scores_pca)
colnames(scores_pca)= c('PC1','PC2','PC3')

subject = c(rep('dup15',7),rep('FMR1',6),rep('Control',14))
scores_pca$subject =subject

write.table(scores_pca, file = "scores_pca.txt", sep = ",",row.names = TRUE,col.names = TRUE)


library(ggalt)
library(ggforce)
ggplot(scores_pca,aes(x=PC2,y=PC3,shape=subject))+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=subject,linetype =subject))+
  scale_fill_grey(
    start = 0.2,
    end = 0.8,
    na.value = "red",
    aesthetics = "fill"
  )+
  geom_point(size=3)+
  labs(y = "PC2", x = "PC3",title = 'PCA using full data set')+
  theme_bw()

## IS estimation using CCPCA ##

# Registering the cluster

J= dim(X)[2]
cardinality= floor(seq(from=100,to=(J-100),length.out=100))

results= matrix(c(0,0,pev_pca,J),nrow = 1)

source('ccpca_functions.R')

pb = txtProgressBar(min = 0, max = length(cardinality), initial = 0)
count = 0
for (i in cardinality) {
  PS = 1-i/J
  Comp=CCPCA(X,K,rho = i,alpha=alpha,Corr=CorrM,W=Xsvd$v, comp=TRUE,
             maxiter= 1000, tolerance = 1e-3,progress = FALSE)
  S = X%*%Comp$W
  Rccpca <- qr.R(qr(S))
  pev_ccpca <- sum(diag(Rccpca^2))/totalvariance
  IS_ccpca = pev_ccpca*pev_pca*PS
  results = rbind(results,t(c(PS,IS_ccpca,pev_ccpca,i)))
  count =count+1
  setTxtProgressBar(pb,count)
}
results = as.data.frame(results)
colnames(results)= c('PS','IS','PEV','card')
View(results)

write.table(results, file = "fulldataseIS-PEV-PS-card.txt", sep = ",",row.names = FALSE,col.names = TRUE)


Ind= which.max(results[,'IS'])

Comp=CCPCA(X,K,rho = cardinality[Ind],alpha=alpha,Corr=CorrM,W=Xsvd$v, comp=TRUE,
           maxiter= 1000, tolerance = 1e-3, progress=TRUE)

scores_ccpca = X%*%Comp$W
Rccpca <- qr.R(qr(scores_ccpca))
pev_ccpca <- sum(diag(Rccpca^2))/totalvariance
IS_ccpca = pev_ccpca*pev_pca*PS

scores_ccpca = data.frame(scores_ccpca)
colnames(scores_ccpca)= c('PC1','PC2','PC3')
label = c(rep('dup15',7),rep('FMR1',6),rep('Control',14))
scores_ccpca$subject =label

write.table(scores_ccpca, file = "scores_ccpca.txt", sep = ",",row.names = TRUE,col.names = TRUE)

# unique variables #
card_W =  colSums( Comp$W != 0)
uniq_var = sum( colSums(t(Comp$W)) !=0)

###-- spca implementation --### 
spca_estimation =  spca(x=X,K = K, sparse = "varnum", para =rep(cardinality[Ind],K) ,type = "predictor")

W_scpa = spca_estimation$loadings
Z_spca = X%*%W_scpa
pev_spca = sum(spca_estimation$pev)

###

pev_summary = as.data.frame(matrix(c(pev_pca,pev_ccpca,pev_spca,
                                     0,results[Ind,'PS'],results[Ind,'PS'],
                                     J,(cardinality[Ind]),(cardinality[Ind])),
                                   nrow = 3,ncol = 3,byrow = TRUE))
colnames(pev_summary)=c('PCA','CCPCA','SPCA')
rownames(pev_summary)=c('pev','PS','card')
View(pev_summary)
write.table(pev_summary, file = "pevsummary.txt", sep = ",",
            row.names = TRUE,col.names = TRUE)


## summary scores ##
scores_summary = cbind(scores_pca[,-4],scores_ccpca[,-4],Z_spca)
colnames(scores_summary) =c('PCA1','PCA3','PCA2',
                           'CCPCA1','CCPCA2','CCPCA3',
                           'SPCA1','SPCA2','SPCA3')
scores_summary$subject = label
View(scores_summary)
write.table(scores_summary, file = "ScoresSummary.txt", sep = ",",
            row.names = TRUE,col.names = TRUE)

## plots ##
# load data if needed ##
library(ggplot2)

## plot IS-PEV ##

p <- ggplot(results, aes(x = card))
p <- p + geom_line(aes(y = IS, linetype = "IS"))
p <- p + geom_line(aes(y = PEV/4, linetype = "PEV"))
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*3, name = "PEV"))
p <- p + scale_linetype_manual(values=c("solid", "dashed"))
# p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "IS",
              x = "PS", linetype = " ")#,title = "Varimax")#+ ggtitle("Simplimax") 
# p <- p + theme(legend.position = c(0.8, 0.9))
p +  theme_bw()+ theme(plot.title = element_text(size=20 ,hjust = 0.5),
                       legend.position = c(0.2, 0.8),
                       legend.text = element_text(size = 16),
                       axis.text=element_text(size=16),
                       axis.title = element_text(size = 16))
# ggsave("Varimax.eps",width = Width, height = Height,
#         limitsize = FALSE)


## ploting the component scores #

library(ggalt)
library(ggforce)

ggplot(scores_ccpca,aes(x=PC3,y=PC2,shape=subject))+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=subject,linetype =subject))+
  scale_fill_grey(
    start = 0.2,
    end = 0.8,
    na.value = "red",
    aesthetics = "fill"
  )+
  geom_point(size=3)+
  labs(y = "Score3", x = "Score2",title = 'C. Scores - CCPCA')+
  theme_bw()

ggplot(scores_summary,aes(x=CCPCA3,y=CCPCA2,shape=subject))+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=subject,linetype =subject))+
  scale_fill_grey(
    start = 0.2,
    end = 0.8,
    na.value = "red",
    aesthetics = "fill"
  )+
  geom_point(size=3)+
  labs(y = "Score3", x = "Score2",title = 'C. Scores - SPCA')+
  theme_bw()

ggplot(scores_summary,aes(x=PCA2,y=PCA3,shape=subject))+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=subject,linetype =subject))+
  scale_fill_grey(
    start = 0.2,
    end = 0.8,
    na.value = "red",
    aesthetics = "fill"
  )+
  geom_point(size=3)+
  labs(y = "Score3", x = "Score2",title = 'C. Scores - PCA')+
  theme_bw()


save.image(file ='finalApplicaiton.RData')
load('finalApplicaiton.RData')
