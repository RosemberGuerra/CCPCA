### Title:    Analysing the results from: CCPCA and SPCA
### Author:   Rosember Guerra-Urzola
### Created:  19-07-2021
### Modified: 



rm(list = ls(all.names = TRUE))

library(data.table)
library(ggplot2)
library(reshape)
# Setting the working directory to the source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

load('../data/SimDesig.RData')
load(paste0("../output/estimation",1,".RData"))
Joint_dataset = output$measures

for (i in 2:SimDesig$Ndatasets) {
  load(paste0("../output/estimation",i,".RData"))
  Joint_dataset = rbind(Joint_dataset,output$measures)
  # View(Joint_dataset)
}
Joint_dataset[,c("Method","I","J","K","PS","Noise")]= lapply(Joint_dataset[,c("Method","I","J","K","PS","Noise")],
                                                  as.factor)
Joint_dataset[,c("Tucker",'Rec',"Zero_Rec","Total_Rec",'PEV')]= lapply(Joint_dataset[,c("Tucker",'Rec',"Zero_Rec","Total_Rec",'PEV')],
                                                             function(x) as.numeric(as.character(x)))

levels(Joint_dataset$I) = c('I = 100')
levels(Joint_dataset$PS) = c('PS = 20%','PS = 80%')
levels(Joint_dataset$Noise) = c('Noise = 5%','Noise = 20%','Noise = 80%')
# Joint_dataset$J = relevel(Joint_dataset$J2,ref = "50")

str(Joint_dataset)
Joint_dataset= setDT(Joint_dataset)
class(Joint_dataset)

### Tucker Congruence ###
# load('SimulationResults.RData')

ggplot(data = Joint_dataset, aes(x = Method, y = Tucker)) +
  geom_boxplot(aes(fill = J), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(PS), cols = vars(Noise)) +
  geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
  # geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="Variables")) +
  xlab("") + 
  # ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Tucker Congruence") + xlab('Method')+
  # coord_fixed(ratio = 2/1) +
  theme_bw()
  
Width = 10.24
Height = 8.85
ggsave("Tucker.eps",path = "../../main/plots",
       width = Width, height = Height ,limitsize = FALSE)

ggplot(data = Joint_dataset, aes(x = Method, y = Total_Rec)) + 
  geom_boxplot(aes(fill = J), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(PS), cols = vars(Noise)) +
  geom_hline(yintercept = 0.6, col = "black", linetype = 2) +
  # geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="Variables")) +
  xlab("") + 
  # ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Total Recovery") + xlab('Method')+
  # coord_fixed(ratio = 4/1) +
  theme_bw()
ggsave("RecoveryWeights.eps",path = "../../main/plots",
       width = Width, height = Height ,limitsize = FALSE)

ggplot(data = Joint_dataset, aes(x = Method, y = PEV)) + 
  geom_boxplot(aes(fill = J), outlier.size=0.3, lwd=0.2) +
  facet_grid(rows = vars(PS), cols = vars(Noise)) +
  #geom_hline(yintercept = 0.6, col = "black", linetype = 2) +
  # geom_vline(xintercept = 3.5, col = "black", linetype = 2) +
  guides(fill=guide_legend(title="Variables")) +
  xlab("") + 
  # ylim(0,1)+
  scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
  ylab("Total Recovery") + xlab('Method')+
  # coord_fixed(ratio = 4/1) +
  theme_bw()


# Joint_dataset2 = melt(Joint_dataset,id.vars = c('Method','I','J','K','PS','Noise'),
#                       measure.vars=c('Rec','Zero_Rec','Total_Rec'))
# 

### second part ###

meanMAB = matrix(rep(0,nrow(SimDesig$index)*2),ncol = 2)
meanVAR = matrix(rep(0,nrow(SimDesig$index)*2),ncol = 2)
meanMSE = matrix(rep(0,nrow(SimDesig$index)*2),ncol = 2)
R = ncol(SimDesig$index) # number of repetitions 
for (i in 1:nrow(SimDesig$index) ) {
  count =0
  K = SimDesig$Matrix[i,'K']
  J = SimDesig$Matrix[i,'J']
  W_aes_ccpca = matrix(rep(0,J*K),nrow = J)
  W_aes_spca = matrix(rep(0,J*K),nrow = J)
  for(r in SimDesig$index[i,]){
    count =count+1
    load(paste0("../output/estimation",r,".RData"))
    W_aes_ccpca = (output$W_ccpca + (count-1)*W_aes_ccpca)/count
    W_aes_spca = (output$W_spca + (count-1)*W_aes_spca)/count
    }
  # aesw_ccpca = c(aesw_ccpca,assign( paste0('Waes',i),W_aes_ccpca))
  # aesw_spca = c(aesw_spca, paste0('Waes',i)=W_aes_spca)
  for (r in SimDesig$index[i,]) {
    load(paste0("../output/estimation",r,".RData"))
    ## mean absolute bias ##
    meanMAB[i,1] = meanMAB[i,1] + sum(abs(W_aes_ccpca- output$W_ccpca)) 
    meanMAB[i,2] = meanMAB[i,2] + sum(abs(W_aes_spca- output$W_spca)) 
    
    ## mean VAR ##
    meanVAR[i,1] = meanVAR[i,1]+ sum((W_aes_ccpca- output$W_ccpca)^2)
    meanVAR[i,2] = meanVAR[i,2]+ sum((W_aes_spca- output$W_spca)^2)
    
    ## mean MSE ##
    meanMSE[i,1] = meanMSE[i,1]+ sum((output$W- output$W_ccpca)^2)
    meanMSE[i,2] = meanMSE[i,2]+ sum((output$W- output$W_spca)^2)
  }# remember to divide by R
  meanMAB[i,] = meanMAB[i,]/(J*K)
  meanVAR[i,] = meanVAR[i,]/(J*K)
  meanMSE[i,] = meanMSE[i,]/(J*K)
  
}

meanMAB = meanMAB/R
colnames(meanMAB) = c('MAB-CCPCA','MAB-SPCA')
meanVAR = meanVAR/R
colnames(meanVAR) = c('VAR-CCPCA','VAR-SPCA')
meanMSE = meanMSE/R
colnames(meanMSE) = c('MSE-CCPCA','MSE-SPCA')
results2 = cbind(SimDesig$Matrix,meanMAB,meanVAR,meanMSE)

write.table(results2,file = 'Bias-variance.csv',sep = ',',col.names = TRUE)

View(results2)
save.image(file = 'SimulationResults.RData')

