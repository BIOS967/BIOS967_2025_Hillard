#Reading data csv
CoMatrix=read.csv("data/Complete PAtable_new.csv")
#Loading needed packages
#I had to manually install the cooccur package that I need as it wasn't on CRAN
#I was able to download it manually as a .tar file and install it manually
install.packages("gmp")
install.packages("reshape2")
library(ggplot2)
library(cowplot)
library(tidyverse)
library(gmp)
library(reshape2)
library(cooccur)


#Swapped columns and rows to match the structure of co-occurence matrices standard in literature
#Source=https://stackoverflow.com/questions/33643181/how-do-i-flip-rows-and-columns-in-r
CoMatrixFix=data.frame(t(CoMatrix[-1]))
colnames(CoMatrixFix)= CoMatrix[,1]

#Filter different columns based on time period which is the first row
#Source: https://stackoverflow.com/questions/49396354/filter-r-data-frame-by-columns-instead-of-by-rows
Pleistocene=CoMatrixFix[,CoMatrixFix[1,]>=12000]


Comatrixnoyear=CoMatrixFix[-1,]
Comatrixnoyearint= Comatrixnoyear %>% mutate_all(as.integer)
bruhswouse=cooccur(Comatrixnoyearint, type="spp_site",thresh=TRUE,spp_names=TRUE)
summary(bruhswouse6)

Timebin1=CoMatrixFix[,CoMatrixFix[1,]<=30000 & CoMatrixFix[1,]>=25000]
Timebin2=CoMatrixFix[,CoMatrixFix[1,]<25000 & CoMatrixFix[1,]>=20000]
Timebin3=CoMatrixFix[,CoMatrixFix[1,]<20000 & CoMatrixFix[1,]>=15000]
Timebin4=CoMatrixFix[,CoMatrixFix[1,]<15000 & CoMatrixFix[1,]>=10000]
Timebin5=CoMatrixFix[,CoMatrixFix[1,]<10000 & CoMatrixFix[1,]>=5000]
Timebin6=CoMatrixFix[,CoMatrixFix[1,]<5000 & CoMatrixFix[1,]>=0]
Timebin1noyear=Timebin1[-1,]
Timebin1noyearint=Timebin1noyear %>% mutate_all(as.integer)
bruhswouse1=cooccur(Timebin1noyearint, type="spp_site",thresh=TRUE, spp_names=TRUE)

coocfun= lapply(list())
Timebin1ProportionSig=(bruhswouse1$positive+bruhswouse1$negative)/(bruhswouse1$positive+bruhswouse1$negative+bruhswouse1$random)
Timebin1ProportionAgg=(bruhswouse1$positive)/(bruhswouse1$positive+bruhswouse1$negative)
