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

Timebin1=CoMatrixFix[,CoMatrixFix[1,]>=30000]
Bruh swouse bruh swouse

