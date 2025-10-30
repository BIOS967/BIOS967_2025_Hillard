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
Timebin1=CoMatrixFix[,CoMatrixFix[1,]<=30000 & CoMatrixFix[1,]>=25000]
Timebin2=CoMatrixFix[,CoMatrixFix[1,]<25000 & CoMatrixFix[1,]>=20000]
Timebin3=CoMatrixFix[,CoMatrixFix[1,]<20000 & CoMatrixFix[1,]>=15000]
Timebin4=CoMatrixFix[,CoMatrixFix[1,]<15000 & CoMatrixFix[1,]>=10000]
Timebin5=CoMatrixFix[,CoMatrixFix[1,]<10000 & CoMatrixFix[1,]>=5000]
Timebin6=CoMatrixFix[,CoMatrixFix[1,]<5000 & CoMatrixFix[1,]>=0]

#Made list of the timebins for using in function
Bins=list(Timebin1,Timebin2,Timebin3,Timebin4,Timebin5,Timebin6)
#Used lapply function to remove the time row from each of the sites. Also turned all of the 0's and 1's to integers for the analysis
Binsnoyearint=lapply(Bins, FUN=function(x)x[-1,] %>% mutate_all(as.integer))

#Cooccur analysis for each time bin. WARNING this will take a while to run because of all the calculations made 
#mat argument is data (list of df in this case), threshold is to remove species that are too rare for an analysis
#type argument is setting the matrix as a site by species matrix. spp_names set to true as we use species names in the data
coocfun= lapply(Binsnoyearint,
                FUN=function(x)cooccur(mat=x, thresh=TRUE, type="spp_site",spp_names=TRUE))
CoocTime1=coocfun[[1]]
CoocTime2=coocfun[[2]]
CoocTime3=coocfun[[3]]
CoocTime4=coocfun[[4]]
CoocTime5=coocfun[[5]]
CoocTime6=coocfun[[6]]

#Putting the results from the co-occur analysis in list. Special object type called cooccur
Cooctotals=list(CoocTime1,CoocTime2,CoocTime3,CoocTime4,CoocTime5,CoocTime6)
#Calculating the proportion of significant pairs by dividing the sum of positive and negative pairs for each time bin by total pairs
TotalProportionsSig=lapply(Cooctotals, FUN=function(x)(x$co_occurrences)/(x$pairs))
#Calculating the proportion of aggregate pairs by dividing the number of aggregate pairs by the sum of significant species pairs in each time bin
TotalProportionsAgg=lapply(Cooctotals,FUN=function(x) (x$positive)/(x$co_occurrences))

#Turning the list of proportions into a dataframe for ggplot
TotalProportionsSig=as.data.frame(TotalProportionsSig)
#Adding a row to the dataframe for the categorical timebins for ggplot, added timebin labels
TotalProportionsSig[nrow(TotalProportionsSig)+1,]=c("Timebin1","Timebin2","Timebin3","Timebin4","Timebin5","Timebin6")
#Transposing the dataframe so that the values for proportions are a column not a row
TotalProportionsSig=data.frame(t(TotalProportionsSig))
#Turned the proportion values to numerics from characters
TotalProportionsSig$X1=as.numeric(TotalProportionsSig$X1)
#Plot for proportion of significant pairs. stat="identity for single values
#Used cowplot theme, added labels and title. Changed element text to make it look better
ggplot(TotalProportionsSig, aes(x=X2,y=X1))+geom_bar(stat="identity", fill="skyblue", width=0.5)+
  theme_cowplot()+labs(x="Time",y="Proportion of Nonrandom Pairs",title="Proportion of Nonrandom Pairs Through Time in North America")+
  theme(plot.title=element_text(size=12), axis.title.y=element_text(size=10),
        axis.title.x=element_blank(), axis.text.y=element_text(size=10))

#Turning the list of proportions into a dataframe for ggplot
TotalProportionsAgg=as.data.frame(TotalProportionsAgg)
#Adding a row to the dataframe for the categorical timebins for ggplot, added timebin labels
TotalProportionsAgg[nrow(TotalProportionsAgg)+1,]=c("Timebin1","Timebin2","Timebin3","Timebin4","Timebin5","Timebin6")
#Transposing the dataframe so that the values for proportions are a column not a row
TotalProportionsAgg=data.frame(t(TotalProportionsAgg))
#Turned the proportion values to numerics from characters
TotalProportionsAgg$X1=as.numeric(TotalProportionsAgg$X1)
#Plot for proportion of aggregated pairs. stat="identity for single values
#Used cowplot theme, added labels and title. Changed element text to make it look better, also changed axis range to illuminate trend
ggplot(TotalProportionsAgg, aes(x=X2,y=X1))+geom_bar(stat="identity", fill="tomato", width=0.5)+
  theme_cowplot()+labs(x="Time",y="Proportion of Aggregated Pairs", title="Proportion of Aggregated Pairs Through Time in North America")+
  theme(plot.title=element_text(size=12), axis.title.y=element_text(size=10),
          axis.title.x=element_blank(), axis.text.y=element_text(size=10))+
  coord_cartesian(ylim=c(0.4,0.9))
  
