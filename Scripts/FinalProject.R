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
save(Cooctotals, file="data/Cooctotals.rdata")
load("data/Cooctotals.rdata")

#Calculating the proportion of significant pairs by dividing the sum of positive and negative pairs for each time bin by total pairs
TotalProportionsSig=lapply(Cooctotals, FUN=function(x)(x$co_occurrences)/(x$pairs))
#Calculating the proportion of aggregate pairs by dividing the number of aggregate pairs by the sum of significant species pairs in each time bin
TotalProportionsAgg=lapply(Cooctotals,FUN=function(x) (x$positive)/(x$co_occurrences))

#Turning the list of proportions into a dataframe for ggplot
TotalProportionsSig=as.data.frame(TotalProportionsSig)
#Adding a row to the dataframe for the categorical timebins for ggplot, added timebin labels
TotalProportionsSig[nrow(TotalProportionsSig)+1,]=c("30000","25000","20000","15000","10000","5000")
#Transposing the dataframe so that the values for proportions are a column not a row
TotalProportionsSig=data.frame(t(TotalProportionsSig))
#Turned the proportion values to numerics from characters
TotalProportionsSig$X1=as.numeric(TotalProportionsSig$X1)
#Plot for proportion of significant pairs. stat="identity for single values
#Used cowplot theme, added labels and title. Changed element text to make it look better
ggplot(TotalProportionsSig, aes(x=X2,y=X1))+geom_bar(stat="identity", fill="skyblue", width=0.5, color="black",linewidth=1)+
  theme_cowplot()+labs(x="Time",y="Proportion of Nonrandom Pairs",title="Proportion of Nonrandom Pairs Through Time in North America")+
  theme(plot.title=element_text(size=11.5), axis.title.y=element_text(size=10),
        axis.title.x=element_blank(), axis.text.y=element_text(size=10))
  

#Turning the list of proportions into a dataframe for ggplot
TotalProportionsAgg=as.data.frame(TotalProportionsAgg)
#Adding a row to the dataframe for the categorical timebins for ggplot, added timebin labels
TotalProportionsAgg[nrow(TotalProportionsAgg)+1,]=c("30000","25000","20000","15000","10000","5000")
#Transposing the dataframe so that the values for proportions are a column not a row
TotalProportionsAgg=data.frame(t(TotalProportionsAgg))
#Turned the proportion values to numerics from characters
TotalProportionsAgg$X1=as.numeric(TotalProportionsAgg$X1)
#Plot for proportion of aggregated pairs. stat="identity for single values
#Used cowplot theme, added labels and title. Changed element text to make it look better, also changed axis range to illuminate trend
ggplot(TotalProportionsAgg, aes(x=X2,y=X1))+geom_bar(stat="identity", fill="tomato", width=0.5, color="black",linewidth=1)+
  theme_cowplot()+labs(x="Time",y="Proportion of Aggregated Pairs", title="Proportion of Aggregated Pairs Through Time in North America")+
  theme(plot.title=element_text(size=11.5), axis.title.y=element_text(size=10),
          axis.title.x=element_blank(), axis.text.y=element_text(size=10))+
  coord_cartesian(ylim=c(0.4,0.9))
  
SpeciesDiets=read.csv("data/TraitdataTestVersion.csv")

probtable1=prob.table(Cooctotals[[1]])
probtable1=probtable1%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable1$Pairtype[probtable1$p_gt<=0.05] = "Positive"
probtable1$Pairtype[probtable1$p_lt<=0.05]="Negative"
probtable1$Pairtype[probtable1$p_gt>0.05 & probtable1$p_lt>0.05]="Random"
probtable1$time=1

probtable2=prob.table(Cooctotals[[2]])
probtable2=probtable2%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable2$Pairtype[probtable2$p_gt<=0.05] = "Positive"
probtable2$Pairtype[probtable2$p_lt<=0.05]="Negative"
probtable2$Pairtype[probtable2$p_gt>0.05 & probtable2$p_lt>0.05]="Random"
probtable2$time=2

probtable3=prob.table(Cooctotals[[3]])
probtable3=probtable3%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable3$Pairtype[probtable3$p_gt<=0.05] = "Positive"
probtable3$Pairtype[probtable3$p_lt<=0.05]="Negative"
probtable3$Pairtype[probtable3$p_gt>0.05 & probtable3$p_lt>0.05]="Random"
probtable3$time=3

probtable4=prob.table(Cooctotals[[4]])
probtable4=probtable4%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable4$Pairtype[probtable4$p_gt<=0.05] = "Positive"
probtable4$Pairtype[probtable4$p_lt<=0.05]="Negative"
probtable4$Pairtype[probtable4$p_gt>0.05 & probtable4$p_lt>0.05]="Random"
probtable4$time =4

probtable5=prob.table(Cooctotals[[5]])
probtable5=probtable5%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable5$Pairtype[probtable5$p_gt<=0.05] = "Positive"
probtable5$Pairtype[probtable5$p_lt<=0.05]="Negative"
probtable5$Pairtype[probtable5$p_gt>0.05 & probtable5$p_lt>0.05]="Random"
probtable5$time=5

probtable6=prob.table(Cooctotals[[6]])
probtable6=probtable6%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable6$Pairtype[probtable6$p_gt<=0.05] = "Positive"
probtable6$Pairtype[probtable6$p_lt<=0.05]="Negative"
probtable6$Pairtype[probtable6$p_gt>0.05 & probtable6$p_lt>0.05]="Random"
probtable6$time=6

Probtablecomb= list(probtable1,probtable2,probtable3,probtable4,probtable5,probtable6)
Probtablecomb= lapply(Probtablecomb,FUN=function(x)(x[x$Pairtype != "Random",]))
Splitprobtable= lapply(Probtablecomb, FUN=function(x)(split(x, x$Pairtype)))
Posprobtable= lapply(Splitprobtable, `[[`, "Positive")
Posprobtable=bind_rows(Posprobtable, .id="source")
Negprobtable= lapply(Splitprobtable, `[[`, "Negative")
Negprobtable=bind_rows(Negprobtable, .id="source")


Posprobtable$sp1_diet = as.character(Posprobtable$sp1_diet)
Posprobtable$sp2_diet = as.character(Posprobtable$sp2_diet)
Posprobtable$diet_pair = ifelse(is.na(Posprobtable$sp1_diet) | is.na(Posprobtable$sp2_diet),
                                      NA_character_, paste(pmin(Posprobtable$sp1_diet, Posprobtable$sp2_diet),
                                                           pmax(Posprobtable$sp1_diet, Posprobtable$sp2_diet),
                                                           sep="-"))
Negprobtable$sp1_diet = as.character(Negprobtable$sp1_diet)
Negprobtable$sp2_diet = as.character(Negprobtable$sp2_diet)
Negprobtable$diet_pair = ifelse(is.na(Negprobtable$sp1_diet) | is.na(Negprobtable$sp2_diet),
                                NA_character_, paste(pmin(Negprobtable$sp1_diet, Negprobtable$sp2_diet),
                                                     pmax(Negprobtable$sp1_diet, Negprobtable$sp2_diet),
                                                     sep="-"))

Posprobtable$species_pair= paste(pmin(as.character(Posprobtable$sp1_name),as.character(Posprobtable$sp2_name)),
                                 pmax(as.character(Posprobtable$sp1_name),as.character(Posprobtable$sp2_name)),
                                 sep="-")
Negprobtable$species_pair= paste(pmin(as.character(Negprobtable$sp1_name),as.character(Negprobtable$sp2_name)),
                                 pmax(as.character(Negprobtable$sp1_name),as.character(Negprobtable$sp2_name)),
                                 sep="-")

CarnPosprobtable=Posprobtable%>% filter(diet_pair %in% c("Browser-Carnivore","Carnivore-Carnivore","Carnivore-Frugivore","Carnivore-Grazer","Carnivore-Insectivore","Carnivore-Omnivore"))
