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
#Making sure the categorical timebins are in order
TotalProportionsSig$X2=factor(TotalProportionsSig$X2, levels=c("30000","25000","20000","15000","10000","5000"))
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
#Making sure the categorical timebins are in order
TotalProportionsAgg$X2=factor(TotalProportionsAgg$X2, levels=c("30000","25000","20000","15000","10000","5000"))
#Plot for proportion of aggregated pairs. stat="identity for single values
#Used cowplot theme, added labels and title. Changed element text to make it look better, also changed axis range to illuminate trend
ggplot(TotalProportionsAgg, aes(x=X2,y=X1))+geom_bar(stat="identity", fill="tomato", width=0.5, color="black",linewidth=1)+
  theme_cowplot()+labs(x="Time",y="Proportion of Aggregated Pairs", title="Proportion of Aggregated Pairs Through Time in North America")+
  theme(plot.title=element_text(size=11.5), axis.title.y=element_text(size=10),
           axis.title.x=element_text(size=10),axis.text.y=element_text(size=10))+
  coord_cartesian(ylim=c(0.4,0.9))
#Now I am inserting the diet data
SpeciesDiets=read.csv("data/TraitdataTestVersion.csv")

#For each timebin I am using the Cooccur function of prob table to get data frames with the species pairs and whether the pair is significant 
probtable1=prob.table(Cooctotals[[1]])
#Adding the column with diet variables to the prob table from the diet data. Made two diet columns for both species in pairs so species names match.
probtable1=probtable1%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
#In the probtable, positive and negative pairs are determined by whether the p_gt(positive) or p_lt(negative) values are significant (below an alpha value of 0.05)
probtable1$Pairtype[probtable1$p_gt<=0.05] = "Positive"
probtable1$Pairtype[probtable1$p_lt<=0.05]="Negative"
probtable1$Pairtype[probtable1$p_gt>0.05 & probtable1$p_lt>0.05]="Random"
#Adding a column for time as I will be combining all the timebins later.
probtable1$time=30000

#Repeated these steps for the other 5 timebins. Couldn't figure out how to make this into a loop.
probtable2=prob.table(Cooctotals[[2]])
probtable2=probtable2%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable2$Pairtype[probtable2$p_gt<=0.05] = "Positive"
probtable2$Pairtype[probtable2$p_lt<=0.05]="Negative"
probtable2$Pairtype[probtable2$p_gt>0.05 & probtable2$p_lt>0.05]="Random"
probtable2$time=25000

probtable3=prob.table(Cooctotals[[3]])
probtable3=probtable3%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable3$Pairtype[probtable3$p_gt<=0.05] = "Positive"
probtable3$Pairtype[probtable3$p_lt<=0.05]="Negative"
probtable3$Pairtype[probtable3$p_gt>0.05 & probtable3$p_lt>0.05]="Random"
probtable3$time=20000

probtable4=prob.table(Cooctotals[[4]])
probtable4=probtable4%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable4$Pairtype[probtable4$p_gt<=0.05] = "Positive"
probtable4$Pairtype[probtable4$p_lt<=0.05]="Negative"
probtable4$Pairtype[probtable4$p_gt>0.05 & probtable4$p_lt>0.05]="Random"
probtable4$time =15000

probtable5=prob.table(Cooctotals[[5]])
probtable5=probtable5%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable5$Pairtype[probtable5$p_gt<=0.05] = "Positive"
probtable5$Pairtype[probtable5$p_lt<=0.05]="Negative"
probtable5$Pairtype[probtable5$p_gt>0.05 & probtable5$p_lt>0.05]="Random"
probtable5$time=10000

probtable6=prob.table(Cooctotals[[6]])
probtable6=probtable6%>%left_join(SpeciesDiets %>% select(sppname, Diet), by=join_by("sp1_name"=="sppname"))%>%rename(sp1_diet=Diet)%>%
  left_join(SpeciesDiets%>% select(sppname, Diet), by=join_by("sp2_name"=="sppname"))%>%rename(sp2_diet=Diet)
probtable6$Pairtype[probtable6$p_gt<=0.05] = "Positive"
probtable6$Pairtype[probtable6$p_lt<=0.05]="Negative"
probtable6$Pairtype[probtable6$p_gt>0.05 & probtable6$p_lt>0.05]="Random"
probtable6$time=5000

#Now that all the diet datas and pair types are assigned I can combine all the prob tables into a list.
Probtablecomb= list(probtable1,probtable2,probtable3,probtable4,probtable5,probtable6)
#First I want to get rid of all the pairs that are random as these are not part of my analysis. 
Probtablecomb= lapply(Probtablecomb,FUN=function(x)(x[x$Pairtype != "Random",]))
#Now that I have just positive and negative pairs, I can split the list into a nested list with separate dataframes for positive and negative pairs
Splitprobtable= lapply(Probtablecomb, FUN=function(x)(split(x, x$Pairtype)))
#Now I am pulling out the Positive and negative dataframes from the nested list and making them into a new individual dataframe with just positive pairs
Posprobtable= lapply(Splitprobtable, `[[`, "Positive")
Posprobtable=bind_rows(Posprobtable, .id="source")
#Same thing for negative pairs
Negprobtable= lapply(Splitprobtable, `[[`, "Negative")
Negprobtable=bind_rows(Negprobtable, .id="source")

#Next what I need to do is combine the diets of each species to a new column that records the diet pair, as in both diets combined
#First I turn the diets to characters from factors
Posprobtable$sp1_diet = as.character(Posprobtable$sp1_diet)
Posprobtable$sp2_diet = as.character(Posprobtable$sp2_diet)
#Then I combine them into a new column. The pmin and pmax organize them alphabetically so that Herbivore-Carnivore would be the same as Carnivore-Herbivore
#I used chatgpt to help here and it showed me how to use pmin and pmax and how to set up that code, I asked it how I could combine the data while preserving different orders being the same.
Posprobtable$diet_pair = ifelse(is.na(Posprobtable$sp1_diet) | is.na(Posprobtable$sp2_diet),
                                      NA_character_, paste(pmin(Posprobtable$sp1_diet, Posprobtable$sp2_diet),
                                                           pmax(Posprobtable$sp1_diet, Posprobtable$sp2_diet),
                                                           sep="-"))
#Then I do the same thing for the negative df.
Negprobtable$sp1_diet = as.character(Negprobtable$sp1_diet)
Negprobtable$sp2_diet = as.character(Negprobtable$sp2_diet)
Negprobtable$diet_pair = ifelse(is.na(Negprobtable$sp1_diet) | is.na(Negprobtable$sp2_diet),
                                NA_character_, paste(pmin(Negprobtable$sp1_diet, Negprobtable$sp2_diet),
                                                     pmax(Negprobtable$sp1_diet, Negprobtable$sp2_diet),
                                                     sep="-"))
#I also need to add species pairs as its own column with the names of each species, so I used pmin and pmax again to order them alphabetically
Posprobtable$species_pair= paste(pmin(as.character(Posprobtable$sp1_name),as.character(Posprobtable$sp2_name)),
                                 pmax(as.character(Posprobtable$sp1_name),as.character(Posprobtable$sp2_name)),
                                 sep="-")
#Repeated the same thing for negative pairs
Negprobtable$species_pair= paste(pmin(as.character(Negprobtable$sp1_name),as.character(Negprobtable$sp2_name)),
                                 pmax(as.character(Negprobtable$sp1_name),as.character(Negprobtable$sp2_name)),
                                 sep="-")
#Because the species names are big, I want to turn them into numbers so that it is not super cluttered on the graph
Posprobtable = Posprobtable%>% mutate(species_pairnum=as.numeric(factor(species_pair)))
Negprobtable=Negprobtable%>% mutate(species_pairnum=as.numeric(factor(species_pair)))


Comprobtable=bind_rows(Posprobtable, Negprobtable)
Comprobtable=Comprobtable%>% group_by(time, diet_pair, Pairtype)%>% summarise(pairnum=n(),.groups="drop")
Comprobtable=Comprobtable%>% complete(time, diet_pair, Pairtype, fill=list(pairnum=0))
Totalpairs=Comprobtable%>% group_by(time) %>% summarise(totalsig=sum(pairnum), .groups="drop")
Comprobtable=Comprobtable%>% left_join(Totalpairs, by="time")%>% mutate(proportion=pairnum/totalsig)
Comprobtable=Comprobtable%>% mutate(time=factor(time, levels=sort(unique(time),decreasing=T)))
Posprobtable=Comprobtable%>% filter(Pairtype=="Positive")
Negprobtable=Comprobtable%>% filter(Pairtype=="Negative")

#Because there are so many diet pairs and species pairs, I decided that I need to split the graphs into different subsections of diet pairs
#I made a carnivore prob table with the carnivore diet pairs and an herbivore prob table with the different herbivore pairs
CarnPosprobtable=Posprobtable%>% filter(diet_pair %in% c("Browser-Carnivore","Carnivore-Frugivore","Carnivore-Grazer"))
CarnNegprobtable=Negprobtable%>% filter(diet_pair%in% c("Browser-Carnivore","Carnivore-Frugivore","Carnivore-Grazer"))
HerbPosprobtable=Posprobtable%>% filter(diet_pair %in% c("Browser-Browser","Browser-Grazer","Grazer-Grazer", "Frugivore-Grazer","Browser-Frugivore","Frugivore-Frugivore"))
HerbNegprobtable=Negprobtable%>% filter(diet_pair %in% c("Browser-Browser","Browser-Grazer","Grazer-Grazer","Frugivore-Grazer","Browser-Frugivore","Frugivore-Frugivore"))
OmPosprobtable=Posprobtable%>% filter(diet_pair %in% c("Carnivore-Carnivore","Omnivore-Omnivore","Insectivore-Insectivore","Carnivore-Omnivore","Carnivore-Insectivore","Insectivore-Omnivore"))
OmNegprobtable=Negprobtable%>% filter(diet_pair %in% c("Carnivore-Carnivore","Omnivore-Omnivore","Insectivore-Insectivore","Carnivore-Omnivore","Carnivore-Insectivore","Insectivore-Omnivore"))

Hb1=ggplot(HerbPosprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("Proportion of Species Pairs")+theme(legend.position="none")+
  scale_x_discrete(drop=FALSE)
Hb2=ggplot(HerbNegprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("")+labs(fill="Trophic Pair")+scale_x_discrete(drop=FALSE)
hbc=plot_grid(Hb1,Hb2, ncol=2,align="hv",axis="tblr", labels=c("Aggregations","Segregations"), label_x=0.2, label_y=1.0)
ggsave("Herbplot.png", hbc, width=14,height=6,dpi=300)

Cb1=ggplot(CarnPosprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("Proportion of Species Pairs")+theme(legend.position="none")
Cb2=ggplot(CarnNegprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("")+labs(fill="Trophic Pair")
cbc=plot_grid(Cb1,Cb2, ncol=2,align="hv",axis="tblr", labels=c("Aggregations","Segregations"), label_x=0.2, label_y=1.0)
ggsave("Carnplot.png", cbc, width=14, height=6, dpi=300)

Ob1=ggplot(OmPosprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("Proportion of Species Pairs")+theme(legend.position="none")
Ob2=ggplot(OmNegprobtable, aes(x=time, y=proportion, fill=diet_pair))+geom_bar(stat="identity")+theme_cowplot()+
  xlab("Time")+ylab("")+labs(fill="Trophic Pair")
obc=plot_grid(Ob1,Ob2, ncol=2,align="hv",axis="tblr", labels=c("Aggregations","Segregations"), label_x=0.2, label_y=1.0)
ggsave("Omplot.png",obc,width=14,height=6,dpi=300)








#Obsolete code for this project, but may come in handy later
Hp1=ggplot(HerbPosprobtable, aes(x=time, y=species_pairnum, color=diet_pair, shape=diet_pair))+geom_point(size=2)+theme_cowplot()+
  xlab("Time")+ ylab("Species Pairs")+theme(legend.position="none")+
  scale_color_discrete(name="Diet Pair")+scale_shape_discrete(name="Diet Pair")+
  guides(color=guide_legend(), shape=guide_legend())+scale_x_continuous(trans="reverse")
Hp2=ggplot(HerbNegprobtable, aes(x=time, y=species_pairnum, color=diet_pair, shape=diet_pair))+geom_point(size=2)+theme_cowplot()+
  xlab("Time")+ylab("")+labs(color="Diet Pair")+ scale_color_discrete(name="Diet Pair")+scale_shape_discrete(name="Diet Pair")+
  guides(color=guide_legend(), shape=guide_legend())+scale_x_continuous(trans="reverse")
plot_grid(Hp1,Hp2,ncol=2,align="hv",axis="tblr", labels=c("Aggregations","Segregations"), label_x=0.15, label_y=1.0)

Cp1=ggplot(CarnPosprobtable, aes(x=time, y=species_pairnum, color=diet_pair, shape=diet_pair))+geom_point(size=2)+theme_cowplot()+
  xlab("Time")+ ylab("Species Pairs")+theme(legend.position="none")+
  scale_color_discrete(name="Diet Pair")+scale_shape_discrete(name="Diet Pair")+
  guides(color=guide_legend(), shape=guide_legend())+scale_x_continuous(trans="reverse")
Cp2=ggplot(CarnNegprobtable, aes(x=time, y=species_pairnum, color=diet_pair, shape=diet_pair))+geom_point(size=2)+theme_cowplot()+
  xlab("Time")+ylab("")+labs(color="Diet Pair")+ scale_color_discrete(name="Diet Pair")+scale_shape_discrete(name="Diet Pair")+
  guides(color=guide_legend(), shape=guide_legend())+scale_x_continuous(trans="reverse")
plot_grid(Cp1,Cp2,ncol=2,align="hv",axis="tblr", labels=c("Aggregations","Segregations"), label_x=0.3, label_y=1.0)
