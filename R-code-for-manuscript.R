library(ape)
library(dplyr)
library(expss)
library(ggplot2)
library(ggtree)
library(lme4)
library(lmerTest)
library(tidyr)
library(vegan)
library(rstudioapi)
library(lsr)
library(MASS)
library(phangorn)
library(phytools)


setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()

##### Input Data #####
InputOverlayMeasures <- read.table(file = 'Data/Final.Overlays.tsv', sep = '\t', header = TRUE) # use a quote character that I know isn't used to avoid some problems
InputAntibioticMeasures <- read.table(file = 'Data/Final.Antibiotics.tsv', sep = '\t', header = TRUE) # use a quote character that I know isn't used to avoid some problems
InputAlignment <- read.FASTA("Final.Seqs.fasta")
LigninUseInput <- read.table(file = 'Data/IsolatesTestedForLigninGrowth.tsv', sep = '\t', header = TRUE,quote = '$')

##### Universal Variables #####
FunctionMeanSE <- function(measurements){ m <- mean(measurements) ; se <- sd(measurements)/sqrt(length(measurements)) ; return(paste(m,"±",se))  }
MinimumInhibition <- 2
AllOverlays <- unique(InputOverlayMeasures$Overlay)
AllDots <- unique(InputOverlayMeasures$Dot)
AllCodes <- unique(InputOverlayMeasures$Code)
AllAntibiotics <- unique(InputAntibioticMeasures$Antibiotic)

##### Percent of Strep which use lignin as a sole carbon source #####
paste(sum(LigninUseInput$Lignin),' of ',nrow(LigninUseInput),' use lignin as a sole carbon source (',round(sum(LigninUseInput$Lignin)/nrow(LigninUseInput)*100,digits=1),'%)',sep='')
paste(sum(LigninUseInput$Lignin[LigninUseInput$NPK=='fert+']),' of ',nrow(LigninUseInput[LigninUseInput$NPK=='fert+',]),' from fert+ use lignin as a sole carbon source (',round(sum(LigninUseInput$Lignin[LigninUseInput$NPK=='fert+'])/nrow(LigninUseInput[LigninUseInput$NPK=='fert+',])*100,digits=1),'%)',sep='')
paste(sum(LigninUseInput$Lignin[LigninUseInput$NPK=='fert-']),' of ',nrow(LigninUseInput[LigninUseInput$NPK=='fert-',]),' from fert- use lignin as a sole carbon source (',round(sum(LigninUseInput$Lignin[LigninUseInput$NPK=='fert-'])/nrow(LigninUseInput[LigninUseInput$NPK=='fert-',])*100,digits=1),'%)',sep='')

# Table 1
glmer(Lignin~NPK*Res+(1|Block),data=LigninUseInput,family=binomial) %>% summary

##### Calculate means of all dot and overlay combinations ##### 
MeanByComboDF <- 
  aggregate(x=list(Inhib=InputOverlayMeasures$rmean),
            by=list(Overlay=InputOverlayMeasures$Overlay,Dot=InputOverlayMeasures$Dot,Code=InputOverlayMeasures$Code,
                    UsesLignin=InputOverlayMeasures$UsesLignin,plot=InputOverlayMeasures$plot,block=InputOverlayMeasures$block,
                    fertilizer=InputOverlayMeasures$fertilizer,residue=InputOverlayMeasures$residue),
            FUN=mean
  )
MeanByComboDF$Inhib <- ifelse(MeanByComboDF$Inhib<MinimumInhibition,0,MeanByComboDF$Inhib)


##### Select those which inhibit #####
Inhibitors.MeanByComboDF <- MeanByComboDF[MeanByComboDF$Inhib>0,]

##### Normalize the inhibition data #####
hist(Inhibitors.MeanByComboDF$Inhib) ; shapiro.test(Inhibitors.MeanByComboDF$Inhib)
shapiro.test(log(Inhibitors.MeanByComboDF$Inhib)) ; hist(log(Inhibitors.MeanByComboDF$Inhib))
shapiro.test(1/sqrt(Inhibitors.MeanByComboDF$Inhib)) ; hist(1/sqrt(Inhibitors.MeanByComboDF$Inhib))
# I also tested exp, squared, cubed, square root, and cube root.
Inhibitors.MeanByComboDF$Log.Inhib <- log(Inhibitors.MeanByComboDF$Inhib) # Not quite normal, but within the bounds of what a linear model can handle
Inhibitors.MeanByComboDF$InvSqrt.Inhib <- 1/sqrt(Inhibitors.MeanByComboDF$Inhib) # normal

##### Zone of Inhibition #####

# Table 2, part 2
lmer(InvSqrt.Inhib~UsesLignin*fertilizer*residue+(1|block),Inhibitors.MeanByComboDF) %>% anova

# Table S3
lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$fertilizer=='NPK+'& Inhibitors.MeanByComboDF$residue=='res+',]) %>% anova
lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$fertilizer=='NPK-'& Inhibitors.MeanByComboDF$residue=='res+',]) %>% anova
lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$fertilizer=='NPK+'& Inhibitors.MeanByComboDF$residue=='res-',]) %>% anova
lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$fertilizer=='NPK-'& Inhibitors.MeanByComboDF$residue=='res-',]) %>% anova

lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res+',]) %>% anova
lmer(InvSqrt.Inhib~UsesLignin+(1|block),Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res-',]) %>% anova



# mean and se for Table S3
FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$fertilizer=='NPK+'& Inhibitors.MeanByComboDF$residue=='res+'])
FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$fertilizer=='NPK-'& Inhibitors.MeanByComboDF$residue=='res+'])
FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$fertilizer=='NPK+'& Inhibitors.MeanByComboDF$residue=='res-'])
FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$fertilizer=='NPK-'& Inhibitors.MeanByComboDF$residue=='res-'])

FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$residue=='res+'])
FunctionMeanSE(Inhibitors.MeanByComboDF$Inhib[Inhibitors.MeanByComboDF$residue=='res-'])


# Calculate Cohen's D for the impact of lignin use on inhibition strength for res+/fert+ and res+/fert-
# res+/fert+
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res+' & Inhibitors.MeanByComboDF$fertilizer=='NPK+',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)

# res+/fert-
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res+' & Inhibitors.MeanByComboDF$fertilizer=='NPK-',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)

# res-/fert+
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res-' & Inhibitors.MeanByComboDF$fertilizer=='NPK+',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)

# res-/fert-
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res-' & Inhibitors.MeanByComboDF$fertilizer=='NPK-',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)

# res+
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res+',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)

# res-
TempDF <- Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$residue=='res-',];cohensD(x=TempDF$Log.Inhib[TempDF$UsesLignin==1],y=TempDF$Log.Inhib[TempDF$UsesLignin==0]);rm(TempDF)


# set up for Figure 2
FacetedLogInhibitSignificanceDF <- data.frame(FullResidueName=c('Residue retained','Residue retained'),FullFertilizerName=c('Fertilizer added','No fertilizer'),label=c('*','***'))

Inhibitors.MeanByComboDF$FullResidueName <- ifelse(Inhibitors.MeanByComboDF$residue=='res+','Residue retained','Residue removed')
Inhibitors.MeanByComboDF$FullFertilizerName <- ifelse(Inhibitors.MeanByComboDF$fertilizer=='NPK+','Fertilizer added','No fertilizer')

# Figure 2
Figure2.QuadDiagram <- 
ggplot(Inhibitors.MeanByComboDF,aes(y=log10(Inhib),x=as.factor(UsesLignin),fill=as.factor(UsesLignin)))+theme_bw()+
  geom_hline(yintercept = 0)+
  ylab('Log strength of inhibition (mm)')+xlab('')+
  coord_cartesian(ylim=c(0,1.25))+
  facet_grid(FullResidueName~FullFertilizerName)+
  scale_fill_manual(values=c('1'='grey','0'='white'),labels=c('1'='Yes','0'='No'),name='Utilizes lignin as\nsole carbon source')+
  stat_summary(geom='bar',fun=mean,color='black')+stat_summary(fun.data = mean_se, geom = "errorbar",width=0.25)+
  theme(panel.grid = element_blank())+
  geom_segment(data=FacetedLogInhibitSignificanceDF,inherit.aes = F,aes(x=1,xend=2,y=1.075,yend=1.075),size=0.8)+
  geom_text(data=FacetedLogInhibitSignificanceDF[1,],inherit.aes = F,aes(x=1.5,y=1.15,label='*'),size=8)+
  geom_text(data=FacetedLogInhibitSignificanceDF[2,],inherit.aes = F,aes(x=1.5,y=1.15,label='***'),size=8)+
  scale_x_discrete(labels=c('1'='Uses\nlignin','0'="Doesn't use\nlignin"))+
  theme(legend.position = 'none')+
  theme(axis.text = element_text(color='black'))+
  theme(text = element_text(face="bold", size=12), strip.text = element_text(face="bold", size=10))

Figure2.QuadDiagram
# ggsave(filename='fig2.jpeg',plot=Figure2.QuadDiagram,dpi=1000,width=5,height=5,units='in') # commented to keep from accidently overwriting

##### Intensity of inhibition x cultural management, facet by lignin use #####

head(Inhibitors.MeanByComboDF)

shapiro.test(Inhibitors.MeanByComboDF$InvSqrt.Inhib[Inhibitors.MeanByComboDF$UsesLignin==1]);hist(Inhibitors.MeanByComboDF$InvSqrt.Inhib[Inhibitors.MeanByComboDF$UsesLignin==1])
shapiro.test(Inhibitors.MeanByComboDF$InvSqrt.Inhib[Inhibitors.MeanByComboDF$UsesLignin==0]);hist(Inhibitors.MeanByComboDF$InvSqrt.Inhib[Inhibitors.MeanByComboDF$UsesLignin==0])
# Doesn't meet assumptions of lm, use a robust model, look at the histogram of UsesLignin==0

# Table S4
MASS::rlm(Inhib~fertilizer*residue,data=Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$UsesLignin==1,]) %>% summary
repmod::rob.pvals(MASS::rlm(Inhib~fertilizer*residue,data=Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$UsesLignin==1,]))

# Table S4
MASS::rlm(Inhib~fertilizer*residue,data=Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$UsesLignin==0,]) %>% summary
repmod::rob.pvals(MASS::rlm(Inhib~fertilizer*residue,data=Inhibitors.MeanByComboDF[Inhibitors.MeanByComboDF$UsesLignin==0,]))


##### Binomial model of the number of standards inhibited #####
MeanByComboDF$Binary.DoesInhibit <- ifelse(MeanByComboDF$Inhib<MinimumInhibition,0,1)

PropInhibitedDF <-
  aggregate(x=list(NumberInhibited=MeanByComboDF$Binary.DoesInhibit),
            by=list(Dot=MeanByComboDF$Dot,UsesLignin=MeanByComboDF$UsesLignin,plot=MeanByComboDF$plot,
                    block=MeanByComboDF$block,fertilizer=MeanByComboDF$fertilizer,residue=MeanByComboDF$residue),
            FUN=sum
  )
head(PropInhibitedDF)

PropInhibitedDF$NumberInhibited %>% min
PropInhibitedDF$NumberInhibited %>% max
sum(PropInhibitedDF$NumberInhibited>=7)
sum(PropInhibitedDF$NumberInhibited<=3)

# Table 2, part 1
glmer(Binary.DoesInhibit~UsesLignin*fertilizer*residue+(1|block),data=MeanByComboDF,family=binomial) %>% summary
  
##### Resistance data #####

head(InputAntibioticMeasures)
SimplifiedAntibioticsDF <- aggregate(x=list(rmean=InputAntibioticMeasures$rmean),
                                     by=list(Isolate=InputAntibioticMeasures$Isolate,
                                             Antibiotics=InputAntibioticMeasures$Antibiotic,
                                             UsesLignin=InputAntibioticMeasures$UsesLignin,
                                             plot=InputAntibioticMeasures$plot,
                                             block=InputAntibioticMeasures$block,
                                             fertilizer=InputAntibioticMeasures$fertilizer,
                                             residue=InputAntibioticMeasures$residue,
                                             FullName=InputAntibioticMeasures$FullName),
                                     FUN=mean,na.rm=T)
SimplifiedAntibioticsDF$rmean <- ifelse(SimplifiedAntibioticsDF$rmean<MinimumInhibition,0,SimplifiedAntibioticsDF$rmean)
SimplifiedAntibioticsDF$BinaryResists <- ifelse(SimplifiedAntibioticsDF$rmean>MinimumInhibition,0,1)
AntibioticResistanceSummary <- aggregate(x=list(NumberResisted=SimplifiedAntibioticsDF$BinaryResists),by=list(Isolate=SimplifiedAntibioticsDF$Isolate),FUN=sum)

head(SimplifiedAntibioticsDF)
head(AntibioticResistanceSummary)

AntibioticResistanceSummary$NumberResisted%>%max # max number of antibiotics completely resisted
AntibioticResistanceSummary$NumberResisted%>%min # min number of antibiotics completely resisted
sum(AntibioticResistanceSummary$NumberResisted<=1) # number of isolates resisting 1 or fewer antibiotics

# Table 4, part 1
glmer(BinaryResists ~UsesLignin*fertilizer*residue+(1|block),data=SimplifiedAntibioticsDF,family=binomial) %>% summary

SimplifiedAntibioticsDF.SensetiveInteractions <- SimplifiedAntibioticsDF[SimplifiedAntibioticsDF$BinaryResists==0,]
shapiro.test(SimplifiedAntibioticsDF.SensetiveInteractions$rmean^(1/2)) ; hist(SimplifiedAntibioticsDF.SensetiveInteractions$rmean^1/2)

# Table 4, part 2
lmer(sqrt(rmean)~as.factor(UsesLignin)*fertilizer*residue+(1|Antibiotics),SimplifiedAntibioticsDF.SensetiveInteractions) %>% anova

# presented in text and informs Figure S1
for (x in unique(SimplifiedAntibioticsDF.SensetiveInteractions$Antibiotics)){
  if(x==first(unique(SimplifiedAntibioticsDF.SensetiveInteractions$Antibiotics))){ SigAntiList <- vector()  }
  TempDF <- SimplifiedAntibioticsDF.SensetiveInteractions[SimplifiedAntibioticsDF.SensetiveInteractions$Antibiotic==x,]
  Bigger <- ifelse(mean(TempDF$rmean[TempDF$UsesLignin==1])>mean(TempDF$rmean[TempDF$UsesLignin==0]),'Lignin users more sens','Lignin users less sens')
  Res <- lm(rmean~as.factor(UsesLignin),data=TempDF) %>% anova
  print(paste(x,round(Res$`Pr(>F)`[1],3),Bigger))
  if(  Res$`Pr(>F)`[1]<0.05){ SigAntiList <- c(SigAntiList,x)  }
  if(x==last(unique(SimplifiedAntibioticsDF.SensetiveInteractions$Antibiotics))){
    print("");print("");print("Significant antibiotics by lignin use:");print('Codes:');print(SigAntiList)
    print('');SimplifiedAntibioticsDF.SensetiveInteractions[SimplifiedAntibioticsDF.SensetiveInteractions$Antibiotics%in%SigAntiList,]%>%dplyr::select(FullName)%>%unique%>%print
    rm(SigAntiList)
  }
  rm(x,TempDF,Bigger,Res)
}

# Figure S1
AntiSigDF <- data.frame(FullName=c('Chloramphenicol','Novobiocin','Rifampin','Tetracycline') ,
                        vert=c(15.5,15.5,11,11),vert2=c(15.5,15.5,11,11))

# 725 x 400 when copied
FigureS1.Anti <-
ggplot(SimplifiedAntibioticsDF.SensetiveInteractions,aes(y=rmean,x=as.factor(UsesLignin),fill=as.factor(UsesLignin)))+theme_bw()+
  geom_hline(yintercept = 0)+facet_wrap(~FullName)+
  scale_fill_manual(values=c('1'='grey','0'='white'),labels=c('1'='Yes','0'='No'),name='Utilizes lignin as\nsole carbon source')+
  scale_x_discrete(labels=c('1'='Uses\nlignin','0'="Doesn't use\nlignin"))+
  theme(legend.position = 'none')+coord_cartesian(ylim=c(0,19))+
  stat_summary(geom='bar',color='black',fun=mean)+
  stat_summary(geom='errorbar',width=0.15,fun.data=mean_se)+
  geom_segment(inherit.aes = F,data=AntiSigDF,aes(x=1,xend=2,y=vert,yend=vert),linewidth=0.75)+
  geom_text(inherit.aes=F,data=AntiSigDF,aes(x=1.5,y=vert+1,label='*'),size=7)+
  theme(axis.text = element_text(color='black',face='bold'))+
  theme(strip.text = element_text(face="bold"))+
  theme(panel.grid = element_blank())+xlab('')+ylab('Average zone of inhibition (mm)')

FigureS1.Anti
# ggsave(filename='figS1.jpeg',plot=FigureS1.Anti,dpi=1000,width=7.3,height=4,units='in') # commented to keep from accidently overwriting


##### Biolog data #####
InputBiolog <- read.table(file = 'Data/Final.Biolog.tsv', sep = '\t', header = TRUE)
WideBiolog <- InputBiolog %>% dplyr::select(Isolate,Well,FinalOD) %>% pivot_wider(names_from=Well,values_from=FinalOD) %>% as.data.frame ; rownames(WideBiolog) <- WideBiolog$Isolate ; WideBiolog$Isolate <- NULL
BiologEucDist <- vegdist(WideBiolog,method='euclidean') %>% as.matrix()
BiologMeta <- InputBiolog %>% dplyr::select(Isolate,plot,block,residue,fertilizer,UsesLignin) %>% unique
if(!all(rownames(BiologEucDist)==BiologMeta$Isolate)){stop('Biolog row names dont match biolog meta names')}
for (x in 1:nrow(BiologMeta)){
  TempDF <- InputBiolog[InputBiolog$Isolate==BiologMeta$Isolate[x],]
  BiologMeta$NicheWidth[x] <- sum(TempDF$FinalOD>0)
  BiologMeta$AvgGrowth[x] <- (sum(TempDF$FinalOD)/sum(TempDF$FinalOD>0))
  rm(TempDF,x)
}
head(BiologMeta)

# table 3, part 1
set.seed(0507251434);adonis2(BiologEucDist~UsesLignin*fertilizer*residue,strata=BiologMeta$block,by='terms',na.action=na.omit,permutations=10000,data=BiologMeta)

# table 3, part 2
lmer(NicheWidth~UsesLignin*fertilizer*residue+(1|block),data=BiologMeta) %>% anova

# table 3, part 3
hist(BiologMeta$AvgGrowth)
lmer(AvgGrowth~UsesLignin*fertilizer*residue+(1|block),data=BiologMeta) %>% anova

# posthoc in text, tests change in OD of each well by lignin use, applies a bonferroni correction n=95
for (x in unique(InputBiolog$Well)){
  if(x==first(unique(InputBiolog$Well))){ BiologPosthocSig <- data.frame(Well=vector(),pval=vector()) }
  if(sum(InputBiolog$FinalOD[InputBiolog$Well==x]>0)<5){next} # minimum number that use to test
  Res <- lm(FinalOD~UsesLignin,InputBiolog[InputBiolog$Well==x,]) %>% anova
  # Res <- lmer(FinalOD~UsesLignin+(1|block),InputBiolog[InputBiolog$Well==x,]) %>% anova
  BiologPosthocSig <- rbind(BiologPosthocSig,data.frame(Well=x,pval=Res$`Pr(>F)`[1]))
  if(x==last(unique(InputBiolog$Well))){ BiologPosthocSig$qval <- p.adjust(BiologPosthocSig$pval,n=95,method='bonferroni');print(BiologPosthocSig[which(BiologPosthocSig$qval<0.05),])}
  rm(x,Res)
}

# Binomial model to determine if any nutrients drive niche width by lignin use
for (x in unique(InputBiolog$Well)){
  if(x==first(unique(InputBiolog$Well))){ BiologPosthocSigBinomial <- data.frame(Well=vector(),pval=vector()) }
  if(sum(InputBiolog$FinalOD[InputBiolog$Well==x]>0)<5){next} # minimum number that use to test
  TempDF <- InputBiolog[InputBiolog$Well==x,] ; TempDF$Binary <- ifelse(TempDF$FinalOD>0,1,0)
  Res <- glm(Binary~UsesLignin,TempDF,family=binomial) %>% summary
  BiologPosthocSigBinomial <- rbind(BiologPosthocSigBinomial,data.frame(Well=x,pval=Res$coefficients[8]))
  if(x==last(unique(InputBiolog$Well))){ BiologPosthocSigBinomial$qval <- p.adjust(BiologPosthocSigBinomial$pval,n=95,method='bonferroni');print(BiologPosthocSigBinomial[which(BiologPosthocSigBinomial$qval<0.05),])}
  rm(x,Res,TempDF)
  # Some of the glm did not converge, from my inspection we can interpert these as not significant
}

# For how many carbon sources is the OD higher, even if not significant
for (x in unique(InputBiolog$Well)){
  if(x==first(unique(InputBiolog$Well))){ CompareUseDF <- data.frame(Well=vector(),Bigger=vector())  }
  TempDF <- InputBiolog[InputBiolog$Well==x,]
  if(sum(TempDF$FinalOD>0)<(nrow(TempDF)*.1)){temp<-'Not testable'}else{
    yes <- TempDF[TempDF$UsesLignin==1,]; no <- TempDF[TempDF$UsesLignin==0,]
    temp <- ifelse(mean(yes$FinalOD)>mean(no$FinalOD),"Users more efficient",'Non-users more efficient')
    rm(yes,no)
  }
  CompareUseDF <- rbind(CompareUseDF, data.frame(Well=x,Bigger=temp) )
  rm(x,TempDF,temp)
}
table(CompareUseDF$Bigger)

# Prep for Figure 4
BiologCalculatePCOA <- pcoa(BiologEucDist)
pcoaDF <- BiologCalculatePCOA$vectors %>% as.data.frame()
RelativeEig <- round(BiologCalculatePCOA$values$Relative_eig[1:2]*100,2)
head(pcoaDF)
if(!all(rownames(pcoaDF)==BiologMeta$Isolate)){stop("Row names do not match while making pcoa")}else{ pcoaDF2 <- cbind(pcoaDF,BiologMeta) }

# Figure 3
Figure3.PCOA <- 
ggplot(pcoaDF2,aes(x=Axis.1,y=Axis.2,fill=as.factor(UsesLignin)))+theme_bw()+
  coord_fixed()+theme(legend.position = 'none')+
  geom_point(pch=21)+stat_ellipse(aes(linetype=as.factor(UsesLignin)))+
  scale_fill_manual(values=c('1'='grey50','0'='white'))+
  scale_linetype_manual(values=c('1'='solid','0'='dashed'))+
  xlab(paste("Axis 1 (",RelativeEig[1],"%)",sep=''))+
  ylab(paste("Axis 2 (",RelativeEig[2],"%)",sep=''))

Figure3.PCOA
# ggsave(filename='fig3.jpeg',plot=Figure3.PCOA,dpi=1000,width=4,height=5,units='in') # commented to keep from accidently overwriting

# Figure S2
head(InputBiolog)

FigureS2.Carbon <-
ggplot(InputBiolog[InputBiolog$Well=='B10',],aes(y=FinalOD,x=as.factor(UsesLignin),fill=as.factor(UsesLignin)))+theme_bw()+
  theme(legend.position = 'none')+ geom_violin()+
  scale_x_discrete(labels=c('1'='Uses\nlignin','0'="Doesn't use\nlignin"))+
  scale_fill_manual(values=c('1'='grey','0'='white'),labels=c('1'='Yes','0'='No'),name='Utilizes lignin as\nsole carbon source')+
  ylab('OD (590 nm)')+xlab('')+ggtitle('D-Gluconic Acid')+
  theme(axis.text.y = element_text(color='black',face='bold'))+
  theme(axis.text.x = element_text(color='black',face='bold',size=10))+
  theme(panel.grid = element_blank())

FigureS2.Carbon
ggsave(filename='figS2.jpeg',plot=FigureS2.Carbon,dpi=1000,width=2.5,height=4,units='in') # commented to keep from accidently overwriting

##### MUSCLE Alignment #####

if (1==2){ # run manually only
  if(!exists('mt')){ set.seed(06051019);mt <- modelTest(InputAlignment) } # Coded to not run again if it was already run, it takes a while to run; TrN is the best
  if(!exists('fit_bb')){ set.seed(06031506); fit_bb <- pml_bb(mt) } # Coded to not run again if it was already run, it takes a long time
  # write.pml(x=fit_bb,file='FinalTree') # commented out so the tree has to be saved by a manual choice, I won't accidently overwrite a tree
}
# Even with seeds set, I am getting a lot of variation run-to-run, especially when moving between computers
# To address this, I have saved the tree I have selected and I am loading the previously formed tree each time
InputTreePML <- read.tree(file='FinalTree_tree.nwk')

tree_rooted <- root(InputTreePML, outgroup = "MG1", resolve.root = TRUE)
TempDF <- data.frame(Isolate=as.vector(tree_rooted$tip.label))
temp <- LigninUseInput %>% dplyr::select(Isolate=Isolate,Res,NPK,Crop,Lignin,Block)
temp2 <- rbind(temp,  data.frame(Isolate='MG1',Res=NA,NPK=NA,Crop=NA,Lignin=NA,Block=NA) )
TreeMeta <- left_join(x=TempDF,y=temp2,by='Isolate') ; rm(TempDF,temp)
if(!all(TreeMeta$MG.Name==tree_rooted$tip.label) ){stop("tip names don't match meta names")}

# If an edge[,2] contains 1:length(tree_rooted$tip.label) then it is terminal
NoZeroTerminalsTree <- tree_rooted
# find the smallest non-zero terminal edge
for (x in 1:nrow(tree_rooted$edge)){
  if(x==1){TerminalLengths <- vector() }
  if(tree_rooted$edge[x,2]>length(tree_rooted$tip.label)){rm(x);next}
  if(tree_rooted$edge.length[x]==0){rm(x);next}
  TerminalLengths <- c(TerminalLengths,tree_rooted$edge.length[x] )
  rm(x)
}

AddBranchLength <- min(TerminalLengths)/(10000)

for (x in 1:nrow(tree_rooted$edge) ){
  if(x==1){ NoZeroTerminalsTree <- tree_rooted }
  if(NoZeroTerminalsTree$edge[,2][x]<=length(NoZeroTerminalsTree$tip.label)){ if(NoZeroTerminalsTree$edge.length[x]==0){NoZeroTerminalsTree$edge.length[x] <- AddBranchLength} }
  rm(x)
}

# Next drop the root because it is not a Streptomyces, it is only here to help the tree come together properly
DropRootTree <- drop.tip(NoZeroTerminalsTree,tip='MG1')
TreeMeta.DropRoot <- TreeMeta[TreeMeta$Isolate!='MG1',]
if(!all(DropRootTree$tip.label == TreeMeta.DropRoot$Isolate)){stop("Tip names do not match")}

Binary.Lignin <- TreeMeta.DropRoot$Lignin ; names(Binary.Lignin) <- TreeMeta.DropRoot$Isolate
Binary.Res <-   ifelse(TreeMeta.DropRoot$Res   == 'res+'  ,1,0); names(Binary.Res)<-TreeMeta.DropRoot$Isolate
Binary.NPK <-   ifelse(TreeMeta.DropRoot$NPK   == 'fert+' ,1,0); names(Binary.NPK)<-TreeMeta.DropRoot$Isolate
Binary.Block <- ifelse(TreeMeta.DropRoot$Block == 'left'  ,1,0); names(Binary.Block)<-TreeMeta.DropRoot$Isolate

# Table S2
phylosig(x=Binary.Res,tree=DropRootTree,nsim = 10000,test=T,method='lambda')
phylosig(x=Binary.NPK,tree=DropRootTree,nsim = 10000,test=T,method='lambda')
phylosig(x=Binary.Lignin,tree=DropRootTree,nsim = 10000,test=T,method='lambda')
phylosig(x=Binary.Block,tree=DropRootTree,nsim = 10000,test=T,method='lambda')

##### Make phylogenetic tree #####

TreeMeta2 <- TreeMeta
TreeMeta2$Lignin <- ifelse(TreeMeta2$Isolate=='MG1',2,TreeMeta2$Lignin)

# Figure 1
Figure1.Tree <-
ggtree(tree_rooted) %<+% TreeMeta2 +  # %<+% joins data to tree
  geom_tippoint(aes(fill = as.factor(Lignin)), size = 2,pch=21) +
  scale_fill_manual(values = c("1" = "green3", "0" = "red2",'2'='black')) +
  theme_tree2()+theme(legend.position = 'none')+ylim((-1),136)+
  theme(axis.line.x=element_blank(),axis.ticks=element_blank(),axis.text=element_blank())

Figure1.Tree
# ggsave(filename='fig1.jpeg',plot=Figure1.Tree,dpi=1000,width=6,height=12,units='in') # commented to keep from accidently overwriting
