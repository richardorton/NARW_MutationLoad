



##################################################
#                                                #  
#                      NARW                      #
#                                                #
##################################################


library(devtools)
library(ggplot2)
library(detectRUNS)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(pegas)
library(vcfR)
library(rCNV)
library(corrplot)
library(wesanderson)
library(ggforce)
library(ggdist)
library(gghalves)
library(qqman)
library(reshape2)
library(plyr)
library("adegenet")
library(ggridges)
library(gghighlight)
library(WebGestaltR)
library(patchwork)


###################################################
####  nucleotide diversity window-based      ######
###################################################

### read pi data for X ###

X_pi_NARW <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/pixy_output/X/pixy_pi.txt", header = T)
head(X_pi_NARW)

X_pi_NARW$CHROM <- 'X chromosome'



### read pi data for autosomes ###

autosomes_pi_NARW <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/pixy_output/autosomes/pixy_pi.txt", header = T)
head(autosomes_pi_NARW)

autosomes_pi_NARW$CHROM <- 'Autosomes'

### box plots ###


pi_NARW <- rbind(autosomes_pi_NARW, X_pi_NARW)
head(pi_NARW)

foura <- ggplot(pi_NARW, aes(x=CHROM, y=avg_pi, fill = pop)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 26, color = "black")) +
  ylab(expression("Nucleotide diveristy"~(pi))) +
  xlab("Genomic region") +
  ylim(0,.001)


NARW_pi_B <- na.omit(pi_NARW)
agg_NARW_pi <- aggregate(NARW_pi_B$avg_pi,  
                                      by=list(NARW_pi_B$pop,
                                              NARW_pi_B$CHROM), FUN=mean)


agg_NARW_pi

#################################
#Group.1   Group.2            x
#1    High Autosomes 0.0003998813
#2     Low Autosomes 0.0004015667
#3    High         X 0.0002443749
#4     Low         X 0.0002461747
#################################

NARW_pi.aov2 <- aov(avg_pi  ~ pop + CHROM, data = pi_NARW)
summary(NARW_pi.aov2)

#Df   Sum Sq   Mean Sq F value Pr(>F)    
#pop             1 0.000000 3.000e-08   0.098  0.754    
#CHROM           1 0.000056 5.643e-05 163.453 <2e-16 ***
#  Residuals   47207 0.016299 3.500e-07     


########################################################################################
#                              F - genome-wide                                         #
########################################################################################

############### read in data #######
## autosomes
NARW_F_autosomes_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.ALL.het", header = T)
NARW_F_autosomes_data$CHR <- 'Autosomes'

## X 
NARW_F_X_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/X.ALL.het", header = T)
NARW_F_X_data$CHR <- 'X chromosome'

NARW_F_data <- rbind(NARW_F_autosomes_data, NARW_F_X_data )
NARW_F_data
#### Fecundity ########
#EGL00252-1	Low
#EGL013-3qa	Low
#EGL183-1 	High
#EGL254-1	  Low
#EGL276-1	  Low
#EGL308-1a	Low
#EGL312-1a	High
#EGL336_1b	High
#SID179132	High
#SID181803	High
######################

NARW_F_data$fecundity_group <- c('Low','Low','High','Low','Low','Low',
                                 'High','High', 'High','High')

head(NARW_F_data)

agg_NARW_F_data <- aggregate(NARW_F_data$F,  
                         by=list(NARW_F_data$fecundity_group,
                                 NARW_F_data$CHR), FUN=mean)

agg_NARW_F_data

###############################
#Group.1   Group.2         x
#1    High Autosomes -0.069198
#2     Low Autosomes -0.074278
#3    High         X -0.314216
#4     Low         X -0.272808
###############################



NARW_F_data$Het <- (NARW_F_data$N_SITES - NARW_F_data$O.HOM.)/NARW_F_data$N_SITES
head(NARW_F_data)


NARW_het.aov2 <- aov(Het  ~ fecundity_group + CHR, data = NARW_F_data)
summary(NARW_het.aov2)

#Df  Sum Sq Mean Sq F value   Pr(>F)    
#fecundity_group  1 0.00018 0.00018   0.917    0.352    
#CHR              1 0.03838 0.03838 196.107 9.17e-11 ***
#  Residuals       17 0.00333 0.00020  


####### plot heterozygosity #########################################################

fourb <- ggplot(NARW_F_data, aes(x=CHR, y=Het, fill = fecundity_group)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Heterozygosity (genome-wide)") +
  xlab("Genomic region")




####### plot Wright's F #########################################################
ggplot(NARW_F_data, aes(x=CHR, y=F, fill = fecundity_group)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Inbreeding (Wright's F)") +
  xlab("Genomic region")
  theme(legend.position = "none") 




####### plot ROH #########################################################
autosomes_roh <- read.table("/Volumes/cetacea/Genetic_load/_roh/autosomes.population.roh", header = T)
NARW_autosomes_roh <- autosomes_roh[which(autosomes_roh$Population == "NARW"),] 
head(NARW_autosomes_roh)
unique(NARW_autosomes_roh$V2)

NARW_autosomes_roh$Fecundity_group <- case_when(
  NARW_autosomes_roh$V2 == "EGL00252-1" ~ "Low",
  NARW_autosomes_roh$V2 == "EGL013-3qa" ~ "Low",
  NARW_autosomes_roh$V2 == "EGL183-1" ~ "High",
  NARW_autosomes_roh$V2 == "EGL254-1" ~ "Low",
  NARW_autosomes_roh$V2 == "EGL276-1" ~ "Low",
  NARW_autosomes_roh$V2 == "EGL308-1a" ~ "Low",
  NARW_autosomes_roh$V2 == "EGL312-1a" ~ "High",
  NARW_autosomes_roh$V2 == "EGL336_1b" ~ "High",
  NARW_autosomes_roh$V2 == "SID179132" ~ "High",
  NARW_autosomes_roh$V2 == "SID181803" ~ "High",
  TRUE ~ "unknown"
)


unique(NARW_autosomes_roh$Fecundity_group)
NARW_autosomes_roh <- NARW_autosomes_roh[which(NARW_autosomes_roh$Fecundity_group == "Low" | NARW_autosomes_roh$Fecundity_group == "High" ),] 
head(NARW_autosomes_roh)

ggplot(NARW_autosomes_roh, aes(x = V6, y = Fecundity_group, fill = Fecundity_group)) + 
  geom_density_ridges() + 
  stat_density_ridges(quantile_lines = TRUE) +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45)) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Population") +
  scale_x_continuous(breaks=seq(0, 2000000, 500000),
                     labels=c("0", "0.5", "1.0", "1.5", "2.0")) +
  #scale_x_continuous(labels=scientific) + 
  #xlim("0", "2000000") +
  xlab("Autosomes - Length of ROH (Mb)") +
  theme(legend.position = "none")  


aggregate(NARW_autosomes_roh$V6,  by=list(NARW_autosomes_roh$Fecundity_group), FUN=mean)
table(NARW_autosomes_roh$V2)



X_roh <- read.table("/Volumes/cetacea/Genetic_load/_roh/X.population.roh", header = T)
X_roh_NARW <- X_roh[which(X_roh$Population == "NARW"),] 

X_roh_NARW$Fecundity_group <- case_when(
  X_roh_NARW$V2 == "EGL00252-1" ~ "Low",
  X_roh_NARW$V2 == "EGL013-3qa" ~ "Low",
  X_roh_NARW$V2 == "EGL183-1" ~ "High",
  X_roh_NARW$V2 == "EGL254-1" ~ "Low",
  X_roh_NARW$V2 == "EGL276-1" ~ "Low",
  X_roh_NARW$V2 == "EGL308-1a" ~ "Low",
  X_roh_NARW$V2 == "EGL312-1a" ~ "High",
  X_roh_NARW$V2 == "EGL336_1b" ~ "High",
  X_roh_NARW$V2 == "SID179132" ~ "High",
  X_roh_NARW$V2 == "SID181803" ~ "High",
)

unique(X_roh_NARW$Fecundity_group)
#X_roh_NARW <- X_roh[which(X_roh$Fecundity_group == "Low" | X_roh$Fecundity_group == "High" ),] 
head(X_roh_NARW)

ggplot(X_roh_NARW, aes(x = V6, y = Fecundity_group, fill = Fecundity_group)) + 
  geom_density_ridges() + 
  stat_density_ridges(quantile_lines = TRUE) +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45)) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Population") +
  scale_x_continuous(breaks=seq(0, 2000000, 500000),
                     labels=c("0", "0.5", "1.0", "1.5", "2.0")) +
  #scale_x_continuous(labels=scientific) + 
  #xlim("0", "2000000") +
  xlab("Autosomes - Length of ROH (Mb)") +
  theme(legend.position = "none")  









##############################################################
#########   Genetic load per site per impact     #############
##############################################################


##### read in X data #####
delta_load_X <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/X.NARW.delta_load", sep = '\t', header = T)
delta_load_X$CHR <- 'X chromosome'
head(delta_load_X)


##### read in autosomes data #####
delta_load_autosomes <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.NARW.delta_load", sep = '\t', header = T)
delta_load_autosomes$CHR <- 'Autosomes'
head(delta_load_autosomes)


delta_load <- rbind(delta_load_autosomes, delta_load_X)
unique(delta_load$estimate)
head(delta_load)



3### Total additive genetic load ###
deltaGenetic <- delta_load[which(delta_load$estimate == "deltaGenetic"),]
deltaGenetic <- deltaGenetic[which(deltaGenetic$value != 0),]
unique(deltaGenetic$Impact)
deltaGenetic$Impact <- factor(deltaGenetic$Impact, levels=c("Low", "Moderate", "High"))
head(deltaGenetic)

#######################################################################

FourC <- ggplot(deltaGenetic, aes(x=CHR, y=value, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "median",geom = "point", color = "black", position = position_dodge(.9)) +
  scale_fill_manual(values=c("gold", "orange3", "darkred")) +
  geom_hline(yintercept=0.06, linetype="dashed", color = "black", size=.25) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  ylab(expression(Delta~ "Total mutation load")) +
  xlab("Genomic region") +
  theme(legend.position = "none")    



### Realized genetic load ###
deltaRealized <- delta_load[which(delta_load$estimate == "deltaRealized"),]
deltaRealized <- deltaRealized[which(deltaRealized$value != 0),]
unique(deltaRealized$Impact)
deltaRealized$Impact <- factor(deltaRealized$Impact, levels=c("Low", "Moderate", "High"))
head(deltaRealized)



#foura <- 
ggplot(deltaRealized, aes(x=CHR, y=value, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=.25) +
  stat_summary(fun = "median",geom = "point", color = "black", position = position_dodge(.9)) +
  scale_fill_manual(values=c("gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  ylab(expression(Delta~ "Realized mutation load")) +
  xlab("Genomic region") +
  theme(legend.position = "none")  

deltaRealized_test <- aov(value ~ Impact + CHR, data = deltaRealized)
summary(deltaRealized_test)


#Df Sum Sq  Mean Sq F value Pr(>F)
#Impact         2  0.005 0.002341   0.514  0.598
#CHR            1  0.001 0.000711   0.156  0.693
#Residuals   6820 31.057 0.004554   



### Delta Masked load ###
deltaMasked <- delta_load[which(delta_load$estimate == "deltaMasked"),]
deltaMasked <- deltaMasked[which(deltaMasked$value != 0),]
unique(deltaMasked$Impact)
deltaMasked$Impact <- factor(deltaMasked$Impact, levels=c("Low", "Moderate", "High"))
head(deltaMasked)

ggplot(deltaMasked, aes(x=CHR, y=value, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "median",geom = "point", color = "black", position = position_dodge(.9)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=.25) +
  scale_fill_manual(values=c("gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  ylab("Delta Masked genetic load") +
  xlab("Fecundity")


### Delta Hom ###
deltaHom <- delta_load[which(delta_load$estimate == "deltaHom"),]
#deltaHom <- deltaHom[which(deltaHom$value != 0),]
deltaHom$Impact <- factor(deltaHom$Impact, levels=c("Modifier", "Low", "Moderate", "High"))
unique(deltaHom$Impact)
head(deltaHom)

ggplot(deltaHom, aes(x=CHR, y=value, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "mean",geom = "point", color = "black", position = position_dodge(.9)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=.25) +
  scale_fill_manual(values=c("gray50", "#E5C23A", "#D9892D", "#DE4530")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  ylab("Delta Derived homozygosity") +
  xlab("Fecundity")




### Delta allele frequency ###
deltaHom <- delta_load[which(delta_load$estimate == "deltaHom"),]
#deltaHom <- deltaHom[which(deltaHom$value != 0),]
deltaHom$Impact <- factor(deltaHom$Impact, levels=c("Modifier", "Low", "Moderate", "High"))
unique(deltaHom$Impact)
head(deltaHom)

ggplot(deltaHom, aes(x=CHR, y=value, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "mean",geom = "point", color = "black", position = position_dodge(.9)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=.25) +
  scale_fill_manual(values=c("gray50", "#E5C23A", "#D9892D", "#DE4530")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  ylab("Delta Derived homozygosity") +
  xlab("Fecundity")











###################################################
#########   Differentiation per site     ##########
###################################################

## autosomes ##
autosomes.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/autosomes.fst", sep = '\t', header = T)
head(autosomes.fst)

autosomes.fst <- autosomes.fst[which(autosomes.fst$Impact != "Impact"),]
autosomes.fst <- autosomes.fst[which(autosomes.fst$WEIR_AND_COCKERHAM_FST != "-nan"),]
autosomes.fst$WEIR_AND_COCKERHAM_FST <- as.numeric(as.character(autosomes.fst$WEIR_AND_COCKERHAM_FST))
autosomes.fst <- autosomes.fst[which(autosomes.fst$WEIR_AND_COCKERHAM_FST >= 0),]
autosomes.fst$CHR <- 'Autosomes'


## X ##
X.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/X.fst", sep = '\t', header = T)
head(autosomes.fst)

X.fst <- X.fst[which(X.fst$Impact != "Impact"),]
X.fst <- X.fst[which(X.fst$WEIR_AND_COCKERHAM_FST != "-nan"),]
X.fst$WEIR_AND_COCKERHAM_FST <- as.numeric(as.character(X.fst$WEIR_AND_COCKERHAM_FST))
X.fst <- X.fst[which(X.fst$WEIR_AND_COCKERHAM_FST >= 0),]
X.fst$CHR <- 'X chromosome'



NARW.fst <- rbind(autosomes.fst, X.fst)
head(NARW.fst)

agg_NARW_Fst_data <- aggregate(NARW.fst$WEIR_AND_COCKERHAM_FST,  
                                    by=list(NARW.fst$CHR, 
                                            NARW.fst$Impact), FUN=mean)

agg_NARW_Fst_data

Fst_test <- aov(WEIR_AND_COCKERHAM_FST ~ Impact + CHR, data = NARW.fst)
summary(Fst_test)

#Df Sum Sq Mean Sq F value Pr(>F)    
#Impact           3     10   3.179   239.9 <2e-16 ***
#  CHR              1     20  19.518  1472.8 <2e-16 ***
#  Residuals   254694   3375   0.013                   
# Signif. codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1


NARW.fst_NS <- NARW.fst[which(NARW.fst$Impact != "Modifier"),]
NARW.fst_NS$Impact <- factor(NARW.fst_NS$Impact , levels=c("Low", "Moderate", "High"))
head(NARW.fst_NS)
tail(NARW.fst_NS)
unique(NARW.fst_NS$Impact)


test <- NARW.fst_NS[which(NARW.fst_NS$CHR == "Autosomes"),]
quantile(test$WEIR_AND_COCKERHAM_FST, probs = .9)

lines <- data.frame(value=c(0.07,0.001),boxplot.nr=c("Autosomes", "X chromosome"))
Fst_threshold <- data.frame(value=0.25)

### violin or density plot ###
FourB <- ggplot(NARW.fst_NS, aes(x=CHR, y=WEIR_AND_COCKERHAM_FST, fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "median",geom = "point", color = "black", position = position_dodge(.9)) +
  scale_fill_manual(values=c("gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  #geom_boxplot(data=lines,aes(factor(boxplot.nr),value),
               #inherit.aes=FALSE,color="black",size=.25, linetype = "dashed") +
  geom_hline(data=Fst_threshold, linetype="dashed", size=.25, aes(yintercept=value)) +
  ylab(expression(paste(italic("Fst")))) +
  xlab("Genomic region") +
  theme(legend.position = "none") +
  ylim(0,1)





###################################################
#########        PCAdapt per site        ##########
###################################################


### X chromosome ###
Xpca <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/nps_pvalues_no_na_X", header = T)
Xpca$CHR <- 'X chromosome'
head(Xpca)

median(Xpca$data_pcadapt.pvalues)
#-log10(0.5) = 0.30103


### Autosomes ###
autopca <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/nps_pvalues_no_na_autosomes", header = T)
autopca$CHR <- 'Autosomes'
head(autopca)

median(autopca$data_pcadapt.pvalues)
#-log10(0.5) = 0.30103

pca_data <- rbind(autopca, Xpca)
pca_data_b <- separate(pca_data, locNames.data3., into = c("CHROM", "POS"), sep = "(?<=[0-9])_")
head(pca_data_b)

Impact_auto <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.NARW.delta_load_unmelted", header = T)
Impact_X <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/X.NARW.delta_load_unmelted", header = T)
Impact <- rbind(Impact_auto, Impact_X)
names(Impact)[names(Impact) == 'CHR'] <- 'CHROM'
head(Impact)

pca_data_D <- merge(Impact, pca_data_b,by=c("CHROM","POS"))

P#CA_outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/PCAdapt.refined.list", header = T)
#PCA_outliers$Outlier <- "candidate"

outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/candidate.variants", header = T)
outliers$Outlier <- "candidate"

NARW_scanPCA <- merge(pca_data_D, PCA_outliers, all = T)
NARW_scanPCA$Outlier[is.na(NARW_scanPCA$Outlier)] <- "non"
Chrom_no <- subset(NARW_scanLoad, select = c(14,15,17))
NARW_scanPCA_b <- merge(NARW_scanPCA, Chrom_no, by = c("CHROM", "POS"))

### only keep candidates with p values below threshold ... multiallelic sites get duplicated while merging  
NARW_scanPCA_c <- NARW_scanPCA_b[which(NARW_scanPCA_b$Outlier == "candidate" & NARW_scanPCA_b$data_pcadapt.pvalues <= 0.04657388 |
                                         NARW_scanPCA_b$Outlier == "non"),] 



NARW_scanPCA_c$condition <- case_when(
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "Low" & NARW_scanPCA_c$Outlier =="candidate" ~ "Low",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "Low" & NARW_scanPCA_c$Outlier =="candidate" ~ "Low",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "Moderate" & NARW_scanPCA_c$Outlier =="candidate" ~ "Moderate",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "Moderate" & NARW_scanPCA_c$Outlier =="candidate" ~ "Moderate",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "High" & NARW_scanPCA_c$Outlier == "candidate" ~ "High",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "High" & NARW_scanPCA_c$Outlier == "candidate" ~ "High",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 ~ "odd",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 ~ "even",
  TRUE ~ "unknown"
)


### violin or density plot ###
pca_lines <- data.frame(value=c(0.30103,0.30103),boxplot.nr=c("Autosomes", "X chromosome"))
pca_threshold <- data.frame(value=1.331858)


FourA <- ggplot(NARW_scanPCA_c, aes(x=CHR, y= -log10(data_pcadapt.pvalues), fill = Impact)) + 
  geom_violin(scale = "width") + theme_classic() +
  stat_summary(fun = "median",geom = "point", color = "black", position = position_dodge(.9)) +
  scale_fill_manual(values=c("gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  #geom_boxplot(data=pca_lines,aes(factor(boxplot.nr),value),
               #inherit.aes=FALSE,color="black",size=.25, linetype = "solid") +
  geom_hline(data=pca_threshold, linetype="dashed", size=.25, aes(yintercept=value)) +
  labs(y = expression("-log"["10"] ~ "(P value)")) +
  #ylab("-log (P)") +
  xlab("Genomic region") +
  theme(legend.position = "none")






##########################################################
#########   Plot genome scans with outliers     ##########
##########################################################

### read in and organize Fst and Mutation Load data ###

## autosomes Fst ##
autosomes.scan.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/autosomes.fst", sep = '\t', header = T)

## X Fst ##
X.scan.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/X.fst", sep = '\t', header = T)

NARW.scan.fst <- rbind(autosomes.scan.fst, X.scan.fst)
NARW.scan.fst <- NARW.scan.fst[which(NARW.scan.fst$Impact != "Impact"),]
head(NARW.scan.fst)
nrow(NARW.scan.fst)

### Total mutation load ###

Autosomes.NARW.homozygosity  <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.NARW.delta_load_unmelted", sep = '\t', header = T)
Autosomes.NARW.homozygosity$CHR <- 'Autosomes'
X.NARW.homozygosity  <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/X.NARW.delta_load_unmelted", sep = '\t', header = T)
X.NARW.homozygosity$CHR <- 'X chromosome'
NARW.homozygosity <- rbind(Autosomes.NARW.homozygosity, X.NARW.homozygosity)
head(NARW.homozygosity)
nrow(NARW.homozygosity)


#### merge Fst and mutation load ... ###
NARW_scan <- cbind(NARW.homozygosity, NARW.scan.fst)
NARW_scan <- subset(NARW_scan, select = -c(2))

## number chromosomes and positions ##
NARW_scan$CHROM_no <- NARW_scan %>% group_indices(CHROM) 
NARW_scan$POS <- as.numeric(as.character(NARW_scan$POS))
NARW_scan$WEIR_AND_COCKERHAM_FST <- as.numeric(as.character(NARW_scan$WEIR_AND_COCKERHAM_FST))
nrow(NARW_scan)

## Mutation load 
NARW_scanLoad <- NARW_scan[which(NARW_scan$Impact != "Modifier"),]
head(NARW_scanLoad)
summary(NARW_scanLoad)

## Fst
NARW_scan <- NARW_scan[which(NARW_scan$Impact != "Impact"),]
NARW_scanFst <- NARW_scan[which(NARW_scan$Impact != "Modifier"),]
nrow(NARW_scanFst)
hist(NARW_scanFst$WEIR_AND_COCKERHAM_FST)

NARW_scanFst$WEIR_AND_COCKERHAM_FST <-ifelse(NARW_scanFst$WEIR_AND_COCKERHAM_FST < 0 ,0 ,NARW_scanFst$WEIR_AND_COCKERHAM_FST)
hist(NARW_scanFst$WEIR_AND_COCKERHAM_FST)
nrow(NARW_scanFst)


#PCAdapt
#PCA_outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/PCAdapt.refined.list", header = T)
#PCA_outliers$Outlier <- "candidate"

outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/candidate.variants", header = T)
outliers$Outlier <- "candidate"



NARW_scanPCA <- merge(pca_data_D, outliers, all = T)
NARW_scanPCA$Outlier[is.na(NARW_scanPCA$Outlier)] <- "non"
Chrom_no <- subset(NARW_scanLoad, select = c(14,15,17))
NARW_scanPCA_b <- merge(NARW_scanPCA, Chrom_no, by = c("CHROM", "POS"))


### Plot scan for PCAdapt ###
### only keep candidates with p values below threshold ... 
NARW_scanPCA_c <- NARW_scanPCA_b[which(NARW_scanPCA_b$Outlier == "candidate" & NARW_scanPCA_b$data_pcadapt.pvalues <= 0.04657388 |
                                NARW_scanPCA_b$Outlier == "non" & NARW_scanPCA_b$data_pcadapt.pvalues >= 0.000000000005),] 


NARW_scanPCA_c$condition <- case_when(
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "Low" & NARW_scanPCA_c$Outlier =="candidate" ~ "Low",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "Low" & NARW_scanPCA_c$Outlier =="candidate" ~ "Low",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "Moderate" & NARW_scanPCA_c$Outlier =="candidate" ~ "Moderate",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "Moderate" & NARW_scanPCA_c$Outlier =="candidate" ~ "Moderate",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 & NARW_scanPCA_c$Impact == "High" & NARW_scanPCA_c$Outlier == "candidate" ~ "High",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 & NARW_scanPCA_c$Impact == "High" & NARW_scanPCA_c$Outlier == "candidate" ~ "High",
  NARW_scanPCA_c$CHROM_no %% 2 != 0 ~ "odd",
  NARW_scanPCA_c$CHROM_no %% 2 == 0 ~ "even",
  TRUE ~ "unknown"
)



unique(NARW_scanPCA_c$condition)
#"darkred", 
COLS = c("gray50", "gray75", "orange3", "gold")
avalues= c(.35, .35, .99, .99)
fillCOLS = sapply(1:5,function(i)alpha(COLS[i],avalues[i]))

NARW_scanPCA_c$condition <- factor(NARW_scanPCA_c$condition , levels=c("odd", "even", "Moderate", "Low"))
ordered_NARW_scanPCA_c <- NARW_scanPCA_c[order(NARW_scanPCA_c$condition),]


FourD <- ggplot(ordered_NARW_scanPCA_c, aes(POS, -log10(data_pcadapt.pvalues), color = condition)) +
  #geom_point(size = 2.5) +
  geom_jitter(height = .05, width = .1, size = 3.5) +
  scale_color_manual("category", values = fillCOLS) + 
  facet_grid(. ~ CHROM_no, scales = "free_x", switch = "x")  +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  theme_classic() +
  geom_hline(data=pca_lines, linetype="solid", size=.25, aes(yintercept=value)) +
  geom_hline(data=pca_threshold, linetype="dashed", size=.25, aes(yintercept=value)) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  labs(y = expression("-log"["10"] ~ "(P value)")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab("Chromosome")




  
## Fst 
  
outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/candidate.variants", header = T)
outliers$Outlier <- "candidate"  
  
#Fst_outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Fst.refined.list", header = T)
#Fst_outliers$Outlier <- "candidate"

NARW_scanFst_b <- merge(NARW_scanFst, outliers, all = T)
nrow(NARW_scanFst_b)
NARW_scanFst_b$Outlier[is.na(NARW_scanFst_b$Outlier)] <- "non"

NARW_scanFst_b$condition <- case_when(
    NARW_scanFst_b$CHROM_no %% 2 != 0 & NARW_scanFst_b$Impact.1 == "Low" & NARW_scanFst_b$Outlier =="candidate" ~ "Low",
    NARW_scanFst_b$CHROM_no %% 2 == 0 & NARW_scanFst_b$Impact.1 == "Low" & NARW_scanFst_b$Outlier =="candidate" ~ "Low",
    NARW_scanFst_b$CHROM_no %% 2 != 0 & NARW_scanFst_b$Impact.1 == "Moderate" & NARW_scanFst_b$Outlier =="candidate" ~ "Moderate",
    NARW_scanFst_b$CHROM_no %% 2 == 0 & NARW_scanFst_b$Impact.1 == "Moderate" & NARW_scanFst_b$Outlier =="candidate" ~ "Moderate",
    NARW_scanFst_b$CHROM_no %% 2 != 0 & NARW_scanFst_b$Impact.1 == "High" & NARW_scanFst_b$Outlier == "candidate" ~ "High",
    NARW_scanFst_b$CHROM_no %% 2 == 0 & NARW_scanFst_b$Impact.1 == "High" & NARW_scanFst_b$Outlier == "candidate" ~ "High",
    NARW_scanFst_b$CHROM_no %% 2 != 0 ~ "odd",
    NARW_scanFst_b$CHROM_no %% 2 == 0 ~ "even",
    TRUE ~ "unknown"
  )
 
NARW_scanFst_c <- NARW_scanFst_b[which(NARW_scanFst_b$Outlier == "candidate" & NARW_scanFst_b$WEIR_AND_COCKERHAM_FST >= 0.25 |
                                         NARW_scanFst_b$Outlier == "non"),] 

test <- NARW_scanFst_c[which(NARW_scanFst_c$Outlier == "candidate"),] 
nrow(test)

unique(NARW_scanFst_c$condition)
#"darkred"
COLS = c("gray50","gray75", "orange3", "gold")
avalues= c(.35, .35, .99, .99)
fillCOLS = sapply(1:5,function(i)alpha(COLS[i],avalues[i]))
 
#"High"
NARW_scanFst_c$condition <- factor(NARW_scanFst_c$condition , levels=c("odd", "even", "Moderate", "Low"))
ordered_NARW_scanFst_c <- NARW_scanFst_c[order(NARW_scanFst_c$condition),]
head(NARW_scanFst_c)

lines_auto <- data.frame(value=c(0.07),CHROM_no=c(1:21))
line_x <- data.frame(value=0.001,CHROM_no=(22))
scan_lines <- rbind(lines_auto, line_x) 
Fst_threshold <- data.frame(value=0.25)
nrow(ordered_NARW_scanFst_c)
  
FourE <- ggplot(ordered_NARW_scanFst_c, aes(POS, WEIR_AND_COCKERHAM_FST, color = condition)) +
    #geom_point(size = 2.5) +
    geom_jitter(height = .02, width = .1, size = 3.5) +
    scale_color_manual("category", values = fillCOLS) + 
    facet_grid(. ~ CHROM_no, scales = "free_x", switch = "x")  +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    theme_classic() +
    geom_hline(data=scan_lines, linetype="solid", size=.25, aes(yintercept=value)) +
    geom_hline(data=Fst_threshold, linetype="dashed", size=.25, aes(yintercept=value)) +
    theme(axis.text = element_text(size = 10, color = "black")) +
    theme(axis.title = element_text(size = 12, color = "black")) +
    ylab(expression(paste(italic("Fst")))) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
    xlab("Chromosome") +
    ylim(0,1)

  
# Mutation load #
  
  
outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/candidate.variants", header = T)
outliers$Outlier <- "candidate" 

#ML_outliers <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/mutationLoad.refined.list", header = T)
#ML_outliers$Outlier <- "candidate"

NARW_scanLoad_b <- merge(NARW_scanLoad, outliers, all = T) 
NARW_scanLoad_b$Outlier[is.na(NARW_scanLoad_b$Outlier)] <- "non"
head(NARW_scanLoad_b)
nrow(NARW_scanLoad_b)


NARW_scanLoad_b$condition <- case_when(
  NARW_scanLoad_b$CHROM_no %% 2 != 0 & NARW_scanLoad_b$Impact == "Low" & NARW_scanLoad_b$Outlier =="candidate" ~ "Low",
  NARW_scanLoad_b$CHROM_no %% 2 == 0 & NARW_scanLoad_b$Impact == "Low" & NARW_scanLoad_b$Outlier =="candidate" ~ "Low",
  NARW_scanLoad_b$CHROM_no %% 2 != 0 & NARW_scanLoad_b$Impact == "Moderate" & NARW_scanLoad_b$Outlier =="candidate" ~ "Moderate",
  NARW_scanLoad_b$CHROM_no %% 2 == 0 & NARW_scanLoad_b$Impact == "Moderate" & NARW_scanLoad_b$Outlier =="candidate" ~ "Moderate",
  NARW_scanLoad_b$CHROM_no %% 2 != 0 & NARW_scanLoad_b$Impact == "High" & NARW_scanLoad_b$Outlier == "candidate" ~ "High",
  NARW_scanLoad_b$CHROM_no %% 2 == 0 & NARW_scanLoad_b$Impact == "High" & NARW_scanLoad_b$Outlier == "candidate" ~ "High",
  NARW_scanLoad_b$CHROM_no %% 2 != 0 ~ "odd",
  NARW_scanLoad_b$CHROM_no %% 2 == 0 ~ "even",
  TRUE ~ "unknown"
)


NARW_scanLoad_c <- NARW_scanLoad_b[which(NARW_scanLoad_b$Outlier == "candidate" & NARW_scanLoad_b$deltaGenetic >= 0.06 |
                                         NARW_scanLoad_b$Outlier == "non"),] 



unique(NARW_scanLoad_c$condition)
#"darkred"
COLS = c("gray50", "gray75", "orange3", "gold")
avalues= c(.35, .35, .99, .99)
fillCOLS = sapply(1:5,function(i)alpha(COLS[i],avalues[i]))

#"High"
NARW_scanLoad_c$condition <- factor(NARW_scanLoad_c$condition , levels=c("odd", "even", "Moderate", "Low"))
ordered_NARW_scanLoad_c <- NARW_scanLoad_c[order(NARW_scanLoad_c$condition),]
nrow(ordered_NARW_scanLoad_c)

FourF <- ggplot(ordered_NARW_scanLoad_c, aes(POS, deltaGenetic, color = condition)) +
  #geom_point(size = 2.5) +
  geom_jitter(height = .02, width = .1, size = 3.5) +
  scale_color_manual("category", values = fillCOLS) + 
  facet_grid(. ~ CHROM_no, scales = "free_x", switch = "x")  +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  theme_classic() +
  geom_hline(yintercept=0 , linetype="solid", color = "black", size=.25) +
  geom_hline(yintercept=0.06 , linetype="dashed", color = "black", size=.25) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  ylab(expression(Delta~ "Total mutation load")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab("Chromosome")




Figure_4 <- (
  (FourA/FourB/FourC) | (FourD/FourE/FourF)) + 
  plot_layout(ncol=2,widths=c(2,9)) + 
  plot_annotation(tag_levels = 'A') 

Figure_4 #view multi-panel figure








######################################################################
#                                                                    #
#    gene set enrichment with Getsalt (or whatever its called)       #
#                 Total mother/offspring burden                      #
######################################################################

######### plotting results #########3


geneontology <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1694542499/weighted_set.txt", 
                     header = T,
                     sep="\t")

geneontology
outlier1 <- geneontology %>%
  mutate(V2 = fct_reorder(description, enrichmentRatio)) %>%
  ggplot(aes(y=V2, x=enrichmentRatio, fill = "goldenrod3")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("goldenrod3")) +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.y = element_text(size = 12, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylab("GO term description (geneontology)") +
  #xlab("Biological process") +
  theme(legend.position = "none")





KEGG_pathway <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1694542398/weighted_set.txt", 
                     header = T,
                     sep="\t")

KEGG_pathway
outlier2 <- KEGG_pathway %>%
  mutate(V2 = fct_reorder(description, enrichmentRatio)) %>%
  ggplot(aes(y=V2, x=enrichmentRatio, fill = "gold4")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("gold4")) +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.y = element_text(size = 12, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylab("GO term description (pathway)") +
  #xlab("Biological process") +
  theme(legend.position = "none") 



Biological_process <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1694542398/goslim_summary_wg_result1694542398_bp.txt", 
                     header = F,
                     sep="\t")
Biological_process
bp <- Biological_process %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "slateblue1")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("slateblue1")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12, color = "black")) +
  ylab("Count") +
  #xlab("Biological process") +
  #scale_y_continuous(breaks= pretty_breaks()) +
  theme(legend.position = "none")  



Cellular_component <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1694542398/goslim_summary_wg_result1694542398_cc.txt", 
                                 header = F,
                                 sep="\t")
Cellular_component
cc <- Cellular_component %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "slateblue3")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("slateblue3")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12, color = "black")) +
  ylab("Count") +
  #xlab("biological process") +
  #scale_y_continuous(breaks= pretty_breaks()) +
  theme(legend.position = "none")  



Molecular_function <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1694542398/goslim_summary_wg_result1694542398_mf.txt", 
                                 header = F,
                                 sep="\t")
Molecular_function
mf <- Molecular_function %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "slateblue4")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("slateblue4")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 10, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 12, color = "black")) +
  ylab("Count") +
  #xlab("Biological process") +
  #scale_y_continuous(breaks= pretty_breaks()) +
  theme(legend.position = "none")  




Figure_5 <- ((outlier1 / outlier2) | (bp / cc / mf)) + 
  #plot_layout(ncol=2,widths =c(5,6)) +
  plot_annotation(tag_levels = 'A') #add figure labels
Figure_5 #view multi-panel figure





