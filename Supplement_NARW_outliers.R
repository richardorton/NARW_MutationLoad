

####################################
#                                  #
#  Chapter one - outlier analyses  #
#   Started: May 15, 2023          #
#   by RW Orton                    #
#                                  #
####################################

library("qvalue")
library("OutFLANK")
library("ggplot2")
library(vcfR)
library(hierfstat)
library(adegenet)
library(radiator)




##### PCAdapt #######


###### X chromosome ################

X_PCA <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/pca_adapt_pop_map_X")
X_PCA

ggplot(X_PCA, aes(x=PC1, y=PC2, color = POPULATION, shape = POPULATION))+
  geom_point(size=5) +
  scale_color_manual(values=c("black", "gray75")) +
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  #theme(legend.title = element_blank()) +
  #theme(legend.position = "none")  +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.2,aes(color=POPULATION)) 



###### Autosomes ################

autosomes_PCA <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/pca_adapt_pop_map_autosomes")

ggplot(autosomes_PCA, aes(x=PC1, y=PC2, color = POPULATION, shape = POPULATION))+
  geom_point(size=5) +
  scale_color_manual(values=c("black", "gray75")) +
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  #theme(legend.title = element_blank()) +
  #theme(legend.position = "none")  +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.2,aes(color=POPULATION)) 


### outliers ###

autosomes_snps <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/nps_pvalues_no_na_autosomes")
X_snps <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/nps_pvalues_no_na_X")
snps <- rbind(autosomes_snps, X_snps)
nrow(snps)

head(snps)
tail(snps)

### quantile approach ###
quantile(snps$data_pcadapt.pvalues, probs = c(.1, 0.9))
pcadapt_outliers_autosomes <- snps[which(snps$data_pcadapt.pvalues <= 0.04657388),] 
pcadapt_outliers_autosomes
nrow(pcadapt_outliers_autosomes)

write.table(pcadapt_outliers_autosomes, "/Volumes/cetacea/Genetic_load/_enrichment/PCAdapt.rawGene.list")


### FDR approach ###
qval <- qvalue(autosomes_snps$data_pcadapt.pvalues)$qvalues
FDR <- cbind(autosomes_snps,qval)
head(FDR)
alpha <- 0.15
outliers_FDRautosomes <- FDR[which(FDR$qval <= alpha),] 
head(outliers_FDRautosomes)

outliers_FDRautosomes


write.table(outliers_FDRautosomes, "/Volumes/cetacea/Genetic_load/_enrichment/PCAdapt.rawGene.autosomeslist")



### Hample filter ###

lower_bound <- median(outliers_FDRautosomes$qval) - 3 * mad(outliers_FDRautosomes$qval, constant = 1)
lower_bound


outliers_FDR <- FDR[which(FDR$qval <= 0.005780537),] 
outliers_FDR



#########################################################################
#########      Determine oultiers per site for Fst and ML      ##########
#########################################################################


## autosomes Fst ##
autosomes.scan.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/autosomes.fst", sep = '\t', header = T)

## X Fst ##
X.scan.fst <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/differentiation/X.fst", sep = '\t', header = T)

NARW.scan.fst <- rbind(autosomes.scan.fst, X.scan.fst)
NARW.scan.fst <- NARW.scan.fst[which(NARW.scan.fst$Impact != "Impact"),]
head(NARW.scan.fst)


### Total mutation load ###

Autosomes.NARW.homozygosity  <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.NARW.delta_load_unmelted", sep = '\t', header = T)
Autosomes.NARW.homozygosity$CHR <- 'Autosomes'
X.NARW.homozygosity  <- read.csv("/Volumes/cetacea/Genetic_load/NARW_females/X.NARW.delta_load_unmelted", sep = '\t', header = T)
X.NARW.homozygosity$CHR <- 'X chromosome'
NARW.homozygosity <- rbind(Autosomes.NARW.homozygosity, X.NARW.homozygosity)
head(NARW.homozygosity)


#### merge Fst and mutation load ... ###
NARW_scan <- cbind(NARW.homozygosity, NARW.scan.fst)
NARW_scan <- subset(NARW_scan, select = -c(2))
head(NARW_scan)

## number chromosomes and positions ##
NARW_scan$CHROM_no <- NARW_scan %>% group_indices(CHROM) 
NARW_scan$POS <- as.numeric(as.character(NARW_scan$POS))
NARW_scan$WEIR_AND_COCKERHAM_FST <- as.numeric(as.character(NARW_scan$WEIR_AND_COCKERHAM_FST))
str(NARW_scan)
head(NARW_scan)

NARW_scan <- NARW_scan[which(NARW_scan$Impact != "Impact"),]
NARW_scanFst <- NARW_scan[which(NARW_scan$WEIR_AND_COCKERHAM_FST >= 0),]
NARW_scanFst <- NARW_scanFst[which(NARW_scanFst$Impact != "Modifier"),]
NARW_scanLoad <- NARW_scan[which(NARW_scan$deltaGenetic > 0),]

head(NARW_scanLoad)

## is the burden on the calf or mother? ##################

######## Select outliers REALIZED load ###################
######## Select outliers REALIZED load ###################
######## Select outliers REALIZED load ###################


#head(NARW_scan)
#unique(NARW_scan$Impact.1)
#of_interest <- NARW_scan[which(NARW_scan$lowFecund_Alt_GF > .7 &
#NARW_scan$highFecund_GF < .3 &
#NARW_scan$deltaHom >= 4),]

fixed <- NARW_scanLoad[which(NARW_scanLoad$deltaHom >= 5),]
fixed
## two variants are fixed - listed in supplemental table ##



######## Select outliers TOTAL load ###################
######## Select outliers TOTAL load ###################
######## Select outliers TOTAL load ###################


head(NARW_scanLoad)
head(NARW_scanFst)
tail(NARW_scanLoad)
tail(NARW_scanFst)
head(pca_data_D)
tail(head(pca_data_D))



unique(NARW_scan$Impact)
unique(NARW_scanLoad$Impact)
unique(NARW_scanFst$Impact)
unique(pca_data_D$Impact)

### Quantile filter genetic load ###
quantile(NARW_scanLoad$deltaGenetic, probs = .9)
quantile(NARW_scanFst$WEIR_AND_COCKERHAM_FST, probs = .9)

### Hample filter genetic load ###
upper_bound <- median(NARW_scanLoad$deltaGenetic) + 3 * mad(NARW_scanLoad$deltaGenetic, constant = 1)
upper_bound

upper_bound <- median(NARW_scanFst$WEIR_AND_COCKERHAM_FST) + 3 * mad(NARW_scanFst$WEIR_AND_COCKERHAM_FST, constant = 1)
upper_bound


### 95% and Hample filter are the same ##

of_interest_geneticLoad <- NARW_scanLoad[which(NARW_scanLoad$deltaGenetic >= .06),]
of_interest_geneticLoad
write.table(of_interest_geneticLoad, "/Volumes/cetacea/Genetic_load/_enrichment/mutationLoad.rawGene.list")

of_interest_Fst <- NARW_scanFst[which(NARW_scanFst$WEIR_AND_COCKERHAM_FST >= 0.25),]
of_interest_Fst
nrow(of_interest_Fst)
write.table(of_interest_Fst, "/Volumes/cetacea/Genetic_load/_enrichment/Fst.rawGene.list")


data <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/total/raw_gene_interest.totalLoad.bed", header = T)
data$STOP <- data$WEIR_AND_COCKERHAM_FST + 1
write.table(data, "/Volumes/cetacea/Genetic_load/_enrichment/total/totalLoad.genelist.bed")









