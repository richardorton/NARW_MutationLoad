

####################################
#                                  #
#   Chapter one - genetic load     #
#   Started: May 15, 2023          #
#   by RW Orton                    #
#                                  #
####################################

library(devtools)
library(ggplot2)
library(detectRUNS)
library(dplyr)
library(tidyverse)
library(tidyr)
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
library(ggridges)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(devtools)
library(treeio)
library(ggtree)
library(scales)
library(plotrix)




########################################################################################
#                                       Color scheme                                   #
########################################################################################





########################################################################
#                                                                      #                                                                   
#                         Interspecific analyses                       #
#                                                                      #
########################################################################





########################################################################################
#                               Whole mtDNA tree with BEAST JC69                       #
########################################################################################



tree <- read.beast("/Users/richardorton/Desktop/mt_genomes.fasta/mcc_mtGenomes_BEAST.tree")
tree


# get branch node numbers
ggtree(tree, branch.length="none") + geom_text(aes(label=node))
tree

#scale_fill_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  
##### tree no. 1 #######
#Onea <- 
  
  ggtree(tree) + geom_tree(cex = .75) + theme_tree2() +  
  #geom_tiplab() + + 
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + 
  geom_hilight(node=59, fill="#8EB6B5", alpha = .75) + 
  geom_hilight(node=50, fill="aquamarine4", alpha = .75) +
  geom_hilight(node=39, fill = "darkslategray", alpha = .75) +
  theme(axis.text.x = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) + 
  xlab("Substitution rate") +
  geom_point2(aes(label=round(as.numeric(posterior),2),
                  subset=as.numeric(posterior)> 0.95), size = 2, color = "black") +
  geom_point2(aes(label=round(as.numeric(posterior),2),
                  subset=as.numeric(posterior) < 0.95), size = 2, shape = 18)


  #geom_hilight(node=59, fill="#8EB6B5", alpha = .75) + 
  #geom_hilight(node=50, fill="#D6C350", alpha = .75) +
  #geom_hilight(node=39, fill = "violetred4", alpha = .75) +
  
  


###### tree no.2 #######
ggtree(tree, layout="circular", branch.length = 'none') +
  #geom_tiplab() + 
  geom_hilight(node= 39, fill="darkslateblue") + 
  geom_hilight(node= 59, fill="lightskyblue4") +
  geom_hilight(node= 50, fill = "deepskyblue4")


###### tree no.3 #######
ggtree(tree, layout="daylight", branch.length = 'none') +
  #geom_tiplab() + 
  geom_hilight(node=39, fill="darkslateblue") + 
  geom_hilight(node=59, fill="cadetblue4") +
  geom_hilight(node= 50, fill = "palevioletred4")





########################################################################################
#                       Nucleotide diversity - species wide                            #
########################################################################################

## read in autosome data (100Mb windows)
int_autosomes_diversity <- read.table("/Volumes/cetacea/Genetic_load/pixy_output/autosomes_out/pixy_pi.txt",sep="\t", header = T)
head(int_autosomes_diversity)
int_autosomes_diversity$CHR <- 'Autosomes' 



## read in X data (100Mb windows)
int_X_diversity <- read.table("/Volumes/cetacea/Genetic_load/pixy_output/X_out/pixy_pi.txt",sep="\t", header = T)
head(int_X_diversity)
int_X_diversity$CHR <- 'X chromosome' 

interspecific_pi <- rbind(int_autosomes_diversity, int_X_diversity)
head(interspecific_pi)

write.table(interspecific_pi, "/Volumes/cetacea/Genetic_load/pixy_output/genome.pi", sep="\t", quote = F)

## box plot ##
interspecific_pi$pop <- factor(interspecific_pi$pop, levels=c("NARW", "SRW", "BH"))
Oneb <- ggplot(interspecific_pi, aes(x=CHR, y=avg_pi, fill = pop)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  xlab("Genomic region") +
  ylab(expression("Nucleotide diveristy"~(pi))) +
  ylim(0, 0.010) +
  theme(legend.position = "none")  


pi.aov2 <- aov(avg_pi  ~ pop + CHR, data = interspecific_pi)
summary(pi.aov2)

#Df  Sum Sq  Mean Sq F value Pr(>F)    
#pop             2 0.02614 0.013071   13762 <2e-16 ***
#  CHR             1 0.00230 0.002303    2425 <2e-16 ***
#  Residuals   70799 0.06724 0.000001  


### Get point estimates of Ne using Pi ###
head(interspecific_pi)
str(interspecific_pi)
### calculate standard error and mean for each group ###
se <- function(x) sd(x) / sqrt(length(x))

interspecific_pi <- na.omit(interspecific_pi)
interspecific_pi$avg_pi <- as.numeric(interspecific_pi$avg_pi)

Pi_summary <- interspecific_pi %>% 
  group_by(pop, CHR) %>% 
  summarise(mean = mean(avg_pi),
            SE = se(avg_pi),
            SD = sd(avg_pi),
            Ne = )

Pi_summary$Ne <- Pi_summary$mean/(4*.00000002)

Pi_summary
#pop   CHR           mean         SE       SD     Ne
#<fct> <chr>        <dbl>      <dbl>    <dbl>  <dbl>
#  1 BH    Autosomes 0.00171  0.00000824 0.00123  21419.
#3 NARW  Autosomes 0.000301 0.00000365 0.000546  3764.
#5 SRW   Autosomes 0.00155  0.00000710 0.00106  19375.

#2 BH    X         0.000546 0.0000162  0.000567  6826.
#4 NARW  X         0.000116 0.00000788 0.000277  1449.
#6 SRW   X         0.000470 0.0000116  0.000406  5878.

Pi_summary$Ne <- Pi_summary$mean/(4*(2*10^-8))

NARW_BH_auto <- 21419-3764
NARW_SRW_auto <- 19375-3764
SRW_BH_auto <- 21419-19375

NARW_BH_X <- 6826-1449
NARW_SRW_X <- 5878-1449
SRW_BH_X <- 6826-5878

NARW_BH_auto
NARW_SRW_auto 
SRW_BH_auto
NARW_BH_X
NARW_SRW_X 
SRW_BH_X 


########################################################################################
#                                   F - species wide                                   #
########################################################################################


## read in (individual) homozygote data

# autosomes
autosomes_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.ALL.het", header = T)
head(autosomes_F_inter_data)
autosomes_F_inter_data$CHR <- 'Autosomes'

autosomes_F_inter_data$pop <- c('NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW',
                                   'SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','NARW','NARW',
                                   'BH','BH','BH','BH','BH','BH','BH','BH','BH','BH')


# X
X_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.ALL.het", header = T)
head(X_F_inter_data)
X_F_inter_data$CHR <- 'X chromosome'

X_F_inter_data$pop <- c('NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW',
                                'SRW','SRW','SRW','SRW','SRW','NARW','NARW',
                                'BH','BH','BH')


F_inter_data <- rbind(autosomes_F_inter_data, X_F_inter_data)
head(F_inter_data)

write.table(F_inter_data, "/Volumes/cetacea/Genetic_load/interspecific/ALL.het", sep="\t", quote = F)


## box plot ##
F_inter_data$pop <- factor(F_inter_data$pop, levels=c("NARW", "SRW", "BH"))
Onec <- ggplot(F_inter_data, aes(x=CHR, y=F, fill = pop)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  theme(axis.title.y = element_text(face = "italic")) +
  xlab("Genomic region") +
  theme(legend.position = "none") + 
  ylab(expression(paste("Inbreeding coefficient (", italic("F"),")"))) +
  ylim(0,1.0)

#  labs(y=expression(paste("Inbreeding coefficient (",italic("F")")")))
#   ylab("Inbreding coeffcient (F)") +

####### plot heterozygosity #########################################################

F_inter_data$Het <- (F_inter_data$N_SITES - F_inter_data$O.HOM.)/F_inter_data$N_SITES
head(F_inter_data)

F_inter_data

F_inter_data$pop <- factor(F_inter_data$pop, levels=c("NARW", "SRW", "BH"))
Onec <- ggplot(F_inter_data, aes(x=CHR, y=Het, fill = pop)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Heterozygosity") +
  xlab("Genomic region") #+
  theme(legend.position = "none")  



Het.aov2 <- aov(Het  ~ pop + CHR, data = F_inter_data)
summary(Het.aov2 )

#Df  Sum Sq Mean Sq F value   Pr(>F)    
#pop          2 0.12482 0.06241  487.17  < 2e-16 ***
#  CHR          1 0.00255 0.00255   19.87 5.28e-05 ***
#  Residuals   46 0.00589 0.00013      







########################################################################################
#                                  ROH distributions                                   #
########################################################################################

autosomes_roh <- read.table("/Volumes/cetacea/Genetic_load/_roh/autosomes.population.roh", header = T)
autosomes_roh$CHR <- 'Autosomes'
head(autosomes_roh)

automsomes_nROH <- count(autosomes_roh, "V2")
automsomes_SROH <- aggregate(autosomes_roh$V6, by=list(V2=autosomes_roh$V2), FUN=sum)
automsomes_meanLROH <- aggregate(autosomes_roh$V6, by=list(V2=autosomes_roh$V2), FUN=mean)

autosome_roh_stats <- cbind(automsomes_nROH, automsomes_SROH, automsomes_meanLROH)
autosome_roh_stats$CHR <- 'Autosomes'
head(autosome_roh_stats)

## plot autosomes ROH distribution ##
autosomes_roh$Population <- factor(autosomes_roh$Population , levels=c("BH", "SRW", "NARW"))

Oned <- ggplot(autosomes_roh, aes(x = V6, y = Population, fill = Population)) + 
  geom_density_ridges() + 
  stat_density_ridges(quantile_lines = TRUE) +
  scale_fill_manual(values=c("#8EB6B5", "aquamarine4", "darkslategray")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45)) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Species") +
  scale_x_continuous(breaks=seq(0, 2000000, 500000),
                     labels=c("0", "0.5", "1.0", "1.5", "2.0")) +
  #scale_x_continuous(labels=scientific) + 
  #xlim("0", "2000000") +
  xlab("Autosomes - Length of ROH (Mb)") +
  theme(legend.position = "none")  


## RoH binned counts ###
auto1Mb <- subset(autosomes_roh, V6 >= 1000000)
auto_nROH_1Mb <- count(auto1Mb, "V2")
auto_nROH_1Mb

auto5Mb <- subset(autosomes_roh, V6 >= 5000000)
auto_nROH_5Mb <- count(auto5Mb, "V2")
auto_nROH_5Mb


auto100Kb <- subset(autosomes_roh, V6 <= 100000)
auto_nROH_1Kb <- count(auto100Kb, "V2")
auto_nROH_1Kb


#### X chromosome ###
X_roh <- read.table("/Volumes/cetacea/Genetic_load/_roh/X.population.roh", header = T)
X_roh$CHR <- 'X chromosome'
head(X_roh)

ROH_inter_data <- rbind(autosomes_roh, X_roh)
head(ROH_inter_data)
tail(ROH_inter_data)

write.table(ROH_inter_data, "/Volumes/cetacea/Genetic_load/interspecific/ALL.ROH", sep="\t", quote = F)


X_nROH <- count(X_roh, "V2")
X_SROH <- aggregate(X_roh$V6, by=list(V2=X_roh$V2), FUN=sum)
X_meanLROH <- aggregate(X_roh$V6, by=list(V2=X_roh$V2), FUN=mean)

X_roh_stats <- cbind(X_nROH, X_SROH, X_meanLROH)
X_roh_stats$CHR <- 'X Chromosome'

ALL.roh_stats <- rbind(autosome_roh_stats, X_roh_stats)
head(ALL.roh_stats)
tail(ALL.roh_stats)

write.table(ALL.roh_stats, "/Volumes/cetacea/Genetic_load/interspecific/ALL.ROH.stats", sep="\t", quote = F)


### plot X ROH distributions ###
X_roh$Population <- factor(X_roh$Population , levels=c("BH", "SRW", "NARW"))
X_roh <- X_roh[which(X_roh$V6 <= 1000000),] 

Onee <- ggplot(X_roh, aes(x = V6, y = Population, fill = Population)) + 
  geom_density_ridges() +
  scale_fill_manual(values=c("#8EB6B5", "aquamarine4", "darkslategray")) +
  theme_classic() +
  stat_density_ridges(quantile_lines = TRUE) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45)) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Species") +
  scale_x_continuous(breaks=seq(0, 1000000, 250000),
                     labels=c("0", "0.25", "0.5", "0.75", "1.0")) +
  xlab("X chromosome - Length of ROH (Mb)") +
  theme(legend.position = "none")  


x1Mb <- subset(X_roh, V6 >= 1000000)
x_nROH_1Mb <- count(x1Mb, "V2")
x_nROH_1Mb

x100Kb <- subset(X_roh, V6 <= 100000)
x_nROH_1Kb <- count(x100Kb, "V2")
x_nROH_1Kb



#### ROH statistics #####
ROH_total <- read.table("/Volumes/cetacea/Genetic_load/_roh/ROH_stats", header = T)
head(ROH_total)

autosomes_ROH_total <- ROH_total[which(ROH_total$CHR == "Autosomes"),] 

aggregate(autosomes_ROH_total$FROH, by=list(autosomes_ROH_total$Species), FUN=mean)
aggregate(autosomes_ROH_total$NROH, by=list(autosomes_ROH_total$Species), FUN=mean)
aggregate(autosomes_ROH_total$SROH, by=list(autosomes_ROH_total$Species), FUN=mean)
aggregate(autosomes_ROH_total$meanLROH, by=list(autosomes_ROH_total$Species), FUN=mean)




X_ROH_total <- ROH_total[which(ROH_total$CHR == "X_Chromosome"),] 

aggregate(X_ROH_total$FROH, by=list(X_ROH_total$Species), FUN=mean)
aggregate(X_ROH_total$NROH, by=list(X_ROH_total$Species), FUN=mean)
aggregate(X_ROH_total$SROH, by=list(X_ROH_total$Species), FUN=mean)
aggregate(X_ROH_total$meanLROH, by=list(X_ROH_total$Species), FUN=mean)


autosomes_ROH_total$Species <- factor(autosomes_ROH_total$Species , levels=c("NARW", "SRW", "BH"))
corr.a <- ggplot(autosomes_ROH_total, aes(x=SROH, y=NROH, color = Species)) + 
  geom_point(size = 7, alpha = 0.75) +
  theme_classic() + 
  theme(axis.text = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  scale_color_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5"))+
  geom_abline(intercept = 0, slope = 1/337718.77, linetype = "dashed", size = 1.5) +
  ylab("NROH") +
  xlab("SROH") +
  theme(legend.position = "none") +
  ylim(0,3000)



X_ROH_total$Species <- factor(X_ROH_total$Species , levels=c("NARW", "SRW", "BH"))
cor.x <- ggplot(X_ROH_total, aes(x=SROH, y=NROH, color = Species)) + 
  geom_point(size = 7, alpha = 0.75) +
  theme_classic() + 
  theme(axis.text = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  scale_color_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5"))+
  geom_abline(intercept = 0, slope = 1/111586.88, linetype = "dashed", size = 1.5) +
  ylab("NROH") +
  xlab("SROH") +
  theme(legend.position = "none") +
  ylim(0,1500)


roh_supp <- (corr.a | cor.x) +
  plot_annotation(tag_levels = 'A') #add figure labels
roh_supp #view multi-panel figure


########################################################################################
#                                 Figure Panel one ............                        #
########################################################################################


nested <- (
  
  Onea|
    (
      (Oneb|Onec)/
        (Oned/Onee)
      )) +
  plot_annotation(tag_levels = 'A') #add figure labels
nested #view multi-panel figure




########################################################################################
#                     Genetic load per individual - species wide                       #
########################################################################################


detach(package:plyr)
## get load data ##

autosomes_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.indv.inter.load.stats")
autosomes_geneticLoad_indv_data$CHR <- 'Autosomes'

X_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.indv.inter.load.stats")
X_geneticLoad_indv_data$CHR <- 'X chromosome'

inter_geneticLoad_indv_data <- rbind(autosomes_geneticLoad_indv_data, X_geneticLoad_indv_data)
tail(inter_geneticLoad_indv_data)

head(inter_geneticLoad_indv_data)

### set selection coefficient ###
inter_geneticLoad_indv_data$s <- case_when(
  inter_geneticLoad_indv_data$impact == "Modifier" ~ "0",
  inter_geneticLoad_indv_data$impact == "Low" ~ ".1",
  inter_geneticLoad_indv_data$impact == "Moderate" ~ ".3",
  inter_geneticLoad_indv_data$impact == "High" ~ ".6",
  TRUE ~ "unknown"
)


inter_geneticLoad_indv_data <- inter_geneticLoad_indv_data %>% mutate_at(c('s'), as.numeric)
inter_geneticLoad_indv_data_NS <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$impact != "Modifier"),]

head(inter_geneticLoad_indv_data_NS)

### calculate total genetic load ###
inter_geneticLoad_indv_data_NS$TotalLoad <- inter_geneticLoad_indv_data_NS$alt_hom_F*inter_geneticLoad_indv_data_NS$s + 
  0.5*inter_geneticLoad_indv_data_NS$s*inter_geneticLoad_indv_data_NS$het_F 

total.aov2 <- aov(TotalLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(total.aov2)

#Df  Sum Sq   Mean Sq F value Pr(>F)
#population    2 0.00061 0.0003058   0.192  0.826
#CHR           1 0.00087 0.0008746   0.548  0.460
#Residuals   146 0.23299 0.0015958 


### calculate realized genetic load ###
inter_geneticLoad_indv_data_NS$realizedLoad <- inter_geneticLoad_indv_data_NS$alt_hom_F*inter_geneticLoad_indv_data_NS$s 

realized.aov2 <- aov(realizedLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(realized.aov2)

#Df  Sum Sq  Mean Sq F value Pr(>F)  
#population    2 0.00713 0.003563   4.552 0.0121 *
#  CHR           1 0.00000 0.000002   0.003 0.9573  
#Residuals   146 0.11428 0.000783                 



### calculate masked genetic load ###
inter_geneticLoad_indv_data_NS$MaskedLoad <- inter_geneticLoad_indv_data_NS$het_F*inter_geneticLoad_indv_data_NS$s*0.5 

masked.aov2 <- aov(MaskedLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(masked.aov2)


#Df  Sum Sq  Mean Sq F value   Pr(>F)    
#population    2 0.00867 0.004334  19.951 2.19e-08 ***
#  CHR           1 0.00079 0.000788   3.628   0.0588 .  
#Residuals   146 0.03171 0.000217



head(inter_geneticLoad_indv_data_NS)
## aggregate to get genome-wide estimate of Total genetic load ##
agg_inter_geneticLoad_indv_data_NS <- aggregate(inter_geneticLoad_indv_data_NS$TotalLoad,  
                                                by=list(inter_geneticLoad_indv_data_NS$sample,
                                                        inter_geneticLoad_indv_data_NS$CHR,
                                                        inter_geneticLoad_indv_data_NS$impact,
                                                        inter_geneticLoad_indv_data_NS$population), FUN=median)




head(agg_inter_geneticLoad_indv_data_NS)

inter_geneticLoad_indv_data_NS_loads <- agg_inter_geneticLoad_indv_data_NS %>% spread(Group.3, x)
inter_geneticLoad_indv_data_NS_loads$sumLoad <- (inter_geneticLoad_indv_data_NS_loads$High + 
                                                   inter_geneticLoad_indv_data_NS_loads$Moderate +  
                                                   inter_geneticLoad_indv_data_NS_loads$Low)
inter_geneticLoad_indv_data_NS_loads



se <- function(x) sd(x) / sqrt(length(x))

sum_Totalload_summary <- inter_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

sum_Totalload_summary
###################################################################

head(inter_geneticLoad_indv_data_NS)
## aggregate to get genome-wide estimate of realized genetic load ##
agg_inter_geneticLoad_indv_data_NS <- aggregate(inter_geneticLoad_indv_data_NS$realizedLoad,  
                                                by=list(inter_geneticLoad_indv_data_NS$sample,
                                                        inter_geneticLoad_indv_data_NS$CHR,
                                                        inter_geneticLoad_indv_data_NS$impact,
                                                        inter_geneticLoad_indv_data_NS$population), FUN=median)




head(agg_inter_geneticLoad_indv_data_NS)
agg_inter_geneticLoad_indv_data_NS

inter_geneticLoad_indv_data_NS_loads <- agg_inter_geneticLoad_indv_data_NS %>% spread(Group.3, x)
inter_geneticLoad_indv_data_NS_loads$sumLoad <- (inter_geneticLoad_indv_data_NS_loads$High + 
                                                   inter_geneticLoad_indv_data_NS_loads$Moderate +  
                                                   inter_geneticLoad_indv_data_NS_loads$Low)




se <- function(x) sd(x) / sqrt(length(x))

sum_realizedload_summary <- inter_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

sum_realizedload_summary
inter_geneticLoad_indv_data_NS_loads
###################################################################

head(inter_geneticLoad_indv_data_NS)
## aggregate to get genome-wide estimate of masked genetic load ##
agg_inter_geneticLoad_indv_data_NS <- aggregate(inter_geneticLoad_indv_data_NS$MaskedLoad,  
                                                by=list(inter_geneticLoad_indv_data_NS$sample,
                                                        inter_geneticLoad_indv_data_NS$CHR,
                                                        inter_geneticLoad_indv_data_NS$impact,
                                                        inter_geneticLoad_indv_data_NS$population), FUN=median)




head(agg_inter_geneticLoad_indv_data_NS)

inter_geneticLoad_indv_data_NS_loads <- agg_inter_geneticLoad_indv_data_NS %>% spread(Group.3, x)
inter_geneticLoad_indv_data_NS_loads$sumLoad <- (inter_geneticLoad_indv_data_NS_loads$High + 
                                                   inter_geneticLoad_indv_data_NS_loads$Moderate +  
                                                   inter_geneticLoad_indv_data_NS_loads$Low)




se <- function(x) sd(x) / sqrt(length(x))

sum_maskedload_summary <- inter_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

sum_maskedload_summary


##################################### PLOT ###############################################################

sum_Totalload_summary$Group.4 <- factor(sum_Totalload_summary$Group.4, levels=c("NARW", "SRW", "BH"))
scaleFUN <- function(x) sprintf("%.2f", x)
  
a1 <- ggplot((sum_Totalload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  scale_y_continuous(labels=scaleFUN) +
  ylab("Total mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none") +
  ylim(0, 0.30)



sum_realizedload_summary$Group.4 <- factor(sum_realizedload_summary$Group.4, levels=c("NARW", "SRW", "BH"))
a2 <- ggplot((sum_realizedload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  ylab("Realized mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none") +
  ylim(0, 0.30)

sum_maskedload_summary$Group.4 <- factor(sum_maskedload_summary$Group.4, levels=c("NARW", "SRW", "BH"))
a3 <- ggplot((sum_maskedload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 16, color = "black")) +
  theme(axis.title = element_text(size = 20, color = "black")) +
  ylab("Masked mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  +
  ylim(0, 0.30)


Figure_2 <- (a1|a2|a3) +
  plot_annotation(tag_levels = 'A') #add figure labels
Figure_2 #view multi-panel figure






###################################################
#########    Genetic purging ######################
###################################################

#############################
#### Rxy (per site) #########
#############################

## read in site genetic load data ##

# autosomes
autosomes_inter_site_load_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.inter.site.load", header = T)
head(autosomes_inter_site_load_data)
autosomes_inter_site_load_data$Chromosome <- 'Autosomes'

# X
X_inter_site_load_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.inter.site.load", header = T)
head(X_inter_site_load_data)
X_inter_site_load_data$Chromosome <- 'X chromosome'


inter_site_load_data <- rbind(autosomes_inter_site_load_data, X_inter_site_load_data)
head(inter_site_load_data)
tail(inter_site_load_data)



agg_Rxy_data <- aggregate(inter_site_load_data$alternate_allele_frequency,  
                             by=list(inter_site_load_data$Population,
                                     inter_site_load_data$Chromosome,
                                     inter_site_load_data$Impact), 
                          FUN = function(x) c(mean = mean(x), se = std.error(x)))

agg_Rxy_data <- aggregate(inter_site_load_data$alternate_allele_frequency,  
                          by=list(inter_site_load_data$Population,
                                  inter_site_load_data$Chromosome,
                                  inter_site_load_data$Impact), 
                          FUN = mean)


agg_Rxy_data



agg_Rxy_data_B <- agg_Rxy_data %>% spread(Group.1, x)
agg_Rxy_data_B$NARW_SRW <- agg_Rxy_data_B$NARW/agg_Rxy_data_B$SRW
agg_Rxy_data_B$NARW_BH <- agg_Rxy_data_B$NARW/agg_Rxy_data_B$BH
agg_Rxy_data_B$SRW_BH <- agg_Rxy_data_B$SRW/agg_Rxy_data_B$BH


agg_Rxy_data_B

agg_Rxy_data_C <- agg_Rxy_data_B %>% gather(Group.1, x)
agg_Rxy_data_C <- melt(agg_Rxy_data_B, id.var = c('Group.2','Group.3'), variable.name = 'estimate')

head(agg_Rxy_data_C)
agg_Rxy_data_C


agg_Rxy_data_C$comparison <- case_when(
  agg_Rxy_data_C$estimate == "NARW_SRW" ~ "NARW:SRW",
  agg_Rxy_data_C$estimate == "NARW_BH" ~ "NARW:BH",
  agg_Rxy_data_C$estimate == "SRW_BH" ~ "SRW:BH",
  TRUE ~ "unknown"
)
###################### Plot data ##########################################
### Autosomes ####
Rxy_autosomes <- agg_Rxy_data_C[which(agg_Rxy_data_C$Group.2 == "Autosomes"),] 
Rxy_autosomes <- Rxy_autosomes[which(Rxy_autosomes$comparison == "NARW:SRW" |
                                       Rxy_autosomes$comparison == "NARW:BH" |
                                       Rxy_autosomes$comparison == "SRW:BH"),] 

unique(Rxy_autosomes$comparison)
head(Rxy_autosomes)

Rxy_autosomes$Group.3 <- factor(Rxy_autosomes$Group.3 , levels=c("Modifier", "Low", "Moderate", "High"))
Rxy_autosomes$comparison <- factor(Rxy_autosomes$comparison , levels=c("NARW:SRW", "NARW:BH", "SRW:BH"))

threea <- ggplot((Rxy_autosomes), aes(y=comparison, x=value, color = Group.3)) + 
  geom_point(size = 6.5) +
  theme_classic() +
  scale_color_manual(values=c("gray", "gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(angle = 45, hjust = 0.5)) +
  xlab("Rxy (Autosomes)") +
  ylab("Species pair") +
  theme(legend.position = "none")    



### X ####
Rxy_X <- agg_Rxy_data_C[which(agg_Rxy_data_C$Group.2 == "X chromosome"),] 
Rxy_X <- Rxy_X[which(Rxy_X$comparison == "NARW:SRW" |
                                       Rxy_X$comparison == "NARW:BH" |
                                       Rxy_X$comparison == "SRW:BH"),] 

unique(Rxy_X$comparison)
Rxy_X

Rxy_X$Group.3 <- factor(Rxy_X$Group.3 , levels=c("Modifier", "Low", "Moderate", "High"))
Rxy_X$comparison <- factor(Rxy_autosomes$comparison , levels=c("NARW:SRW", "NARW:BH", "SRW:BH"))

threeb <- ggplot((Rxy_X), aes(y=comparison, x=value, color = Group.3)) + 
  geom_point(size = 6.5) +
  theme_classic() +
  scale_color_manual(values=c("gray", "gold", "orange3", "darkred")) +
  theme(axis.text = element_text(size = 10, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(angle = 45, hjust = .5)) +
  xlab("Rxy (X chromosome)") +
  ylab("Species pair") +
  theme(legend.position = "none")    





rxy_aov <- aov()


#################################################################################################


## 2) is inbreeding precluding increased selection (ie genetic purging through inbreeding)


#################################################################################################
# Test relationship between excess homozygosity and realized genetic load for segregating sites #
#################################################################################################

# autosomes
autosomes_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.ALL.het", header = T)
head(autosomes_F_inter_data)
autosomes_F_inter_data$CHR <- 'Autosomes'

autosomes_F_inter_data$pop <- c('NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW',
                                'SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','SRW','NARW','NARW',
                                'BH','BH','BH','BH','BH','BH','BH','BH','BH','BH')


# X
X_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.ALL.het", header = T)
head(X_F_inter_data)
X_F_inter_data$CHR <- 'X chromosome'

X_F_inter_data$pop <- c('NARW','NARW','NARW','NARW','NARW','NARW','NARW','NARW',
                        'SRW','SRW','SRW','SRW','SRW','NARW','NARW',
                        'BH','BH','BH')


F_inter_data <- rbind(autosomes_F_inter_data, X_F_inter_data)
head(F_inter_data)


head(inter_geneticLoad_indv_data_NS_loads)
#OR
head(inter_geneticLoad_indv_data_NS)
unique(inter_geneticLoad_indv_data_NS$impact)
inter_geneticLoad_indv_data_NS



## align names between data frames ##
names(F_inter_data)[names(F_inter_data) == 'INDV'] <- 'sample'
names(F_inter_data)[names(F_inter_data) == 'pop'] <- 'population'



purge_cor <- merge(inter_geneticLoad_indv_data_NS, F_inter_data, by=c("sample", "CHR"))
head(purge_cor)

purge_cor$derived_F <- purge_cor$alt_hom_F/purge_cor$F
purge_cor$relative_F <- purge_cor$realizedLoad/purge_cor$F


head(purge_cor)


t <- aov(derived_F ~ CHR + population.x + impact, data = purge_cor)
summary(t)

#f Sum Sq Mean Sq F value   Pr(>F)    
#pop           2 0.2521  0.1260  44.742 7.72e-16 ***
#  impact        2 0.8483  0.4241 150.550  < 2e-16 ***
#  CHR           1 0.0271  0.0271   9.634   0.0023 ** 
#  Residuals   144 0.4057  0.0028   



purge_cor$population.x <- factor(purge_cor$population.x, levels=c("NARW", "SRW", "BH"))

summary(purge_cor$F)

threec <- ggplot(purge_cor, aes(x=CHR, y=derived_F, fill = population.x)) + 
  geom_boxplot() + theme_classic() +
  #scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Relative realized mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  +
  ylim(0,1)

  



#############################
#### ROH enrichment #########
#############################


## for differences in deleterious enrichment within ROH ###

enrichment_data <- read.table("/Volumes/cetacea/Genetic_load/_roh/enrichment.txt", header = T)
head(enrichment_data)


enrichment_data$impact <- case_when(
  enrichment_data$Impact == "modifier" ~ "Modifier",
  enrichment_data$Impact == "low" ~ "Low",
  enrichment_data$Impact == "moderate" ~ "Moderate",
  enrichment_data$Impact == "high" ~ "High",
  TRUE ~ "unknown"
)

head(enrichment_data)

enriched_autosomes <- enrichment_data[which(enrichment_data$CHROM == "autosomes"),] 
enriched_autosomes$impact <- factor(enriched_autosomes$impact , levels=c("Modifier", "Low", "Moderate", "High"))
enriched_autosomes$Population <- factor(enriched_autosomes$Population , levels=c("NARW", "SRW", "BH"))

scaleFUN <- function(x) sprintf("%.2f", x)

threed <- ggplot(enriched_autosomes, aes(x=impact, y=freq_enriched, fill = Population)) + 
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Proportion enriched (Autosomes)") +
  scale_y_continuous(labels=scaleFUN) +
  xlab("Mutation impact") +
  theme(legend.position = "none") + 
  ylim(0,1)





t <- aov(freq_enriched ~ Population*Impact, data = enriched_autosomes)
summary(t)

TukeyHSD(t)


enriched_X <- enrichment_data[which(enrichment_data$CHROM == "X"),] 
enriched_X$impact <- factor(enriched_X$impact, levels=c("Modifier", "Low", "Moderate", "High"))
enriched_X$Population <- factor(enriched_X$Population , levels=c("NARW", "SRW", "BH"))

threee <- ggplot(enriched_X, aes(x=impact, y=freq_enriched, fill = Population)) + 
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Proportion enriched (X chromosome)") +
  xlab("Mutation impact") +
  theme(legend.position = "none")  + 
  ylim(0,1)


t <- aov(freq_enriched ~ Population*Impact, data = enriched_X)
summary(t)

TukeyHSD(t)



######### Figure 3 panel ###############3

Figure_3 <- (
  ((threea|threeb) / (threec) + plot_layout(ncol=1,heights=c(1,3)))
             
             | threed|threee) + plot_layout(ncol=3,widths =c(2,1,1)) +
  
  plot_annotation(tag_levels = 'A') 

Figure_3 #view multi-panel figure





###########################################################
#                                                         #             
#                      SUPPLEMENT                         #
#                                                         #
###########################################################

###############################################################################
######### allele/genotype frequencies form per individual data ######################
###############################################################################

# autosomes
autosomes_inter_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.indv.inter.load.stats", header = T)
head(autosomes_inter_geneticLoad_indv_data)
autosomes_inter_geneticLoad_indv_data$CHR <- 'Autosomes'


# X
X_inter_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.indv.inter.load.stats", header = T)
head(X_inter_geneticLoad_indv_data)
X_inter_geneticLoad_indv_data$CHR <- 'X chromosome'

inter_geneticLoad_indv_data <- rbind(autosomes_inter_geneticLoad_indv_data, X_inter_geneticLoad_indv_data)
head(inter_geneticLoad_indv_data)



X_alt_AF <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$CHR == "X chromosome"),]
X_alt_AF$Population <- factor(X_alt_AF$population, levels=c("NARW", "SRW", "BH"))
X_alt_AF$impact <- factor(X_alt_AF$impact, levels=c("Modifier", "Low", "Moderate", "High"))
head(X_alt_AF)

saf_a <- 
  ggplot(X_alt_AF, aes(x=impact, alt_AF , fill = Population)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Derived allele frequency (X chomosome)") +
  xlab("Predicted mutation impact") +
  theme(legend.position = "none")  + 
    ylim(0,0.4)



Autosomes_alt_AF <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$CHR == "Autosomes"),]
Autosomes_alt_AF$Population <- factor(Autosomes_alt_AF$population, levels=c("NARW", "SRW", "BH"))
Autosomes_alt_AF$impact <- factor(Autosomes_alt_AF$impact, levels=c("Modifier", "Low", "Moderate", "High"))
head(Autosomes_alt_AF)

saf_b <- 
  ggplot(Autosomes_alt_AF, aes(x=impact, alt_AF , fill = Population)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Derived allele frequency (Autosomes)") +
  xlab("Predicted mutation impact") +
    theme(legend.position = "none")  + 
    ylim(0,0.4)



Supplemental_allele_Frequencies <- (
  (saf_b|saf_a)) + plot_annotation(tag_levels = 'A') 

Supplemental_allele_Frequencies



############ Plot genotype frequencies ######################

X_alt_AF <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$CHR == "X chromosome"),]
X_alt_AF$Population <- factor(X_alt_AF$population, levels=c("NARW", "SRW", "BH"))
X_alt_AF$impact <- factor(X_alt_AF$impact, levels=c("Modifier", "Low", "Moderate", "High"))
head(X_alt_AF)

gtf_a <- 


ggplot(X_alt_AF, aes(x=impact, alt_hom_F , fill = Population)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Derived homozygote frequency (X chromosome)") +
  xlab("Predicted mutation impact") + theme(legend.position = "none")  + 
  ylim(0,0.3)




Autosomes_alt_AF <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$CHR == "Autosomes"),]
Autosomes_alt_AF$Population <- factor(Autosomes_alt_AF$population, levels=c("NARW", "SRW", "BH"))
Autosomes_alt_AF$impact <- factor(Autosomes_alt_AF$impact, levels=c("Modifier", "Low", "Moderate", "High"))
head(Autosomes_alt_AF)


gtf_b <- 
  
  ggplot(Autosomes_alt_AF, aes(x=impact, alt_hom_F , fill = Population)) + 
  geom_boxplot() +
  theme_classic() +
    scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
    theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Derived homozygote frequency (Autosomes)") +
  xlab("Predicted mutation impact") + theme(legend.position = "none")  + 
  ylim(0, 0.3)


Supplemental_homozyogte_Frequencies <- (
  (gtf_b|gtf_a)) + plot_annotation(tag_levels = 'A') 
  
Supplemental_homozyogte_Frequencies #view multi-panel figure













########################################################################################
#                                        SFS Autosomes                                 #
########################################################################################


sfs <- read.table("/Volumes/cetacea/Genetic_load/SFS/SFS.all.txt", header = T)
sfs <- sfs[which(sfs$Chromosome == "autosomes"),]
sfs

modifier <- sfs[which(sfs$Impact == "modifier"),]
NARW <- sfs[which(sfs$Population == "NARW"),]
modifier$Population <- factor(modifier$Population , levels=c("NARW", "SRW", "BH"))
NARW$Impact <- factor(test$Impact , levels=c("modifier", "low", "moderate", "high"))


  
NARW_sfs <- ggplot(NARW) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
    theme_classic() +
    theme(axis.text = element_text(size = 14, color = "black")) +
    theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
    scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
    xlab("Derived allele frequency") +
    ylab("Denisty") +
    theme(legend.position = "none")  +
  ylim(0,8)


  
SRW <- sfs[which(sfs$Population == "SRW"),]
SRW$Impact <- factor(SRW$Impact , levels=c("modifier", "low", "moderate", "high"))

SRW_sfs <- ggplot(SRW) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
  xlab("Derived allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  +
  ylim(0,8)


BH<- sfs[which(sfs$Population == "BH"),]
BH$Impact <- factor(BH$Impact , levels=c("modifier", "low", "moderate", "high"))

BH_sfs <- ggplot(BH) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
  xlab("Derived allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none") +
  ylim(0,8)

  
supplement_sfs <- (NARW_sfs|SRW_sfs|BH_sfs) +  
  plot_annotation(tag_levels = 'A') 

supplement_sfs









low <- sfs[which(sfs$Impact == "low"),]
low$Population <- factor(low$Population , levels=c("NARW", "SRW", "BH"))
sfs_b1 <-  ggplot(low) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  

moderate <- sfs[which(sfs$Impact == "moderate"),]
moderate$Population <- factor(moderate$Population , levels=c("NARW", "SRW", "BH"))
sfs_c1 <- ggplot(moderate) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  


high <- sfs[which(sfs$Impact == "high"),]
high$Population <- factor(high$Population , levels=c("NARW", "SRW", "BH"))
sfs_d1 <- ggplot(high) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  



supplement_sfs <- (
  (sfs_a1|sfs_b1)/
    (sfs_c1|sfs_d1)) +  
  plot_annotation(tag_levels = 'A') 

supplement_sfs


inset <- high[which(high$Frequency >= "0.5"),]
inset$Population <- factor(inset$Population , levels=c("NARW", "SRW", "BH"))
ggplot(inset) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  ylab("Allele frequency") +
  xlab("Denisty") +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())
  



########################################################################################
#                                    SFS X chromosome                                 #
########################################################################################


sfs <- read.table("/Volumes/cetacea/Genetic_load/SFS/SFS.all.txt", header = T)
sfs <- sfs[which(sfs$Chromosome == "X"),]
sfs

NARWX <- sfs[which(sfs$Population == "NARW"),]
NARWX$Impact <- factor(NARWX$Impact , levels=c("modifier", "low", "moderate", "high"))



NARWX_sfs <- ggplot(NARWX) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
  xlab("Derived allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  +
  ylim(0,8)



SRWX <- sfs[which(sfs$Population == "SRW"),]
SRWX$Impact <- factor(SRWX$Impact , levels=c("modifier", "low", "moderate", "high"))

SRWX_sfs <- ggplot(SRWX) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
  xlab("Derived allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  +
  ylim(0,8)


BHX<- sfs[which(sfs$Population == "BH"),]
BHX$Impact <- factor(BHX$Impact , levels=c("modifier", "low", "moderate", "high"))

BHX_sfs <- ggplot(BHX) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Impact), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("gray", "gold", "orange3", "darkred")) +
  xlab("Derived allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none") +
  ylim(0,8)


supplement_sfsX <- (NARWX_sfs|SRWX_sfs|BHX_sfs) +  
  plot_annotation(tag_levels = 'A') 

supplement_sfsX






######## by impact ###############

modifier <- sfs[which(sfs$Impact == "modifier"),]
modifier$Population <- factor(modifier$Population , levels=c("NARW", "SRW", "BH"))

sfs_a1 <- ggplot(modifier) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  




low <- sfs[which(sfs$Impact == "low"),]
low$Population <- factor(low$Population , levels=c("NARW", "SRW", "BH"))
sfs_b1 <-  ggplot(low) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  

moderate <- sfs[which(sfs$Impact == "moderate"),]
moderate$Population <- factor(moderate$Population , levels=c("NARW", "SRW", "BH"))
#sfs_c1 <- 

  ggplot(moderate) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  


high <- sfs[which(sfs$Impact == "high"),]
high$Population <- factor(high$Population , levels=c("NARW", "SRW", "BH"))
#sfs_d1 <- 
  
  ggplot(high) + geom_histogram(aes(x=Frequency, y=..density..,weight=Count, fill = Population), position="dodge", binwidth = 0.1) +
  theme(legend.position = "none")  + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  scale_fill_manual(values=c("violetred4", "#D6C350", "#8EB6B5")) +
  xlab("Allele frequency") +
  ylab("Denisty") +
  theme(legend.position = "none")  



supplement_sfs <- (
  (sfs_a1|sfs_b1|sfs_c1|sfs_d1)) +  
  plot_annotation(tag_levels = 'A') 

supplement_sfs




########################################################################################
#                                    Site statistics                                   #
########################################################################################


segregate <- read.table("/Volumes/cetacea/Genetic_load/SFS/site_stats.all.txt", header = T)
segregate
site_stats <- melt(data = segregate, id.vars = c("Population", "Chromosome", "Impact"), 
                   variable.name = "Estimate", 
                   variable.value = "value")

site_stats_autosomes <- site_stats[which(site_stats$Chromosome == "autosomes"),]

site_stats_autosomes <- site_stats_autosomes[which(site_stats_autosomes$Impact != "modifier"),]
site_stats_autosomes <- site_stats_autosomes[which(site_stats_autosomes$Estimate == "proportion_fixed_reference" | 
                                           site_stats_autosomes$Estimate == "proportion_fixed_derived" |
                                           site_stats_autosomes$Estimate == "proportion_segregating_derived" |
                                           site_stats_autosomes$Estimate == "proportion_sites_with_derived_allele"),]

site_stats_autosomes$Population <- factor(site_stats_autosomes$Population , levels=c("NARW", "SRW", "BH"))
site_stats_autosomes$Estimate <- factor(site_stats_autosomes$Estimate , levels=c("proportion_fixed_reference", 
                                                                                 "proportion_fixed_derived", 
                                                                                 "proportion_segregating_derived",
                                                                                 "proportion_sites_with_derived_allele"))

auto_counts <- 
ggplot((site_stats_autosomes), aes(x=Estimate, y=value, fill = Population)) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  ylab("Proportion of sites (Autosomes)") +
  xlab("Class") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("proportion_fixed_reference" = "Fixed for ancestral allele", 
                            "proportion_fixed_derived" = "Fixed for derived allele",
                            "proportion_segregating_derived" = "Derived allele segregating",
                            "proportion_sites_with_derived_allele" = "Derived allele present"))





site_stats_X <- site_stats[which(site_stats$Chromosome == "X"),]

site_stats_X <- site_stats_X[which(site_stats_X$Impact != "modifier"),]
site_stats_X <- site_stats_X[which(site_stats_X$Estimate == "proportion_fixed_reference" | 
                                                     site_stats_X$Estimate == "proportion_fixed_derived" |
                                                     site_stats_X$Estimate == "proportion_segregating_derived" |
                                                     site_stats_X$Estimate == "proportion_sites_with_derived_allele"),]

site_stats_X$Population <- factor(site_stats_X$Population , levels=c("NARW", "SRW", "BH"))
site_stats_X$Estimate <- factor(site_stats_X$Estimate , levels=c("proportion_fixed_reference", 
                                                                                 "proportion_fixed_derived", 
                                                                                 "proportion_segregating_derived",
                                                                                 "proportion_sites_with_derived_allele"))

X_counts <- 
ggplot((site_stats_X), aes(x=Estimate, y=value, fill = Population)) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("darkslategray", "aquamarine4", "#8EB6B5")) +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  ylab("Proportion of sites (Autosomes)") +
  xlab("Class") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("proportion_fixed_reference" = "Fixed for ancestral allele", 
                            "proportion_fixed_derived" = "Fixed for derived allele",
                            "proportion_segregating_derived" = "Derived allele segregating",
                            "proportion_sites_with_derived_allele" = "Derived allele present"))


supplement_sites <- (
  (auto_counts/X_counts)) +  
  plot_annotation(tag_levels = 'A') 

supplement_sites
