########################################## NARW ##################################3

## read in (individual) genetic load data ##
## for all sites (nonsynonymous and synonymous) ##

# autosomes
autosomes_NARW_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.indv.NARW.load.stats", header = T)
head(autosomes_NARW_geneticLoad_indv_data)
autosomes_NARW_geneticLoad_indv_data$CHR <- 'Autosomes'


# X
X_NARW_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/X.indv.NARW.load.stats", header = T)
head(X_NARW_geneticLoad_indv_data)
X_NARW_geneticLoad_indv_data$CHR <- 'X'

geneticLoad_indv_NARW_data <- rbind(autosomes_NARW_geneticLoad_indv_data, X_NARW_geneticLoad_indv_data)
geneticLoad_indv_NARW_data

### set selection coefficient ###
geneticLoad_indv_NARW_data$s <- case_when(
  geneticLoad_indv_NARW_data$impact == "Modifier" ~ "0",
  geneticLoad_indv_NARW_data$impact == "Low" ~ ".1",
  geneticLoad_indv_NARW_data$impact == "Moderate" ~ ".3",
  geneticLoad_indv_NARW_data$impact == "High" ~ "6",
  TRUE ~ "unknown"
)



head(geneticLoad_indv_NARW_data)
### pull dominance coefficient ###
#geneticLoad_indv_NARW_data$h <- by(geneticLoad_indv_NARW_data, seq_len(nrow(geneticLoad_indv_NARW_data)), function(row) rexp(1, 10))
#geneticLoad_indv_NARW_data$h <- 0.02

geneticLoad_indv_NARW_data <- geneticLoad_indv_NARW_data %>% mutate_at(c('s'), as.numeric)

head(geneticLoad_indv_NARW_data)


### calculate total mutation load ###
geneticLoad_indv_NARW_data$TotalLoad <- geneticLoad_indv_NARW_data$alt_hom_F*geneticLoad_indv_NARW_data$s + 
  0.5*geneticLoad_indv_NARW_data$s*geneticLoad_indv_NARW_data$het_F 

### calculate realized mutation load ###
geneticLoad_indv_NARW_data$realizedLoad <- geneticLoad_indv_NARW_data$alt_hom_F*geneticLoad_indv_NARW_data$s #+ geneticLoad_indv_NARW_data$het_F*geneticLoad_indv_NARW_data$h*geneticLoad_indv_NARW_data$s


### calculate masked mutation load ###
geneticLoad_indv_NARW_data$MaskedLoad <- geneticLoad_indv_NARW_data$het_F*geneticLoad_indv_NARW_data$s*(.5)

### remove modifier sites ###
geneticLoad_indv_NARW_data_NS <- geneticLoad_indv_NARW_data[which(geneticLoad_indv_NARW_data$impact != "Modifier"),]
head(geneticLoad_indv_NARW_data_NS)




### calculate standard error and mean for each group ###
se <- function(x) sd(x) / sqrt(length(x))

NARW_total_load_summary <- geneticLoad_indv_NARW_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(TotalLoad),
            SE = se(TotalLoad),
            SD = sd(TotalLoad))

NARW_total.aov2 <- aov(TotalLoad  ~ population + CHR, data = geneticLoad_indv_NARW_data_NS)
summary(NARW_total.aov2)


#Df Sum Sq  Mean Sq F value Pr(>F)
#population   1 0.0006 0.000551   0.036  0.850
#CHR          1 0.0079 0.007892   0.520  0.474
#Residuals   57 0.8657 0.015188

NARW_realized_load_summary <- geneticLoad_indv_NARW_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(realizedLoad),
            SE = se(realizedLoad),
            SD = sd(realizedLoad))



NARW_realized.aov2 <- aov(realizedLoad  ~ population + CHR, data = geneticLoad_indv_NARW_data_NS)
summary(NARW_realized.aov2)

#Df  Sum Sq  Mean Sq F value  Pr(>F)   
#population   1 0.00300 0.003001   0.865 0.35619   
#CHR          1 0.02789 0.027892   8.041 0.00632 **
#  Residuals   57 0.19771 0.003469 


NARW_masked_load_summary <- geneticLoad_indv_NARW_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(MaskedLoad),
            SE = se(MaskedLoad),
            SD = sd(MaskedLoad))

NARW_masked.aov2 <- aov(MaskedLoad  ~ population + CHR, data = geneticLoad_indv_NARW_data_NS)
summary(NARW_masked.aov2)

#Df Sum Sq Mean Sq F value Pr(>F)   
#population   1 0.0061 0.00612   0.806 0.3732   
#CHR          1 0.0655 0.06546   8.613 0.0048 **
#  Residuals   57 0.4332 0.00760 

######## Plot data #############

## Total additive genetic load ##
#d1 <- 
  ggplot((NARW_total_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#65A3B6","#B2BC80")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Total mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")   

## Realized genetic load ##
#d2 <- 
  ggplot((NARW_realized_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#65A3B6","#B2BC80")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Realized mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  

## Masked genetic load ##
#d3 <- 
  ggplot((NARW_masked_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#65A3B6","#B2BC80")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Masked mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")    



ggarrange(d1, d2, d3, ncol = 3, nrow = 1)  







## read in (individual) genetic load data ##
## for all sites (nonsynonymous and synonymous) ##

# autosomes
autosomes_inter_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.indv.inter.load.stats", header = T)
head(autosomes_inter_geneticLoad_indv_data)
autosomes_inter_geneticLoad_indv_data$CHR <- 'Autosomes'


# X
X_inter_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.indv.inter.load.stats", header = T)
head(X_inter_geneticLoad_indv_data)
X_inter_geneticLoad_indv_data$CHR <- 'X'

inter_geneticLoad_indv_data <- rbind(autosomes_inter_geneticLoad_indv_data, X_inter_geneticLoad_indv_data)


### set selection coefficient ###
inter_geneticLoad_indv_data$s <- case_when(
  inter_geneticLoad_indv_data$impact == "Modifier" ~ "0",
  inter_geneticLoad_indv_data$impact == "Low" ~ ".2",
  inter_geneticLoad_indv_data$impact == "Moderate" ~ ".5",
  inter_geneticLoad_indv_data$impact == "High" ~ "1",
  TRUE ~ "unknown"
)

### pull dominance coefficient ###

#inter_geneticLoad_indv_data$h <- by(inter_geneticLoad_indv_data, seq_len(nrow(inter_geneticLoad_indv_data)), function(row) rexp(1, 10))
inter_geneticLoad_indv_data$h <- 0.05

inter_geneticLoad_indv_data <- inter_geneticLoad_indv_data %>% mutate_at(c('h', 's'), as.numeric)
inter_geneticLoad_indv_data_NS <- inter_geneticLoad_indv_data[which(inter_geneticLoad_indv_data$impact != "Modifier"),]

head(inter_geneticLoad_indv_data_NS)

### calculate total genetic load ###
inter_geneticLoad_indv_data_NS$TotalLoad <- inter_geneticLoad_indv_data_NS$alt_hom_F*inter_geneticLoad_indv_data_NS$s + 
  0.5*inter_geneticLoad_indv_data_NS$s*inter_geneticLoad_indv_data_NS$het_F 

total.aov2 <- aov(TotalLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(total.aov2)



### calculate realized genetic load ###
inter_geneticLoad_indv_data_NS$realizedLoad <- inter_geneticLoad_indv_data_NS$alt_hom_F*inter_geneticLoad_indv_data_NS$s +
  inter_geneticLoad_indv_data_NS$het_F*inter_geneticLoad_indv_data_NS$h*inter_geneticLoad_indv_data_NS$s

realized.aov2 <- aov(realizedLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(realized.aov2)
realized.aov2


### calculate masked genetic load ###
inter_geneticLoad_indv_data_NS$MaskedLoad <- inter_geneticLoad_indv_data_NS$het_F*(0.5-inter_geneticLoad_indv_data_NS$h)*inter_geneticLoad_indv_data_NS$s 

masked.aov2 <- aov(MaskedLoad  ~ population + CHR, data = inter_geneticLoad_indv_data_NS)
summary(masked.aov2)

head(inter_geneticLoad_indv_data_NS)




### calculate standard error and mean for each group ###
se <- function(x) sd(x) / sqrt(length(x))

total_load_summary <- inter_geneticLoad_indv_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(TotalLoad),
            SE = se(TotalLoad),
            SD = sd(TotalLoad))

realized_load_summary <- inter_geneticLoad_indv_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(realizedLoad),
            SE = se(realizedLoad),
            SD = sd(realizedLoad))

masked_load_summary <- inter_geneticLoad_indv_data_NS %>% 
  group_by(population, CHR) %>% 
  summarise(mean = mean(MaskedLoad),
            SE = se(MaskedLoad),
            SD = sd(MaskedLoad))




######## Plot data #############

## Total additive genetic load ##
total_load_summary$population <- factor(total_load_summary$population, levels=c("BH", "SRW", "NARW"))
a1 <- ggplot((total_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Total additive mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")   

## Realized genetic load ##
realized_load_summary$population <- factor(realized_load_summary$population, levels=c("BH", "SRW", "NARW"))
a2 <- ggplot((realized_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Realized mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  

## Masked genetic load ##
masked_load_summary$population <- factor(masked_load_summary$population, levels=c("BH", "SRW", "NARW"))
a3 <- ggplot((masked_load_summary), aes(x=CHR, y=mean, fill = population)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=CHR, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Masked mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")    



ggarrange(a1, a2, a3, ncol = 3, nrow = 1)







#################################################################################################
# Test relationship between excess homozygosity and realized genetic load for segregating sites #
#################################################################################################


head(geneticLoad_indv_NARW_data_NS)
unique(geneticLoad_indv_NARW_data_NS$impact)


## get F data ##
head(NARW_F_data)

## align names between data frames ##
names(NARW_F_data)[names(NARW_F_data) == 'INDV'] <- 'sample'



purge_cor <- merge(geneticLoad_indv_NARW_data_NS, NARW_F_data, by=c("sample", "CHR"))
head(purge_cor)

purge_cor$derived_F <- purge_cor$alt_hom_F/purge_cor$F
purge_cor$relative_F <- purge_cor$realizedLoad/purge_cor$F


head(purge_cor)


t <- aov(relative_F ~ population + impact + CHR, data = purge_cor)
summary(t)

#df Sum Sq Mean Sq F value   Pr(>F)    
#pop           2 0.2521  0.1260  44.742 7.72e-16 ***
#  impact        2 0.8483  0.4241 150.550  < 2e-16 ***
#  CHR           1 0.0271  0.0271   9.634   0.0023 ** 
#  Residuals   144 0.4057  0.0028   



purge_cor$pop <- factor(purge_cor$pop, levels=c("BH", "SRW", "NARW"))

summary(purge_cor$F)

ggplot(purge_cor, aes(x=CHR, y=derived_F, fill = population)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("#8EB6B5", "#D6C350")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Relative realized mutation load") +
  xlab("Genomic region")




###################################################
#################      DAPC       #################
###################################################


## read in autosomes gVCF ##
vcf <- read.vcfR("/Volumes/cetacea/Genetic_load/NARW_females/forDAPC/all.ns.vcf.gz")


## convert to genlight object ##
vcf_gi <- vcfR2genlight(vcf)
print(vcf_gi)
vcf_gi$ind.names
vcf_gi$loc.names[1:10]


## assign fecundity groups ##
individuals <- vcf_gi$ind.names

pop_map <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/forDAPC/popmap.txt", header = T)
head(pop_map)

vcf_gi$pop <- pop_map$group
################################



## plot g plot ##
#glPlot(vcf_gi, posi="bottomleft")

# compute the PCA on the SNP genotypes and plot it:
pca1 <- glPca(vcf_gi, nf=4) # nf = number of PC axes to retain (here, 4)
pca1 # prints summary


# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=vcf_gi$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2")
legend("topleft", 
       legend=unique(vcf_gi$pop), 
       pch=20, 
       col=c("black", "red"))


# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999))
# Get their locus names
vcf_gi$loc.names[which(quantile(abs(pca1$loadings))>0.999)]

threshold<-quantile(abs(pca1$loadings),0.999)

vcf_gi$loc.names[which(abs(pca1$loadings)>threshold)]

vcf_gi$loc.names[which(quantile(abs(pca1$loadings),0.999)>0.0770)]


# Run the DAPC using disease status to group samples
disease.dapc <- dapc(vcf_gi, pop=vcf_gi$pop, n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T)

# Scatterplot of results
scatter.dapc(disease.dapc, grp=vcf_gi$pop, legend=T)

# Plot the posterior assignment probabilities to each group
compoplot(disease.dapc)

# Which loci contribute the most to distinguishing Healthy vs. Sick individuals?
loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))



############################################################################
################################## Scratch #################################
############################################################################

###################################################
######### Rxy per individual ######################
###################################################

head(geneticLoad_indv_NARW_data)

######## Prepare data #############



high_geneticLoad_indv_NARW_data <- geneticLoad_indv_NARW_data[which(geneticLoad_indv_NARW_data$population == "HIGH"),]
low_geneticLoad_indv_NARW_data <- geneticLoad_indv_NARW_data[which(geneticLoad_indv_NARW_data$population == "LOW"),]


agg_geneticLoad_indv_NARW_data <- aggregate(geneticLoad_indv_NARW_data$alt_hom_F,
                                            by=list(geneticLoad_indv_NARW_data$population,
                                                    geneticLoad_indv_NARW_data$CHR,
                                                    geneticLoad_indv_NARW_data$impact), FUN=mean)
agg_geneticLoad_indv_NARW_data


geneticLoad_indv_NARW_data_loads <- agg_geneticLoad_indv_NARW_data %>% spread(Group.1, x)
geneticLoad_indv_NARW_data_loads

geneticLoad_indv_NARW_data_loads$Rxy <- (geneticLoad_indv_NARW_data_loads$LOW/geneticLoad_indv_NARW_data_loads$HIGH)
geneticLoad_indv_NARW_data_loads



ggarrange(a2, a3, a1, ncol = 3, nrow = 1)


###################################################
######### Rxy per site ######################
###################################################

## autosomes 
autosome_NARW_site_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.NARW.delta_load", head = T)
head(autosome_NARW_site_data)
autosome_NARW_site_data$CHR <- 'Autosomes'

## X
X_NARW_site_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/X.NARW.delta_load", head = T)
X_NARW_site_data$CHR <- 'X chromosome'

NARW_site_data <- rbind(autosome_NARW_site_data, X_NARW_site_data)
head(NARW_site_data)


unique(NARW_site_data$estimate)

Rxy_NARW_site_data_low <- NARW_site_data[which(NARW_site_data$estimate == "lowFecund_Alt_AF"),]
Rxy_NARW_site_data_high <- NARW_site_data[which(NARW_site_data$estimate == "highFecund_Alt_AF"),]

head(Rxy_NARW_site_data_high)

Rxy_NARW_site_data <- cbind(Rxy_NARW_site_data_low, Rxy_NARW_site_data_high)
Rxy_NARW_site_data_b = subset(Rxy_NARW_site_data, select = -c(1,2,3) )
head(Rxy_NARW_site_data_b)

na.omit(Rxy_NARW_site_data_b)
Rxy_NARW_site_data_b$Rxy <- (Rxy_NARW_site_data_b$value/Rxy_NARW_site_data_b$value.1)

Rxy_NARW_site_data_b$Impact <- factor(Rxy_NARW_site_data_b$Impact , levels=c("Modifier", "Low", "Moderate", "High"))
ggplot((Rxy_NARW_site_data_b), aes(x=CHR, y=Rxy, fill = Impact)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("gray", "yellow", "orange", "red")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Masked mutation load") +
  xlab("Genomic region")



########################################################################################
#                       F - genome-wide  (TEST)                                        #
########################################################################################


## read in (individual) homozygote data

# autosomes
autosomes_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.ALL.het", header = T)
head(autosomes_F_inter_data)
autosomes_F_inter_data$CHR <- 'Autosomes'

# X
X_F_inter_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.ALL.het", header = T)
head(X_F_inter_data)
X_F_inter_data$CHR <- 'X'

F_inter_data <- rbind(autosomes_F_inter_data, X_F_inter_data)

F_inter_data_NARW_test <- F_inter_data[ which(F_inter_data$INDV =='EGL00252-1' |
                                                F_inter_data$INDV == 'EGL013-3qa'|
                                                F_inter_data$INDV == 'EGL183-1'  |
                                                F_inter_data$INDV == 'EGL254-1'	 |
                                                F_inter_data$INDV == 'EGL276-1'	 |
                                                F_inter_data$INDV == 'EGL308-1a' |
                                                F_inter_data$INDV == 'EGL312-1a' |
                                                F_inter_data$INDV == 'EGL336_1b' |
                                                F_inter_data$INDV == 'SID179132' |
                                                F_inter_data$INDV == 'SID181803'), ] 

F_inter_data_NARW_test


F_inter_data_NARW_test$fecundity_group <- c('Low','Low','High','Low','Low','Low',
                                            'High','High', 'High','High')



## box plot ##
ggplot(F_inter_data_NARW_test, aes(x=CHR, y=F, fill = fecundity_group)) + 
  geom_boxplot() + theme_classic() +
  scale_fill_manual(values=c("deepskyblue4", "darkslateblue")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Wright's F") +
  xlab("Genomic region")


agg_F_inter_data_NARW_test <- aggregate(F_inter_data_NARW_test$F,  
                                        by=list(F_inter_data_NARW_test$fecundity_group,
                                                F_inter_data_NARW_test$CHR), FUN=mean)

agg_F_inter_data_NARW_test













#################################################################################

########################## Supplemental #######################################

#################################################################################

### Het:Hom scratch analysis.

###### new test with individual data ############

head(inter_geneticLoad_indv_data)
unique(inter_geneticLoad_indv_data_NS$impact)


## get F data ##
head(F_inter_data)



agg_Het <- aggregate(inter_geneticLoad_indv_data$het_F,  by=list(inter_geneticLoad_indv_data$population,
                                                                 inter_geneticLoad_indv_data$impact,
                                                                 inter_geneticLoad_indv_data$CHR), FUN=mean)

agg_Het


agg_Hom <- aggregate(inter_geneticLoad_indv_data$alt_hom_F,  by=list(inter_geneticLoad_indv_data$population,
                                                                     inter_geneticLoad_indv_data$impact,
                                                                     inter_geneticLoad_indv_data$CHR), FUN=mean)

agg_Hom



########################################################################################
#                        Genetic load by site - interspecific                           #
########################################################################################


## read in site genetic load data ##

# autosomes
autosomes_inter_site_load_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/autosomes.inter.site.load", header = T)
head(autosomes_inter_site_load_data)
autosomes_inter_site_load_data$Chromosome <- 'Autosomes'


## box plot ##
ggplot(autosomes_inter_site_load_data, aes(x=Impact, y=alternate_allele_frequency, fill = Population)) + 
  geom_boxplot() + theme_classic() +
  #scale_fill_manual(values=c("deepskyblue4", "darkslateblue")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Derived allele freq") +
  xlab("mpact")

# X
X_inter_site_load_data <- read.table("/Volumes/cetacea/Genetic_load/interspecific/X.inter.site.load", header = T)
head(X_inter_site_load_data)
X_inter_site_load_data$Chromosome <- 'X'



## box plot ##
ggplot(X_inter_site_load_data, aes(x=Impact, y=alternate_allele_frequency, fill = Population)) + 
  geom_boxplot() + theme_classic() +
  #scale_fill_manual(values=c("deepskyblue4", "darkslateblue")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Wright's F") +
  xlab("Genomic region")




## Derived allele frequency for each population ##
inter_site_load_data <- rbind(autosomes_inter_site_load_data, X_inter_site_load_data)
head(inter_site_load_data)


agg_inter_site_load_data <- aggregate(inter_site_load_data$alternate_allele_frequency,  
                                      by=list(inter_site_load_data$Population,
                                              inter_site_load_data$Impact,
                                              inter_site_load_data$Chromosome), FUN=mean)


agg_inter_site_load_data







## check input data ##
head(inter_site_load_data)


## plot p values for chi-square tests ##

test_B <- na.omit(inter_site_load_data)


X_HWE_ns <- test_B[which(test_B$Chromosome == "X" & test_B$Impact != "Modifier"),]
head(X_HWE_ns)


## correlation matrix ##

BH_cor_test <- test_B[which(test_B$Population == "BH" & test_B$Chromosome == "Autosomes"),]
my_data <- BH_cor_test[,c(11,12,13,14,15,16,17,18,19,20,21,22)]
head(my_data)

res <- cor(my_data)
round(res, 2)

library(corrplot)

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)




ggplot(X_HWE_ns, aes(x=-log10(P_HET_DEFICIT), y=alternate_genotype_frequency, color = Population)) + 
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm) +
  scale_color_manual(values=c("deepskyblue4", "darkslateblue", "cadetblue4")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Frequency of derived allele") +
  xlab("Genomic region")



Autosomes_HWE_ns <- test_B[which(test_B$Chromosome == "Autosomes" & test_B$Impact != "Modifier"),]
head(Autosomes_HWE_ns)

ggplot(Autosomes_HWE_ns, aes(ChiSq_HWE, y=alternate_genotype_frequency, color = Population)) + 
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm) +
  scale_color_manual(values=c("deepskyblue4", "darkslateblue", "cadetblue4")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Frequency of derived allele") +
  xlab("Genomic region")








##### Rxy
#######################################################################################
## divide data by species ##
Rxy_BH <- inter_site_load_data[which(inter_site_load_data$Population == "BH"),] %>% 
  subset(select = c(1,2,4,7,17,23)) %>%
  data.table::setnames('alternate_allele_frequency','BH_alternate_allele_frequency') %>%
  data.table::setnames('HOM2.','BH_HOM')


Rxy_NARW <- inter_site_load_data[which(inter_site_load_data$Population == "NARW"),] %>% 
  subset(select = c(7,17)) %>%
  data.table::setnames('alternate_allele_frequency','NARW_alternate_allele_frequency') %>%
  data.table::setnames('HOM2.','NARW_HOM')


Rxy_SRW <- inter_site_load_data[which(inter_site_load_data$Population == "SRW"),] %>% 
  subset(select = c(7,17)) %>%
  data.table::setnames('alternate_allele_frequency','SRW_alternate_allele_frequency') %>%
  data.table::setnames('HOM2.','SRW_HOM')


Rxy_BH_NARW_SRW <- cbind(Rxy_BH, Rxy_NARW, Rxy_SRW)
head(Rxy_BH_NARW_SRW)

Rxy_BH_NARW_SRW$NARW_SRW <- Rxy_BH_NARW_SRW$NARW_alternate_allele_frequency / Rxy_BH_NARW_SRW$SRW_alternate_allele_frequency

Rxy_BH_NARW_SRW$Impact <- factor(Rxy_BH_NARW_SRW$Impact , levels=c("Modifier", "Low", "Moderate", "High"))
Rxy_BH_NARW_SRW <- Rxy_BH_NARW_SRW[which(Rxy_BH_NARW_SRW$NARW_SRW != 0),] 
ggplot((Rxy_BH_NARW_SRW), aes(x=Chromosome, y=NARW_SRW, fill = Impact)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values=c("gray", "yellow", "orange", "red")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Masked mutation load") +
  xlab("Genomic region") +
  ylim(-2,2)












### get correlations between F and mutation load ###
Autosomes_BH_purge_cor <- purge_cor[which(purge_cor$CHR  == "Autosomes" & purge_cor$pop == "BH" & purge_cor$F > 0.35),]
cor.test(Autosomes_BH_purge_cor$alt_AF, Autosomes_BH_purge_cor$F)

Autosomes_SRW_purge_cor <- purge_cor[which(purge_cor$CHR  == "Autosomes" & purge_cor$pop == "SRW" & purge_cor$F > 0.3),]
cor.test(Autosomes_SRW_purge_cor$alt_AF, Autosomes_SRW_purge_cor$F)

Autosomes_NARW_purge_cor <- purge_cor[which(purge_cor$CHR  == "Autosomes" & purge_cor$pop == "NARW"),]
cor.test(Autosomes_NARW_purge_cor$alt_AF, Autosomes_NARW_purge_cor$F)




## plot values from Autosomes ##

## BH ##
ggplot(Autosomes_BH_purge_cor, aes(x=F, y=alt_AF, color = impact)) + 
  geom_point(size = 1.5) +
  theme_classic() + 
  geom_smooth(method = lm) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Realized genetic load") +
  xlab("Wright's F")

## NARW ##
ggplot(Autosomes_NARW_purge_cor, aes(x=F, y=alt_hom_F, color = impact)) + 
  geom_point(size = 1.5) +
  theme_classic() + 
  geom_smooth(method = lm) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Realized genetic load") +
  xlab("Wright's F")

## SRW ##
ggplot(Autosomes_SRW_purge_cor, aes(x=F, y=alt_hom_F, color = impact)) + 
  geom_point(size = 1.5) +
  theme_classic() + 
  geom_smooth(method = lm) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Realized genetic load") +
  xlab("Wright's F")


ggarrange(d1, d3, d2, ncol = 3, nrow = 1)






#####################################
#### get ratio of ns HWE to s HWE ###
####################################


### 1) Is selection acting on putatively purged sites ? ###

head(inter_site_load_data)

test <- na.omit(inter_site_load_data)

agg_Chi_HWE <- aggregate(test$P_HWE,  by=list(test$Population,
                                              test$Impact,
                                              test$Chromosome), FUN=mean)



agg_Chi_HWE$p <- -log10(agg_Chi_HWE$x)
agg_Chi_HWE<-  subset(agg_Chi_HWE, select=-(x))
head(agg_Chi_HWE)

Chi_HWE_ratio <- agg_Chi_HWE %>% spread(Group.2, p)
head(Chi_HWE_ratio)

Chi_HWE_ratio$ratio_high <- Chi_HWE_ratio$High/Chi_HWE_ratio$Modifier
Chi_HWE_ratio$ratio_moderate <- Chi_HWE_ratio$Moderate/Chi_HWE_ratio$Modifier
Chi_HWE_ratio$ratio_low <- Chi_HWE_ratio$Low/Chi_HWE_ratio$Modifier
head(Chi_HWE_ratio)

Chi_HWE_ratio_for_plot <- Chi_HWE_ratio %>% gather("estimate", "value", -Group.1, -Group.3)
head(Chi_HWE_ratio_for_plot)

Chi_HWE_ratio_for_plot <- Chi_HWE_ratio_for_plot[which(Chi_HWE_ratio_for_plot$estimate  == "ratio_low" |
                                                         Chi_HWE_ratio_for_plot$estimate  == "ratio_moderate" | 
                                                         Chi_HWE_ratio_for_plot$estimate  == "ratio_high"),]

head(Chi_HWE_ratio_for_plot)

Autosomes_agg_Chi_HWE <- Chi_HWE_ratio_for_plot[which(Chi_HWE_ratio_for_plot$Group.3  == "Autosomes"),]
Autosomes_agg_Chi_HWE$estimate <- factor(Autosomes_agg_Chi_HWE$estimate , levels=c("ratio_low", "ratio_moderate", "ratio_high"))
head(Autosomes_agg_Chi_HWE)
###### Plot data ##############

Autosomes_agg_Chi_HWE$Group.1 <- factor(Autosomes_agg_Chi_HWE$Group.1, levels=c("BH", "SRW", "NARW"))
c1 <- ggplot(Autosomes_agg_Chi_HWE, aes(x=estimate, y=value, color = Group.1)) + 
  geom_point(size = 8) +
  theme_classic() + 
  scale_color_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Chi-square N/Chi-square S") +
  xlab("Mutation impact") #+
theme(legend.position = "none")   



X_agg_Chi_HWE <- Chi_HWE_ratio_for_plot[which(Chi_HWE_ratio_for_plot$Group.3  == "X"),]
X_agg_Chi_HWE$estimate <- factor(X_agg_Chi_HWE$estimate , levels=c("ratio_low", "ratio_moderate", "ratio_high"))

X_agg_Chi_HWE$Group.1 <- factor(X_agg_Chi_HWE$Group.1, levels=c("BH", "SRW", "NARW"))
c2 <- ggplot(X_agg_Chi_HWE, aes(x=estimate, y=value, color = Group.1)) + 
  geom_point(size = 8) +
  theme_classic() + 
  scale_color_manual(values=c("#8EB6B5", "#D6C350", "#5599B0")) +
  theme(axis.text = element_text(size = 14, color = "black")) +
  theme(axis.title = element_text(size = 18, color = "black", face = "bold")) +
  ylab("Chi-square N/Chi-square S") +
  xlab("Mutation impact") +
  theme(legend.position = "none")  


ggarrange(c1, c2, ncol = 2, nrow = 1)







###################################################
################ NARW Specific ####################
###################################################
######## Genetic load individual-based ############
###################################################


# autosomes
autosomes_NARW_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/autosomes.indv.NARW.load.stats", header = T)
tail(autosomes_NARW_geneticLoad_indv_data)
autosomes_NARW_geneticLoad_indv_data$CHR <- 'Autosomes'


# X
X_NARW_geneticLoad_indv_data <- read.table("/Volumes/cetacea/Genetic_load/NARW_females/X.indv.NARW.load.stats", header = T)
head(X_NARW_geneticLoad_indv_data)
X_NARW_geneticLoad_indv_data$CHR <- 'X chromosome'

geneticLoad_indv_NARW_data <- rbind(autosomes_NARW_geneticLoad_indv_data, X_NARW_geneticLoad_indv_data)
geneticLoad_indv_NARW_data

### set selection coefficient ###
geneticLoad_indv_NARW_data$s <- case_when(
  geneticLoad_indv_NARW_data$impact == "Modifier" ~ "0",
  geneticLoad_indv_NARW_data$impact == "Low" ~ ".1",
  geneticLoad_indv_NARW_data$impact == "Moderate" ~ ".3",
  geneticLoad_indv_NARW_data$impact == "High" ~ ".6",
  TRUE ~ "unknown"
)



head(geneticLoad_indv_NARW_data)

geneticLoad_indv_NARW_data <- geneticLoad_indv_NARW_data %>% mutate_at(c('s'), as.numeric)

head(geneticLoad_indv_NARW_data)


### calculate total genetic load ###
geneticLoad_indv_NARW_data$TotalLoad <- geneticLoad_indv_NARW_data$alt_hom_F*geneticLoad_indv_NARW_data$s + 
  0.5*geneticLoad_indv_NARW_data$s*geneticLoad_indv_NARW_data$het_F 

### calculate realized genetic load ###
geneticLoad_indv_NARW_data$realizedLoad <- geneticLoad_indv_NARW_data$alt_hom_F*geneticLoad_indv_NARW_data$s #+ geneticLoad_indv_NARW_data$het_F*geneticLoad_indv_NARW_data$h*geneticLoad_indv_NARW_data$s


### calculate masked genetic load ###
geneticLoad_indv_NARW_data$MaskedLoad <- geneticLoad_indv_NARW_data$het_F*geneticLoad_indv_NARW_data$s*0.5

### remove modifier sites ###
geneticLoad_indv_NARW_data_NS <- geneticLoad_indv_NARW_data[which(geneticLoad_indv_NARW_data$impact != "Modifier"),]







head(geneticLoad_indv_NARW_data_NS)
######### Total mutation load #####################################################################
## aggregate to get genome-wide estimate of Realized genetic load ##
agg_NARW_geneticLoad_indv_data_NS <- aggregate(geneticLoad_indv_NARW_data_NS$TotalLoad,  
                                               by=list(geneticLoad_indv_NARW_data_NS$sample,
                                                       geneticLoad_indv_NARW_data_NS$CHR,
                                                       geneticLoad_indv_NARW_data_NS$impact,
                                                       geneticLoad_indv_NARW_data_NS$population), FUN=mean)

agg_NARW_geneticLoad_indv_data_NS

NARW_geneticLoad_indv_data_NS_loads <- agg_NARW_geneticLoad_indv_data_NS %>% spread(Group.3, x)
NARW_geneticLoad_indv_data_NS_loads$sumLoad <- (NARW_geneticLoad_indv_data_NS_loads$High + 
                                                  NARW_geneticLoad_indv_data_NS_loads$Moderate +  
                                                  NARW_geneticLoad_indv_data_NS_loads$Low)


se <- function(x) sd(x) / sqrt(length(x))

NARW_sum_Totalload_summary <- NARW_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

NARW_sum_Totalload_summary

NARW.aov.total.sum <- aov(sumLoad ~ Group.4*Group.2, data = NARW_geneticLoad_indv_data_NS_loads)
summary(NARW.aov.total.sum)

#Df   Sum Sq  Mean Sq F value   Pr(>F)    
#Group.4          1 0.001652 0.001652   2.693   0.1203    
#Group.2          1 0.023677 0.023677  38.598 1.24e-05 ***
#  Group.4:Group.2  1 0.003670 0.003670   5.982   0.0264 *  
#  Residuals       16 0.009815 0.000613  



head(geneticLoad_indv_NARW_data_NS)
######### Realized mutation load #####################################################################
## aggregate to get genome-wide estimate of Realized genetic load ##
agg_NARW_geneticLoad_indv_data_NS <- aggregate(geneticLoad_indv_NARW_data_NS$realizedLoad,  
                                               by=list(geneticLoad_indv_NARW_data_NS$sample,
                                                       geneticLoad_indv_NARW_data_NS$CHR,
                                                       geneticLoad_indv_NARW_data_NS$impact,
                                                       geneticLoad_indv_NARW_data_NS$population), FUN=median)

NARW_geneticLoad_indv_data_NS_loads <- agg_NARW_geneticLoad_indv_data_NS %>% spread(Group.3, x)
NARW_geneticLoad_indv_data_NS_loads$sumLoad <- (NARW_geneticLoad_indv_data_NS_loads$High + 
                                                  NARW_geneticLoad_indv_data_NS_loads$Moderate +  
                                                  NARW_geneticLoad_indv_data_NS_loads$Low)


se <- function(x) sd(x) / sqrt(length(x))

NARW_sum_realizedload_summary <- NARW_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

NARW_sum_realizedload_summary

NARW.aov.realized.sum <- aov(sumLoad ~ Group.4*Group.2, data = NARW_geneticLoad_indv_data_NS_loads)
summary(NARW.aov.realized.sum)

#Df  Sum Sq Mean Sq F value   Pr(>F)    
# Group.4          1 0.00020 0.00020   0.110    0.745    
# Group.2          1 0.05534 0.05534  29.692 5.35e-05 ***
#  Group.4:Group.2  1 0.00064 0.00064   0.344    0.566    
# Residuals       16 0.02982 0.00186   



head(geneticLoad_indv_NARW_data_NS)
######### Masked mutation load #####################################################################
## aggregate to get genome-wide estimate of Realized genetic load ##
agg_NARW_geneticLoad_indv_data_NS <- aggregate(geneticLoad_indv_NARW_data_NS$MaskedLoad,  
                                               by=list(geneticLoad_indv_NARW_data_NS$sample,
                                                       geneticLoad_indv_NARW_data_NS$CHR,
                                                       geneticLoad_indv_NARW_data_NS$impact,
                                                       geneticLoad_indv_NARW_data_NS$population), FUN=median)

NARW_geneticLoad_indv_data_NS_loads <- agg_NARW_geneticLoad_indv_data_NS %>% spread(Group.3, x)
NARW_geneticLoad_indv_data_NS_loads$sumLoad <- (NARW_geneticLoad_indv_data_NS_loads$High + 
                                                  NARW_geneticLoad_indv_data_NS_loads$Moderate +  
                                                  NARW_geneticLoad_indv_data_NS_loads$Low)


se <- function(x) sd(x) / sqrt(length(x))

NARW_sum_maskedload_summary <- NARW_geneticLoad_indv_data_NS_loads %>% 
  group_by(Group.4, Group.2) %>% 
  summarise(mean = mean(sumLoad),
            SE = se(sumLoad),
            SD = sd(sumLoad))

NARW_sum_maskedload_summary

NARW.aov.realized.sum <- aov(sumLoad ~ Group.4*Group.2, data = NARW_geneticLoad_indv_data_NS_loads)
summary(NARW.aov.realized.sum)

#Df  Sum Sq Mean Sq F value   Pr(>F)    
#Group.4          1 0.00096 0.00096   0.404    0.534    
#Group.2          1 0.21095 0.21095  88.275 6.49e-08 ***
#  Group.4:Group.2  1 0.00179 0.00179   0.749    0.400    
#Residuals       16 0.03823 0.00239  


##################################### PLOT ###############################################################
scaleFUN <- function(x) sprintf("%.2f", x)
b1 <- ggplot((NARW_sum_Totalload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 11, color = "black")) +
  scale_y_continuous(labels=scaleFUN) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  ylab("Total mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none") 


b2 <- ggplot((NARW_sum_realizedload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 11, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  scale_y_continuous(labels=scaleFUN) +
  ylab("Realized mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  

scaleFUN <- function(x) sprintf("%.2f", x)
b3 <- ggplot((NARW_sum_maskedload_summary), aes(x=Group.2, y=mean, fill = Group.4)) + 
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(x=Group.2, ymin=mean-SE, ymax=mean+SE), width=0.4, colour="black", alpha=0.9, size=.5, position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c("palevioletred4", "peru")) +
  theme(axis.text = element_text(size = 11, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black")) +
  scale_y_continuous(labels=scaleFUN) +
  ylab("Masked mutation load") +
  xlab("Genomic region") +
  theme(legend.position = "none")  


Figure_4 <- (b1|b2|b3) +
  plot_annotation(tag_levels = 'A') #add figure labels
Figure_4 #view multi-panel figure








######################################################################
#                                                                    #
#    gene set enrichment with Getsalt                                #
#                 REALIZED female burden                             #
######################################################################



enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                            enrichDatabase="network_PPI_BIOGRID", interestGeneFile="/Volumes/cetacea/Genetic_load/_enrichment/gene.prelim.txt",
                            interestGeneType="genesymbol", sigMethod="top", topThr=10,
                            outputDirectory=getwd(), highlightSeedNum=10,
                            networkConstructionMethod="Network_Retrieval_Prioritization")

listOrganism()


refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
outputDirectory <- "/Volumes/cetacea/Genetic_load/_enrichment/"
enrichResult <- WebGestaltR(enrichMethod="NTA", organism="btaurus",
                            enrichDatabase="network_PPI_BIOGRID", interestGeneFile="/Volumes/cetacea/Genetic_load/_enrichment/gene.prelim.txt",
                            interestGeneType="genesymbol", sigMethod="top", topThr=10,
                            outputDirectory="/Volumes/cetacea/Genetic_load/_enrichment/", 
                            highlightSeedNum=10,
                            isOutpu = T,
                            gseaPlotFormat = c("png", "svg"),
                            networkConstructionMethod="Network_Expansion")


######### plotting results #########3


ORA_bp <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1690904587/goslim_summary_wg_result1690904587_bp.txt", 
                     header = F,
                     sep="\t")


eone <- ORA_bp %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "navy")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("navy")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 8, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, color = "black")) +
  ylab("Count (Biological prosesses)") +
  #xlab("Biological process") +
  theme(legend.position = "none")  + 
  scale_y_continuous(breaks= pretty_breaks())





ORA_cc <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1690904587/goslim_summary_wg_result1690904587_cc.txt", 
                     header = F,
                     sep="\t")


etwo <- ORA_cc %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "navy")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("dodgerblue4")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 8, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, color = "black")) +
  ylab("Count (Cellular components)") +
  #xlab("Biological process") +
  theme(legend.position = "none") + 
  scale_y_continuous(breaks= pretty_breaks())




ORA_mf <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1690904587/goslim_summary_wg_result1690904587_mf.txt", 
                     header = F,
                     sep="\t")


ethree <- ORA_mf %>%
  mutate(V2 = fct_reorder(V2, desc(V3))) %>%
  ggplot(aes(x=V2, y=V3, fill = "navy")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("slateblue")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  #theme(axis.text.x = element_text(size = 10, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, size = 8, color = "black", hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, color = "black")) +
  ylab("Count (Molecular functions)") +
  #xlab("Biological process") +
  scale_y_continuous(breaks= pretty_breaks()) +
  theme(legend.position = "none")  





ORA_enrichment <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/Project_wg_result1690904587/enrichment_results_wg_result1690904587.txt", 
                             header = T,
                             sep="\t")

efour <- ORA_enrichment %>%
  mutate(description = fct_reorder(description, enrichmentRatio)) %>%
  ggplot(aes(x=enrichmentRatio, y=description, fill = "darkgoldenrod")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("darkgoldenrod")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 10, color = "black")) +
  theme(axis.title = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_blank()) +
  #ylab("Gene Ontology (ORA)") +
  ggtitle("Gene Ontology (Over representation)") +
  xlab("Enrichment ratio") +
  theme(legend.position = "none")  +
  theme(plot.title.position = "plot")
#theme(axis.title.y = element_text(angle = 0, hjust = 0, margin = margin(r = -125)))






NTA_enrichment <- read.table("/Volumes/cetacea/Genetic_load/_enrichment/wg_result1690922182/wg_result1690922182.network_PPI_BIOGRID.Network_Expansion_enrichedResult.txt", 
                             header = T,
                             sep="\t")


efive <- NTA_enrichment %>%
  mutate(description = fct_reorder(description, enrichmentRatio)) %>%
  ggplot(aes(x=enrichmentRatio, y=description, fill = "darkgolenrod2")) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c("darkgoldenrod2")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  #theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 10, color = "black")) +
  theme(axis.title.x = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_blank()) +
  ggtitle("Gene Ontology (Network topology)") +
  #ylab("Gene Ontology (NTA)") +
  xlab("Enrichment ratio") +
  theme(legend.position = "none") +
  theme(plot.title.position = "plot")
#theme(axis.title.y = element_text(angle = 0, hjust = 0, margin = margin(r = -125)))





Figure_5 <- (((eone | etwo | ethree) ) | (efour / efive)) + 
  plot_layout(ncol=4,widths =c(3,3,2.5,2)) +
  plot_annotation(tag_levels = 'A') #add figure labels
Figure_5 #view multi-panel figure












#########################################
#### Rxy (per site) version two #########
#########################################

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


inter_site_load_data_B <- melt(data = inter_site_load_data, id.vars = c("Chromosome", "POS", "Population", "Impact"), 
                   variable.name = "Estimate", 
                   variable.value = "value")

head(inter_site_load_data_B)
tail(inter_site_load_data_B)

inter_site_load_data_C <- inter_site_load_data_B[which(inter_site_load_data_B$Estimate == "alternate_allele_frequency"),]
head(inter_site_load_data_C)

inter_site_load_data_D <- subset(inter_site_load_data_C, select = -c(Estimate, POS))
head(inter_site_load_data_D)

inter_site_load_data_E <- inter_site_load_data_D %>% spread(Population, value)
head(inter_site_load_data_E)

agg_Rxy_data_B$NARW_SRW <- agg_Rxy_data_B$NARW/agg_Rxy_data_B$SRW
agg_Rxy_data_B$NARW_BH <- agg_Rxy_data_B$NARW/agg_Rxy_data_B$BH
agg_Rxy_data_B$SRW_BH <- agg_Rxy_data_B$SRW/agg_Rxy_data_B$BH


