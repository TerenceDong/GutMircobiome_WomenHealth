setwd("~/TRAINING/curatedMetagenomicData/disease_code_recheck")
source("code_library.R")
library(dplyr)
library(openxlsx)

cohortID <- read.xlsx("CohortID.xlsx",sheet = 1)
individualID <- read.xlsx("CohortID.xlsx",sheet = 2)
ageID <- read.xlsx("CohortID.xlsx",sheet = 3)

IBDCountryCohort <- as.character(cohortID$IBDCountryCohort[!is.na(cohortID$IBDCountryCohort)])
IBDContinentCohort <- as.character(cohortID$IBDContinentCohort[!is.na(cohortID$IBDContinentCohort)])
T2DCountryCohort <- as.character(cohortID$T2DCountryCohort[!is.na(cohortID$T2DCountryCohort)])
T2DContinentCohort <- as.character(cohortID$T2DContinentCohort[!is.na(cohortID$T2DContinentCohort)])
AVCDCountryCohort <- as.character(cohortID$AVCDCountryCohort[!is.na(cohortID$AVCDCountryCohort)])
AVCDContinentCohort <- as.character(cohortID$AVCDContinentCohort[!is.na(cohortID$AVCDContinentCohort)])
IGTCountryCohort <- as.character(cohortID$IGTCountryCohort[!is.na(cohortID$IGTCountryCohort)])
IGTContinentCohort <- as.character(cohortID$IGTContinentCohort[!is.na(cohortID$IGTContinentCohort)])

IBDIndividuals <- as.character(individualID$IBDIndividuals[!is.na(individualID$IBDIndividuals)])
T2DIndividuals <- as.character(individualID$T2DIndividuals[!is.na(individualID$T2DIndividuals)])
AVCDIndividuals <- as.character(individualID$AVCDIndividuals[!is.na(individualID$AVCDIndividuals)])
IGTIndividuals <- as.character(individualID$IGTIndividuals[!is.na(individualID$IGTIndividuals)])
SelectControls <- as.character(individualID$SelectControls[!is.na(individualID$SelectControls)])

child_adolescent <- as.character(ageID$child_adolescent[!is.na(ageID$child_adolescent)])
young <- as.character(ageID$young[!is.na(ageID$young)])
middle <- as.character(ageID$middle[!is.na(ageID$middle)])
old <- as.character(ageID$old[!is.na(ageID$old)])

CountryCohort <- unique(c(IBDCountryCohort,T2DCountryCohort,AVCDCountryCohort,IGTCountryCohort))

# batch effect before adjusting------------------------------------------------------------

d1_disease <- read.xlsx("metadata_disease.xlsx",rowNames = T)
d2_disease <- read.csv("RelativeAbundance_disease.csv",row.names = 1)
d1_disease <- d1_disease[CountryCohort,]
d2_disease <- d2_disease[CountryCohort,]

d1_disease$DNA_extraction_kit[is.na(d1_disease$DNA_extraction_kit)==TRUE] <- "na" 
d1_disease$sample_id <- rownames(d1_disease)
d1_disease$sample_id <- gsub("\\.","-",d1_disease$sample_id)

d2_reAbund_feature <- as.data.frame(matrix(NA,ncol = 3))
colnames(d2_reAbund_feature) <- c("species","prevalence","proportion")
Species <- colnames(d2_disease)
for (i in 1:ncol(d2_disease)){
  d2_reAbund_feature[i,1] = Species[i]
  d2_reAbund_feature[i,2] = 1-sum(d2_disease[,Species[i]]==0)/nrow(d2_disease)
  d2_reAbund_feature[i,3] = sum(d2_disease[,Species[i]])/sum(d2_disease)
}
d2_reAbund_select <- d2_reAbund_feature$species[d2_reAbund_feature$prevalence >= 0.1 & d2_reAbund_feature$proportion >= 0.00005]

d2_disease_select <- d2_disease[,d2_reAbund_select]
d2_disease_select$sample_id <- row.names(d2_disease_select)
d2_disease_select$sample_id <- gsub("\\.","-",d2_disease_select$sample_id)
d2_disease_select <- merge(d2_disease_select,d1_disease,by="sample_id",all=FALSE)

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(d2_disease_select[,d2_reAbund_select],method = "bray",diag = TRUE,upper = TRUE))

pcoa <- cmdscale(bray_dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) <- c("x", "y", "z") 
eig <- pcoa$eig

library(reshape2)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- detectCores()
registerDoParallel(cl)

library(patchwork)
library(ape) 
library(vegan) 
library(ggsci)
library(ggplot2) 
library(openxlsx)
library(ggpubr)

# sequencing_platform, Figure S4c -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select$sequencing_platform)
colnames(points_2)[4] <- "sequencing_platform"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ sequencing_platform, data = d2_disease_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ sequencing_platform, data = d2_disease_select,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=sequencing_platform)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=sequencing_platform), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#db5b38","#8eb853","#5f9ccb"))
ggsave("bray_sequencing_platform.pdf", p, width = 8, height = 6)

# DNA_extraction_kit, Figure S4a -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select$DNA_extraction_kit)
colnames(points_2)[4] <- "DNA_extraction_kit"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ DNA_extraction_kit, data = d2_disease_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ DNA_extraction_kit, data = d2_disease_select,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=DNA_extraction_kit)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=DNA_extraction_kit), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#91bde2", "#83c5b0", "#8fb79d","#eadfc3", "#eec16c","#d75032", "#ce9c83","#b9b7a2", "#f66f69","#289ebc", "#e8edf3","#8171b5","#4677b7"))
ggsave("bray_DNA_extraction_kit.pdf", p, width = 8, height = 6)

# median_read_length, Figure S4b -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select$median_read_length)
colnames(points_2)[4] <- "median_read_length"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ median_read_length, data = d2_disease_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ median_read_length, data = d2_disease_select,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=median_read_length)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=median_read_length), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_gradientn(colors=colorRampPalette(colors = c("#f2f6f9","#d7e3ee","#a4bfd9","#3372ab"))(500))
ggsave("bray_median_read_length.pdf", p, width = 8, height = 6)

# batch effect calibration------------------------------------------------------------

# BiocManager::install("MMUPHin")
library(MMUPHin)

d2_disease_select2 <- d2_disease_select[,d2_reAbund_select]
rownames(d2_disease_select2) <- d2_disease_select$sample_id
d2_disease_select2 <- t(d2_disease_select2) # feature-by-sample
d2_disease_select2 <- d2_disease_select2/100 # divided by 100 to be compositional data

d2_disease_select2_meta <- d2_disease_select[,colnames(d1_disease)] # sample-by-metadata
rownames(d2_disease_select2_meta) <- d2_disease_select$sample_id

d2_disease_select2_meta$DNA_extraction_kit <- factor(d2_disease_select2_meta$DNA_extraction_kit) # assume the variable adjusted to be factor variable
d2_disease_select2_adjust1 <- adjust_batch(feature_abd = d2_disease_select2,
                                           batch = c("DNA_extraction_kit"),
                                           # covariates = "study_condition",
                                           data = d2_disease_select2_meta)$feature_abd_adj

d2_disease_select2_meta$sequencing_platform <- factor(d2_disease_select2_meta$sequencing_platform) # assume the variable adjusted to be factor variable
d2_disease_select2_adjust2 <- adjust_batch(feature_abd = d2_disease_select2_adjust1,
                                           batch = c("sequencing_platform"),
                                           # covariates = "study_condition",
                                           data = d2_disease_select2_meta)$feature_abd_adj

d2_disease_select2_meta$median_read_length <- factor(d2_disease_select2_meta$median_read_length) # assume the variable adjusted to be factor variable
d2_disease_select2_adjust3 <- adjust_batch(feature_abd = d2_disease_select2_adjust2,
                                           batch = c("median_read_length"),
                                           # covariates = "study_condition",
                                           data = d2_disease_select2_meta)$feature_abd_adj

# batch effect after adjusting------------------------------------------------------------

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(t(d2_disease_select2_adjust3),method = "bray",diag = TRUE,upper = TRUE))

pcoa <- cmdscale(bray_dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) <- c("x", "y", "z") 
eig <- pcoa$eig

library(reshape2)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- detectCores()
registerDoParallel(cl)

# sequencing_platform, Figure S4c -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select2_meta$sequencing_platform)
colnames(points_2)[4] <- "sequencing_platform"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ sequencing_platform, data = d2_disease_select2_meta,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ sequencing_platform, data = d2_disease_select2_meta,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=sequencing_platform)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=sequencing_platform), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#db5b38","#8eb853","#5f9ccb"))
ggsave("bray_sequencing_platform_after.pdf", p, width = 8, height = 6)

# DNA_extraction_kit, Figure S4a -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select2_meta$DNA_extraction_kit)
colnames(points_2)[4] <- "DNA_extraction_kit"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ DNA_extraction_kit, data = d2_disease_select2_meta,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ DNA_extraction_kit, data = d2_disease_select2_meta,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=DNA_extraction_kit)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=DNA_extraction_kit), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_manual(values = c("#91bde2", "#83c5b0", "#8fb79d","#eadfc3", "#eec16c","#d75032", "#ce9c83","#b9b7a2", "#f66f69","#289ebc", "#e8edf3","#8171b5","#4677b7"))
ggsave("bray_DNA_extraction_kit_after.pdf", p, width = 8, height = 6)


# median_read_length, Figure S4b -----------------------------------------------------------
points_2 <- cbind(points, d2_disease_select$median_read_length)
colnames(points_2)[4] <- "median_read_length"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ median_read_length, data = d2_disease_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ median_read_length, data = d2_disease_select,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "")

p <- ggplot(points_2, aes(x=x, y=y, color=median_read_length)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=median_read_length), alpha = 0.1, show.legend = F) +
  theme_bw()+
  scale_color_gradientn(colors=colorRampPalette(colors = c("#f2f6f9","#d7e3ee","#a4bfd9","#3372ab"))(500))
ggsave("bray_median_read_length_after.pdf", p, width = 8, height = 6)


d2_disease_select2_adjust3_t <- data.frame(t(d2_disease_select2_adjust3))
d2_disease_select2_adjust3_t$sample_id <- rownames(d2_disease_select2_adjust3_t)
d2_disease_select2_adjust3_t_merge <- merge(d2_disease_select2_adjust3_t,d2_disease_select2_meta,by="sample_id",all = T)
write.xlsx(d2_disease_select2_adjust3_t_merge,"RelAbund_batchCorrect.xlsx")

# conclusion: DNA_extraction_kit were used as covariates to curb their impacts on microbial abundance

# following analysis of diseased women ----------------------------------

dat_disease <- read.xlsx("RelAbund_batchCorrect.xlsx")
continent <- data.frame(
  country=c("CAN","CHN","DEU","DNK","EST","FIN","FRA","GBR","IND","IRL","ISL","ITA","LUX","MDG","NLD","RUS","SWE","TZA","USA","ESP","HUN","NOR","SVK"),
  continent=c("North America","Asia","Europe","Europe","Europe","Europe","Europe","Europe","Asia","Europe","Europe","Europe","Europe","Africa","Europe","Europe","Europe","Asia","North America","Europe","Europe","Europe","Europe")
)
dat_disease <- merge(dat_disease,continent,by="country") # use merge function to simulate vlookup
bac_species <- c(3:161)

short_phylum <- read.xlsx("short_phylum_disease.xlsx",sheet = 1)
phylum <- unique(short_phylum$phylum)
for (i in phylum){
  dat=subset(short_phylum,phylum==i)
  dat_disease[,i]=tryCatch(rowSums(dat_disease[,dat$short]),error=function(e){"error"})
}
dat_disease$Euryarchaeota <- dat_disease$Methanobrevibacter.smithii
dat_disease$Verrucomicrobia <- dat_disease$Akkermansia.muciniphila
bac_phylum <- c(247:252)

dat_disease$age_group <- ifelse(dat_disease$age <=1,"infant",
                                ifelse(dat_disease$age >=2 & dat_disease$age <=17,"child_adolescent",
                                       ifelse(dat_disease$age >=18 & dat_disease$age <=39,"young",
                                              ifelse(dat_disease$age >=40 & dat_disease$age <=59,"middle",
                                                     ifelse(dat_disease$age >=60,"old",NA)))))
name_all <- colnames(dat_disease)
rownames(dat_disease) <- dat_disease$sample_id
bac_species_name <- name_all[bac_species]
bac_phylum_name <- name_all[bac_phylum]

# Table S6 ----------------------------------------------------------------

dat_disease_2 <- subset(dat_disease,disease!="healthy")
group_by(dat_disease_2[,c("age","BMI","country","disease")],disease,country) %>% summarise_each(funs(length))

dat_disease_2 <- subset(dat_disease,!is.na(age)&disease!="healthy")
group_by(dat_disease_2[,c("age","country","disease")],disease,country) %>% summarise_each(funs(mean,sd))

dat_disease_2 <- subset(dat_disease,!is.na(BMI)&disease!="healthy")
group_by(dat_disease_2[,c("BMI","country","disease")],disease,country) %>% summarise_each(funs(mean,sd))

dat_disease_2 <- subset(dat_disease,disease=="healthy")
group_by(dat_disease_2[,c("age","BMI","country","disease")],disease,country) %>% summarise_each(funs(length))

dat_disease_2 <- subset(dat_disease,!is.na(age)&disease=="healthy")
group_by(dat_disease_2[,c("age","country","disease")],disease,country) %>% summarise_each(funs(mean,sd))

dat_disease_2 <- subset(dat_disease,!is.na(BMI)&disease=="healthy")
group_by(dat_disease_2[,c("BMI","country","disease")],disease,country) %>% summarise_each(funs(mean,sd))
# Westernization和PMID信息在RelAbund_batchCorrect.xlsx对应non_westernized、PMID两列，统计分析手工进行

# proportion of microbial variation among diseased and healthy women in selected population ----------------------------------

dat_disease_2 <- dat_disease[CountryCohort,]

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(dat_disease_2[,name_all[bac_species]],method = "bray",diag = TRUE,upper = TRUE))

pcoa <- cmdscale(bray_dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) <- c("x", "y", "z") 
eig <- pcoa$eig

library(reshape2)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- detectCores()
registerDoParallel(cl)

# country -----------------------------------------------------------
points_2 <- cbind(points, dat_disease_2$country)
colnames(points_2)[4] <- "country"

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ country+DNA_extraction_kit, data = dat_disease_2,permutations = 999)}
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.085 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=country)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=country), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("bray_country_after_selected.pdf", p, width = 8, height = 6)

# continent -----------------------------------------------------------
points_2 <- cbind(points, dat_disease_2$continent)
colnames(points_2)[4] <- "continent"

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ continent+DNA_extraction_kit, data = dat_disease_2,permutations = 999)}
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.019 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=continent)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=continent), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("bray_continent_after_selected.pdf", p, width = 8, height = 6)

# age -----------------------------------------------------------
points_2 <- cbind(points, dat_disease_2$age)
colnames(points_2)[4] <- "age"
points_2 <- subset(points_2,!is.na(age))
idx <- rownames(bray_dis) %in% rownames(points_2)

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis[idx,idx] ~ age+DNA_extraction_kit, data = dat_disease_2[idx,],permutations = 999)}
# group_permanova <- adonis(bray_dis ~ age, data = dat_disease_2,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.072 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=age)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=age), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("bray_age_after_selected.pdf", p, width = 8, height = 6)

# BMI -----------------------------------------------------------
points_2 <- cbind(points, dat_disease_2$BMI)
colnames(points_2)[4] <- "BMI"
points_2 <- subset(points_2,!is.na(BMI))
idx <- rownames(bray_dis) %in% rownames(points_2)

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis[idx,idx] ~ BMI+DNA_extraction_kit, data = dat_disease_2[idx,],permutations = 999)}
# group_permanova <- adonis(bray_dis ~ BMI, data = dat_disease_2,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.004 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=BMI)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=BMI), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("bray_BMI_after_select.pdf", p, width = 8, height = 6)

# summary of country, continent, BMI and age, Figure S5a----------------------------------------
result_barchart <- read.xlsx("result_barchart.xlsx",sheet = 1) # including R2 and p from above PERMANOVA analysis
result_barchart$Variable <- factor(result_barchart$Variable,levels = c("BMI","Age","Continent","Country"))
p <- ggplot(result_barchart,aes(Variable,R2))+
  geom_col(aes(fill=Variable,color=Variable))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 10))+
  labs(x="",y="R2")+
  scale_color_manual(values=rev(c("#336ba2","#33adc3","#f13333","#a87eb2")))+
  scale_fill_manual(values=rev(c("#336ba2","#33adc3","#f13333","#a87eb2")))+
  coord_flip()
p
ggsave("DiseasedCohort_R2.pdf",width = 4,height = 4)

# bootstrapped_PERMANOVA, Figure S5b --------------------------------------------------

bootstrap_PERMANOVA <- function(data,metadata,variable,iter,prob){
  library(foreach)
  library(doParallel)
  library(vegan)
  cl <- detectCores()
  registerDoParallel(cl)
  results=as.data.frame(matrix(NA,ncol = 3))
  colnames(results)=c("ID","r2","p")
  for (i in 1:iter){
    set.seed(i)
    idx=sample(c(0,1),nrow(data),replace = T,prob = prob)
    data_sample=data[idx==0,idx==0]
    metadata_sample=metadata[idx==0,]
    variable_select=metadata_sample[,variable]
    set.seed(i)
    group_permanova <- foreach(
      n=1,.combine=cbind,.packages = 'vegan') %dopar% {
        adonis(data_sample ~ variable_select+DNA_extraction_kit, data = metadata_sample,permutations = 999)
      }
    results[i,1]=i
    results[i,2]=round(group_permanova$aov.tab$R2[1],digits = 3)
    results[i,3]=group_permanova$aov.tab$`Pr(>F)`[1]
  }
  return(results)
}

result_country <- bootstrap_PERMANOVA(bray_dis,dat_disease_2,"country",100,c(0.2,0.8))
write.xlsx(result_country,"result_country.xlsx")

result_continent <- bootstrap_PERMANOVA(bray_dis,dat_disease_2,"continent",100,c(0.2,0.8))
write.xlsx(result_continent,"result_continent.xlsx")

points_2 <- cbind(points, dat_disease_2$age)
colnames(points_2)[4] <- "age"
points_2 <- subset(points_2,!is.na(age))
idx <- rownames(bray_dis) %in% rownames(points_2)
result_age <- bootstrap_PERMANOVA(bray_dis[idx,idx],dat_disease_2[idx,],"age",100,c(0.2,0.8))
write.xlsx(result_age,"result_age.xlsx")

points_2 <- cbind(points, dat_disease_2$BMI)
colnames(points_2)[4] <- "BMI"
points_2 <- subset(points_2,!is.na(BMI))
idx <- rownames(bray_dis) %in% rownames(points_2)
result_BMI <- bootstrap_PERMANOVA(bray_dis[idx,idx],dat_disease_2[idx,],"BMI",100,c(0.2,0.8))
write.xlsx(result_BMI,"result_BMI.xlsx")

bootstrapped_permanova <- read.xlsx("bootstrapped_permanova_disease.xlsx",sheet = 1)
bootstrapped_permanova_reshape <- melt(bootstrapped_permanova,id.vars="ID",variable.name = "group",value.name = "r2")
bootstrapped_permanova_reshape$group <- factor(bootstrapped_permanova_reshape$group,levels = c("Country","Continent","Age","BMI"))

p <- ggplot(bootstrapped_permanova_reshape,aes(group,r2))+
  stat_boxplot(geom = "errorbar",aes(color=group),width=0.2)+
  geom_boxplot(aes(color=group),outlier.alpha = 0,width=0.6)+
  geom_jitter(aes(fill=group),width = 0.2,size=1.4,shape=21,color="black",alpha=0.6)+
  theme_bw()+
  scale_color_manual(values = c("#336ba2","#33adc3","#f13333","#a87eb2"))+
  scale_fill_manual(values = c("#336ba2","#33adc3","#f13333","#a87eb2"))+
  theme(axis.text = element_text(color = "black",size = 8))+
  labs(x="",y="R2")+
  theme(axis.title = element_text(color="black",size=15),
        axis.text.x = element_text(color = "black",size = 12,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(color = "black",size = 12))
ggsave("bootstrapped_permanova_disease_selected.pdf",width=6,height=7)
# geom_jitter添加抖动效果，每次成图时点的水平位置会有随机性差异，但垂直位置保持一致

# iterative RF disease country cohort, Figure 4a & Figure S6a------------------------------------------------------------
library(compositions)
set.seed(2023)
idx <- sample(c(0,1),length(IBDCountryCohort),replace = T,prob = c(0.5,0.5))
IBDCountryCohortTrain <- IBDCountryCohort[idx==0]
IBDCountryCohortTest <- IBDCountryCohort[idx==1]
rfPermutedMultipleNull_IBDCountryCohort <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                        IBDCountryCohortTrain,IBDCountryCohortTest,IBDIndividuals,SelectControls,200,10,round(0.3*length(IBDIndividuals),digits = 0),round(0.3*length(IBDIndividuals),digits = 0))
# trainsize和testsize保持一致，且为训练集或验证集中样本量最小集合中关注结局（此处为Disease）样本量的0.6

set.seed(2023)
idx <- sample(c(0,1),length(T2DCountryCohort),replace = T,prob = c(0.5,0.5))
T2DCountryCohortTrain <- T2DCountryCohort[idx==0]
T2DCountryCohortTest <- T2DCountryCohort[idx==1]
rfPermutedMultipleNull_T2DCountryCohort <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                        T2DCountryCohortTrain,T2DCountryCohortTest,T2DIndividuals,SelectControls,200,10,round(0.3*length(T2DIndividuals),digits = 0),round(0.3*length(T2DIndividuals),digits = 0))

set.seed(2023)
idx <- sample(c(0,1),length(AVCDCountryCohort),replace = T,prob = c(0.5,0.5))
AVCDCountryCohortTrain <- AVCDCountryCohort[idx==0]
AVCDCountryCohortTest <- AVCDCountryCohort[idx==1]
rfPermutedMultipleNull_AVCDCountryCohort <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                         AVCDCountryCohortTrain,AVCDCountryCohortTest,AVCDIndividuals,SelectControls,200,10,round(0.3*length(AVCDIndividuals),digits = 0),round(0.3*length(AVCDIndividuals),digits = 0))

set.seed(2023)
idx <- sample(c(0,1),length(IGTCountryCohort),replace = T,prob = c(0.5,0.5))
IGTCountryCohortTrain <- IGTCountryCohort[idx==0]
IGTCountryCohortTest <- IGTCountryCohort[idx==1]
rfPermutedMultipleNull_IGTCountryCohort <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                        IGTCountryCohortTrain,IGTCountryCohortTest,IGTIndividuals,SelectControls,200,10,round(0.3*length(IGTIndividuals),digits = 0),round(0.3*length(IGTIndividuals),digits = 0))

# IBDCountryCohort 
rfTestAUCs_IBDCountryCohort <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_IBDCountryCohort) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_IBDCountryCohort$AUC <- c(apply(rfPermutedMultipleNull_IBDCountryCohort$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IBDCountryCohort$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_IBDCountryCohort$AUCNull1,1,median))
rfTestAUCs_IBDCountryCohort$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_IBDCountryCohort$Group <- c(rep("IBD",600))
rfTestAUCs_IBDCountryCohort$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# T2DCountryCohort 
rfTestAUCs_T2DCountryCohort <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_T2DCountryCohort) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_T2DCountryCohort$AUC <- c(apply(rfPermutedMultipleNull_T2DCountryCohort$AUCSameAge,1,median),apply(rfPermutedMultipleNull_T2DCountryCohort$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_T2DCountryCohort$AUCNull1,1,median))
rfTestAUCs_T2DCountryCohort$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_T2DCountryCohort$Group <- c(rep("T2D",600))
rfTestAUCs_T2DCountryCohort$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# AVCDCountryCohort 
rfTestAUCs_AVCDCountryCohort <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_AVCDCountryCohort) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_AVCDCountryCohort$AUC <- c(apply(rfPermutedMultipleNull_AVCDCountryCohort$AUCSameAge,1,median),apply(rfPermutedMultipleNull_AVCDCountryCohort$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_AVCDCountryCohort$AUCNull1,1,median))
rfTestAUCs_AVCDCountryCohort$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_AVCDCountryCohort$Group <- c(rep("AVCD",600))
rfTestAUCs_AVCDCountryCohort$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# IGTCountryCohort 
rfTestAUCs_IGTCountryCohort <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_IGTCountryCohort) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_IGTCountryCohort$AUC <- c(apply(rfPermutedMultipleNull_IGTCountryCohort$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IGTCountryCohort$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_IGTCountryCohort$AUCNull1,1,median))
rfTestAUCs_IGTCountryCohort$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_IGTCountryCohort$Group <- c(rep("IGT",600))
rfTestAUCs_IGTCountryCohort$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# merge
rfTestAUCs_CountryDisease <- rbind(rfTestAUCs_AVCDCountryCohort,rfTestAUCs_IBDCountryCohort,rfTestAUCs_IGTCountryCohort,rfTestAUCs_T2DCountryCohort)
write.csv(rfTestAUCs_CountryDisease,"rfTestAUCs_CountryDisease.csv",quote = F,row.names = F)
group_by(rfTestAUCs_CountryDisease[,c("AUC","Group","ClassificationGroup")],Group,ClassificationGroup) %>% summarise_each(funs(mean,sd))

# plot
library(reshape2)
library(ggforce)
library(dplyr)
library(ggpubr)

# Figure S6a---------------------------------------------------------------------
p <- ggplot(rfTestAUCs_CountryDisease)+
  # 箱线图：
  geom_boxplot(aes(x = ClassificationGroup, y = AUC,color=ClassificationGroup),
               width=0.7,position = position_dodge(0),size=0.5,outlier.size = 0) + 
  # 散点：
  geom_jitter(aes(x = ClassificationGroup, y = AUC,fill=ClassificationGroup),
              width=0.2,shape=21,colour="black",size=1)+
  # 设置填充颜色
  scale_fill_manual(values=c("#619cff","#00ba38","#f8766d"))+
  scale_color_manual(values=c("#619cff","#00ba38","#f8766d"))+
  # 添加显著性检验：
  # stat_compare_means(aes(x = ClassificationGroup, y = AUC),
  #                   label="p.signif", method = "anova")+
  # 分面：
  facet_grid(.~Group, labeller = label_both)+
  # X轴y轴标题：
  labs(y="AUC",x="")+
  # 设置主题：
  theme_minimal()+
  theme(strip.background = element_rect(fill="#eaeae0", color = "#dcddcf"),  # 分面颜色
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text.y = element_text(size=11,face="plain",color="black"),
        axis.text.x = element_blank(),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),
        # legend.position="none",
        # legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks.x=element_blank())
ggsave("boxplot_line_country_disease.pdf", height = 6, width = 14)
# geom_jitter添加抖动效果，每次成图时点的水平位置会有随机性差异，但垂直位置保持一致

# Figure 4a---------------------------------------------------------------------
p <- ggplot(rfTestAUCs_CountryDisease[rfTestAUCs_CountryDisease$Group=="IBD",])+
  # 箱线图：
  geom_boxplot(aes(x = ClassificationGroup, y = AUC,color=ClassificationGroup),
               width=0.7,position = position_dodge(0),size=0.5,outlier.size = 0) + 
  # 散点：
  geom_jitter(aes(x = ClassificationGroup, y = AUC,fill=ClassificationGroup),
              width=0.2,shape=21,colour="black",size=1)+
  # 设置填充颜色
  scale_fill_manual(values=c("#619cff","#00ba38","#f8766d"))+
  scale_color_manual(values=c("#619cff","#00ba38","#f8766d"))+
  # 添加显著性检验：
  # stat_compare_means(aes(x = ClassificationGroup, y = AUC),
  #                   label="p.signif", method = "anova")+
  # 分面：
  facet_grid(.~Group, labeller = label_both)+
  # X轴y轴标题：
  labs(y="AUC",x="")+
  # 设置主题：
  theme_minimal()+
  theme(strip.background = element_rect(fill="#eaeae0", color = "#dcddcf"),  # 分面颜色
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text.y = element_text(size=11,face="plain",color="black"),
        axis.text.x = element_blank(),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),
        # legend.position="none",
        # legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks.x=element_blank())
ggsave("boxplot_line_country_IBD.pdf", height = 6, width = 5)
# geom_jitter添加抖动效果，每次成图时点的水平位置会有随机性差异，但垂直位置保持一致

# threshold determination, Figure S6b ----------------------------------------------

rfThreshold <- function(data,name){
  threshold <- seq(0,0.99,0.01)
  rf_threshold <- as.data.frame(matrix(NA,ncol = 3))
  colnames(rf_threshold) <-c("threshold", paste0("mean_featureScore",name),paste0("bac_number",name))
  for (i in 1:length(threshold)){
    rfBacteria=names(which(rank_scale(colMeans(data$featureProfile))>=threshold[i]))
    rf_threshold[i,1] = threshold[i]
    rf_threshold[i,2] = mean(colMeans(data$featureProfile[,rfBacteria]))
    rf_threshold[i,3] = length(rfBacteria)
  }
  return(rf_threshold)
}

rfThreshold_AVCD <- rfThreshold(rfPermutedMultipleNull_AVCDCountryCohort,"AVCD")
rfThreshold_IBD <- rfThreshold(rfPermutedMultipleNull_IBDCountryCohort,"IBD")
rfThreshold_IGT <- rfThreshold(rfPermutedMultipleNull_IGTCountryCohort,"IGT")
rfThreshold_T2D <- rfThreshold(rfPermutedMultipleNull_T2DCountryCohort,"T2D")

rfThreshold_country <- Reduce(function(x,y) merge(x,y,by="threshold",all=T),list(rfThreshold_AVCD,rfThreshold_IBD,rfThreshold_IGT,rfThreshold_T2D))
rfThreshold_country_melt <- melt(rfThreshold_country[,c(1,2,4,6,8)],id.vars = "threshold",value.name = "mean",variable.name = "country")

sample_select <- rownames(rfThreshold_country_melt) %in% rownames(rfThreshold_country_melt[rfThreshold_country_melt$country=="mean_featureScoreIBD",])
p1 <- ggplot(data = rfThreshold_country_melt[sample_select,], aes(x = threshold, y = mean, color = country)) +
  geom_line(data = rfThreshold_country_melt[sample_select,], aes(x = threshold, y = mean),linewidth=1.5) +
  scale_color_manual(values = c("#7ea3b5"))+
  labs(x="Percentile Rank of Mean Decrease in Gini ",y="The Mean of Mean Decrease in Gini of Remaining Species")+
  scale_x_continuous(breaks = seq(0,1,0.1),
                     label = c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))+
  geom_vline(xintercept = 0.85,linetype="dashed",color="black")+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12))
ggsave("rfThreshold_country_IBD.pdf",width = 6,height = 6)

# species marker of disease

rfSpeciesMarker_IBD<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IBDCountryCohort$featureProfile))>=0.85))
rfSpeciesMarker_IGT<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IGTCountryCohort$featureProfile))>=0.85))
rfSpeciesMarker_T2D<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_T2DCountryCohort$featureProfile))>=0.85))
rfSpeciesMarker_country <- data.frame(rfSpeciesMarker_AVCD,rfSpeciesMarker_IBD,rfSpeciesMarker_IGT,rfSpeciesMarker_T2D)
write.csv(rfSpeciesMarker_country,"rfSpeciesMarker_country_disease.csv",quote = F,row.names = F)

# species marker comparisons

rfSpeciesMarker_country_unique <- unique(c(rfSpeciesMarker_AVCD,rfSpeciesMarker_IBD,rfSpeciesMarker_IGT,rfSpeciesMarker_T2D))
rfFeatureProfile_country <- rbind(rfPermutedMultipleNull_AVCDCountryCohort$featureProfile,rfPermutedMultipleNull_IBDCountryCohort$featureProfile,rfPermutedMultipleNull_IGTCountryCohort$featureProfile,rfPermutedMultipleNull_T2DCountryCohort$featureProfile)
rfMeanMarkerScores_country <- cbind(colMeans(rfFeatureProfile_country[1:200,]),colMeans(rfFeatureProfile_country[201:400,]),colMeans(rfFeatureProfile_country[401:600,]),colMeans(rfFeatureProfile_country[601:800,]))
colnames(rfMeanMarkerScores_country) <- c("MeanFeatureProfile_AVCDCountry","MeanFeatureProfile_IBDCountry","MeanFeatureProfile_IGTCountry","MeanFeatureProfile_T2DCountry")
rfComparison_country <- batch_dunns(rfFeatureProfile_country[,rfSpeciesMarker_country_unique],as.factor(c(rep("AVCD",200),rep("IBD",200),rep("IGT",200),rep("T2D",200))))
rfComparisonSummary_country <- cbind(rfMeanMarkerScores_country[rfSpeciesMarker_country_unique,],rfComparison_country)
rfComparisonSummary_country <- as.data.frame(rfComparisonSummary_country)
rfComparisonSummary_country$AVCD <- ifelse(rownames(rfComparisonSummary_country) %in% rfSpeciesMarker_AVCD,1,0)
rfComparisonSummary_country$IBD <- ifelse(rownames(rfComparisonSummary_country) %in% rfSpeciesMarker_IBD,1,0)
rfComparisonSummary_country$IGT <- ifelse(rownames(rfComparisonSummary_country) %in% rfSpeciesMarker_IGT,1,0)
rfComparisonSummary_country$T2D <- ifelse(rownames(rfComparisonSummary_country) %in% rfSpeciesMarker_T2D,1,0)
write.csv(rfComparisonSummary_country,"rfComparisonSummary_country_disease.csv",quote = F,row.names = T)

# bar plot, Figure 4b------------------------------------------------------------------------
IBD_species_mean_feature <- read.xlsx("IBD_species_mean_feature.xlsx",sheet = 1)
IBD_species_mean_feature$rfSpeciesMarker_IBD <- factor(IBD_species_mean_feature$rfSpeciesMarker_IBD,levels = IBD_species_mean_feature$rfSpeciesMarker_IBD)
p <- ggplot(IBD_species_mean_feature,aes(rfSpeciesMarker_IBD,MeanFeatureProfile_IBDCountry))+
  geom_col(color="#7ea3b5",fill="#7ea3b5")+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 10))+
  labs(x="",y="Average Mean Decrease in Gini")+
  coord_flip()
p
ggsave("IBD_species_mean_feature.pdf",width = 6,height = 6)


# IBD-specific species relative abundance comparisons across age groups, Table S7--------------------------------------------------------------
dat_disease_2 <- dat_disease[dat_disease$disease=="IBD",]
AgeGroupsComparisons <- as.data.frame(matrix(NA,ncol = 2))
colnames(AgeGroupsComparisons) <- c("species","kruskal_p")
for (i in 1:length(rfSpeciesMarker_IBD)){
  comparison=kruskal.test(dat_disease_2[,rfSpeciesMarker_IBD[i]]~age_group,data = dat_disease_2)
  AgeGroupsComparisons[i,1]=rfSpeciesMarker_IBD[i]
  AgeGroupsComparisons[i,2]=comparison$p.value
}
AgeGroupsComparisons$p.adjust <- p.adjust(AgeGroupsComparisons$kruskal_p,method = "BH")
write.csv(AgeGroupsComparisons,"AgeGroupsComparisonsIBD.csv",row.names = F,quote = F)

# SHAP analysis, Figure 4c,d,e --------------------------------------------

IBD <- read.xlsx("shap_summary.xlsx",sheet = 1)
AVCD <- read.xlsx("shap_summary.xlsx",sheet = 2)
T2D <- read.xlsx("shap_summary.xlsx",sheet = 3)
IGT <- read.xlsx("shap_summary.xlsx",sheet = 4)

# age - disease probability, Figure 4c, d -----------------------------------------------
library(ggstatsplot)

# IBD, Figure 4d
glmModel <- glm(probability ~ age,data = IBD,family = "gaussian")
summary(glmModel)
p <- ggscatterstats(
  data = IBD,
  x = age,
  y = probability,
  type = "nonparametric", # type of test that needs to be run
  conf.level = 0.95, # confidence level
  xlab = "Age (year)", # label for x axis
  ylab = "IBD probability", # label for y axis
  ggtheme = theme_bw(), # choosing a different theme
  point.args = list(size = 3, alpha = 0.4, stroke = 0, color="#7ea3b5"),
  smooth.line.args = list(size = 1.5, color = "#7ea3b5", method = "lm", formula = y ~ x),
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  margins = "both",
  marginal.type = "density", # type of marginal distribution to be displayed
  xfill = "#3378ae", # color fill for x-axis marginal distribution
  yfill = "#3378ae", # color fill for y-axis marginal distribution
  xalpha = 0.4, # transparency for x-axis marginal distribution
  yalpha = 0.4, # transparency for y-axis marginal distribution
  centrality.para = "median", # central tendency lines to be displayed
  messages = FALSE, # turn off messages and notes
  results.subtitle = TRUE
)+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 16))
p
ggsave("age_IBD_probability.pdf",p,width = 8,height = 8)

# AVCD, Figure 4c
glmModel <- glm(probability ~ age,data = AVCD,family = "gaussian")
summary(glmModel)
p <- ggscatterstats(
  data = AVCD,
  x = age,
  y = probability,
  type = "nonparametric", # type of test that needs to be run
  conf.level = 0.95, # confidence level
  xlab = "Age (year)", # label for x axis
  ylab = "AVCD probability", # label for y axis
  ggtheme = theme_bw(), # choosing a different theme
  point.args = list(size = 3, alpha = 0.4, stroke = 0,color="#7ea3b5"),
  smooth.line.args = list(size = 1.5, color = "#7ea3b5", method = "lm", formula = y ~ x),
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  margins = "both",
  marginal.type = "density", # type of marginal distribution to be displayed
  xfill = "#3378ae", # color fill for x-axis marginal distribution
  yfill = "#3378ae", # color fill for y-axis marginal distribution
  xalpha = 0.4, # transparency for x-axis marginal distribution
  yalpha = 0.4, # transparency for y-axis marginal distribution
  centrality.para = "median", # central tendency lines to be displayed
  messages = FALSE, # turn off messages and notes
  results.subtitle = TRUE
)+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 16))
p
ggsave("age_AVCD_probability.pdf",p,width = 8,height = 8)

# T2D, not shown
glmModel <- glm(probability ~ age,data = T2D,family = "gaussian")
summary(glmModel)
p <- ggscatterstats(
  data = T2D,
  x = age,
  y = probability,
  type = "nonparametric", # type of test that needs to be run
  conf.level = 0.95, # confidence level
  xlab = "Age (year)", # label for x axis
  ylab = "T2D probability", # label for y axis
  ggtheme = theme_bw(), # choosing a different theme
  point.args = list(size = 3, alpha = 0.4, stroke = 0,color="#7ea3b5"),
  smooth.line.args = list(size = 1.5, color = "#7ea3b5", method = "lm", formula = y ~ x),
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  margins = "both",
  marginal.type = "density", # type of marginal distribution to be displayed
  xfill = "#3378ae", # color fill for x-axis marginal distribution
  yfill = "#3378ae", # color fill for y-axis marginal distribution
  xalpha = 0.4, # transparency for x-axis marginal distribution
  yalpha = 0.4, # transparency for y-axis marginal distribution
  centrality.para = "median", # central tendency lines to be displayed
  messages = FALSE, # turn off messages and notes
  results.subtitle = TRUE
)+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 16))
p
ggsave("age_T2D_probability.pdf",p,width = 8,height = 8)

# IGT,not shown
glmModel <- glm(probability ~ age,data = IGT,family = "gaussian")
summary(glmModel)
p <- ggscatterstats(
  data = IGT,
  x = age,
  y = probability,
  type = "nonparametric", # type of test that needs to be run
  conf.level = 0.95, # confidence level
  xlab = "Age (year)", # label for x axis
  ylab = "IGT probability", # label for y axis
  ggtheme = theme_bw(), # choosing a different theme
  point.args = list(size = 3, alpha = 0.4, stroke = 0,color="#7ea3b5"),
  smooth.line.args = list(size = 1.5, color = "#7ea3b5", method = "lm", formula = y ~ x),
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  margins = "both",
  marginal.type = "density", # type of marginal distribution to be displayed
  xfill = "#3378ae", # color fill for x-axis marginal distribution
  yfill = "#3378ae", # color fill for y-axis marginal distribution
  xalpha = 0.4, # transparency for x-axis marginal distribution
  yalpha = 0.4, # transparency for y-axis marginal distribution
  centrality.para = "median", # central tendency lines to be displayed
  messages = FALSE, # turn off messages and notes
  results.subtitle = TRUE
)+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 16))
p
ggsave("age_IGT_probability.pdf",p,width = 8,height = 8)


# age - top 15 global feature importance ----------------------------------
library(psych)
# IBD, Figure 4e rho and FDR
IBD_species <- colnames(IBD)[9:23]
correlation <- corr.test(IBD$age,IBD[,IBD_species],method = "spearman",adjust="BH")
corr_r <- correlation$r
corr_p <- correlation$p.adj
corr <- t(rbind(corr_r,corr_p))
colnames(corr) <- c("r","p.adj")
write.csv(corr,"age_IBD_species.csv",quote = F,row.names = T)

# AVCD (non-significant)
AVCD_species <- colnames(AVCD)[9:23]
correlation <- corr.test(AVCD$age,AVCD[,AVCD_species],method = "spearman",adjust="BH")
corr_r <- correlation$r
corr_p <- correlation$p.adj
corr <- t(rbind(corr_r,corr_p))
colnames(corr) <- c("r","p.adj")
write.csv(corr,"age_AVCD_species.csv",quote = F,row.names = T)

# T2D (non-significant)
T2D_species <- colnames(T2D)[9:23]
correlation <- corr.test(T2D$age,T2D[,T2D_species],method = "spearman",adjust="BH")
corr_r <- correlation$r
corr_p <- correlation$p.adj
corr <- t(rbind(corr_r,corr_p))
colnames(corr) <- c("r","p.adj")
write.csv(corr,"age_T2D_species.csv",quote = F,row.names = T)

# IGT (non-significant)
IGT_species <- colnames(IGT)[9:23]
correlation <- corr.test(IGT$age,IGT[,IGT_species],method = "spearman",adjust="BH")
corr_r <- correlation$r
corr_p <- correlation$p.adj
corr <- t(rbind(corr_r,corr_p))
colnames(corr) <- c("r","p.adj")
write.csv(corr,"age_IGT_species.csv",quote = F,row.names = T)

# age - top 15 global feature importance of IBD, Figure 4e----------------------------------
# IBD
library(reshape2)
IBD_new <- melt(IBD,measure.vars = IBD_species,variable.name = "species",value.name = "shap_value")
IBD_new$species <- factor(IBD_new$species,levels = unique(IBD_new$species))
p <- ggplot(IBD_new, aes(age, shap_value, group=species))+
  geom_point(alpha=0.4, size=1,shape=20,color="#3372ab")+
  theme_bw()+
  theme(axis.title = element_text(color = "black",size = 10),
        axis.text =  element_text(color = "black",size = 10))+
  geom_smooth(color='#B9BEBA', alpha=0.5)+
  facet_wrap(.~species,scale = "free_y",ncol = 5)+
  xlab('Age')+
  ylab('SHAP value')
ggsave("age_IBD_species.pdf",height = 8,width = 12)


# iterative RF IBD country cohort across age groups, Figure 5a,b------------------------------------------------------------

rfPermutedMultipleNull_IBD_child_adolescent <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                        intersect(IBDCountryCohort,child_adolescent),intersect(IBDCountryCohortTest,c(young,middle,old)),IBDIndividuals,SelectControls,200,10,
                                                                                        round(0.6*length(intersect(IBDIndividuals,child_adolescent)),digits = 0),round(0.6*length(intersect(IBDIndividuals,child_adolescent)),digits = 0))
# trainsize和testsize保持一致，且为训练集或验证集中样本量最小集合中关注结局（此处为Disease）样本量的0.6

rfPermutedMultipleNull_IBD_young <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                            intersect(IBDCountryCohort,young),intersect(IBDCountryCohortTest,c(child_adolescent,middle,old)),IBDIndividuals,SelectControls,200,10,
                                                                                            round(0.6*length(intersect(IBDIndividuals,young)),digits = 0),round(0.6*length(intersect(IBDIndividuals,young)),digits = 0))

rfPermutedMultipleNull_IBD_middle <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                 intersect(IBDCountryCohort,middle),intersect(IBDCountryCohortTest,c(child_adolescent,young,old)),IBDIndividuals,SelectControls,200,10,
                                                                                 round(0.6*length(intersect(IBDIndividuals,middle)),digits = 0),round(0.6*length(intersect(IBDIndividuals,middle)),digits = 0))

rfPermutedMultipleNull_IBD_old <- iterative_rf_multiple_with_null_distribution(clr(dat_disease[,bac_species_name]),
                                                                                 intersect(IBDCountryCohort,old),intersect(IBDCountryCohortTest,c(child_adolescent,middle,young)),IBDIndividuals,SelectControls,200,10,
                                                                                 round(0.6*length(intersect(IBDIndividuals,old)),digits = 0),round(0.6*length(intersect(IBDIndividuals,old)),digits = 0))

rfTestAUCs_IBD_age <- as.data.frame(matrix(NA,1600,4))
colnames(rfTestAUCs_IBD_age) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_IBD_age$AUC <- c(apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCDiffAge,1,median),
                                         apply(rfPermutedMultipleNull_IBD_young$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IBD_young$AUCDiffAge,1,median),
                                         apply(rfPermutedMultipleNull_IBD_middle$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IBD_middle$AUCDiffAge,1,median),
                                         apply(rfPermutedMultipleNull_IBD_old$AUCSameAge,1,median),apply(rfPermutedMultipleNull_IBD_old$AUCDiffAge,1,median))
rfTestAUCs_IBD_age$ClassificationGroup <- factor(c(rep("Same Age",200),rep("Different Age",200),rep("Same Age",200),rep("Different Age",200),rep("Same Age",200),rep("Different Age",200),rep("Same Age",200),rep("Different Age",200)),levels=c("Same Age","Different Age"))
rfTestAUCs_IBD_age$Group <- factor(c(rep("Child&Adolescent",400),rep("Young",400),rep("Middle-aged",400),rep("Elderly",400)),levels = c("Child&Adolescent","Young","Middle-aged","Elderly"))
rfTestAUCs_IBD_age$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

write.csv(rfTestAUCs_IBD_age,"rfTestAUCs_IBD_age.csv",row.names = F,quote = F)

group_by(rfTestAUCs_IBD_age[,c("AUC","ClassificationGroup","Group")],Group,ClassificationGroup) %>% summarise_each(funs(mean,sd))

# Figure 5a
p <- ggplot(rfTestAUCs_IBD_age)+
  # 箱线图：
  geom_boxplot(aes(x = ClassificationGroup, y = AUC,color=ClassificationGroup),
               width=0.7,position = position_dodge(0),size=0.5,outlier.size = 0) + 
  # 散点：
  geom_jitter(aes(x = ClassificationGroup, y = AUC,fill=ClassificationGroup),
              width=0.2,shape=21,colour="black",size=1)+
  # 设置填充颜色
  scale_fill_manual(values=c("#00ba38","#619cff"))+
  scale_color_manual(values=c("#00ba38","#619cff"))+
  # 添加显著性检验：
  # stat_compare_means(aes(x = ClassificationGroup, y = AUC),
  #                   label="p.signif", method = "anova")+
  # 分面：
  facet_grid(.~Group, labeller = label_both)+
  # X轴y轴标题：
  labs(y="AUC",x="Age Group")+
  # 设置主题：
  theme_minimal()+
  theme(strip.background = element_rect(fill="#eaeae0", color = "#dcddcf"),  # 分面颜色
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text.y = element_text(size=11,face="plain",color="black"),
        axis.text.x = element_blank(),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),
        # legend.position="none",
        # legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks.x=element_blank())
ggsave("boxplot_line_IBD_age.pdf", height = 6, width = 12)
# geom_jitter添加抖动效果，每次成图时点的水平位置会有随机性差异，但垂直位置保持一致

wilcox.test(AUC~ClassificationGroup,data=subset(rfTestAUCs_IBD_age,Group=="Child&Adolescent",paired=TRUE))
wilcox.test(AUC~ClassificationGroup,data=subset(rfTestAUCs_IBD_age,Group=="Young",paired=TRUE))
wilcox.test(AUC~ClassificationGroup,data=subset(rfTestAUCs_IBD_age,Group=="Middle-aged",paired=TRUE))
wilcox.test(AUC~ClassificationGroup,data=subset(rfTestAUCs_IBD_age,Group=="Elderly",paired=TRUE))


# permutation test
rfPermutationTest_IBD_age <- as.data.frame(matrix(NA,1600,4))
colnames(rfPermutationTest_IBD_age) <- c("AUCDiff","TestType","Training","Iter")
rfPermutationTest_IBD_age$AUCDiff <- c(apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCSameAge,1,median)-apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCDiffAge,1,median),
                                       apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCNull2,1,median)-apply(rfPermutedMultipleNull_IBD_child_adolescent$AUCNull1,1,median),
                                       apply(rfPermutedMultipleNull_IBD_young$AUCSameAge,1,median)-apply(rfPermutedMultipleNull_IBD_young$AUCDiffAge,1,median),
                                       apply(rfPermutedMultipleNull_IBD_young$AUCNull2,1,median)-apply(rfPermutedMultipleNull_IBD_young$AUCNull1,1,median),
                                       apply(rfPermutedMultipleNull_IBD_middle$AUCSameAge,1,median)-apply(rfPermutedMultipleNull_IBD_middle$AUCDiffAge,1,median),
                                       apply(rfPermutedMultipleNull_IBD_middle$AUCNull2,1,median)-apply(rfPermutedMultipleNull_IBD_middle$AUCNull1,1,median),
                                       apply(rfPermutedMultipleNull_IBD_old$AUCSameAge,1,median)-apply(rfPermutedMultipleNull_IBD_old$AUCDiffAge,1,median),
                                       apply(rfPermutedMultipleNull_IBD_old$AUCNull2,1,median)-apply(rfPermutedMultipleNull_IBD_old$AUCNull1,1,median))
rfPermutationTest_IBD_age$TestType <- factor(c(rep("Actual Distribution",200),rep("Random Distribution",200),rep("Actual Distribution",200),rep("Random Distribution",200),rep("Actual Distribution",200),rep("Random Distribution",200),rep("Actual Distribution",200),rep("Random Distribution",200)),levels=c("Actual Distribution","Random Distribution"))
rfPermutationTest_IBD_age$Training <- factor(c(rep("Child&Adolescent",400),rep("Young",400),rep("Middle-aged",400),rep("Elderly",400)),levels = c("Child&Adolescent","Young","Middle-aged","Elderly"))
rfPermutationTest_IBD_age$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

group_by(rfPermutationTest_IBD_age[,c("AUCDiff","TestType","Training")],Training,TestType) %>% summarise_each(funs(mean,sd))

# Figure 5b
p <- ggplot(rfPermutationTest_IBD_age)+
  # 箱线图：
  geom_boxplot(aes(x = TestType, y = AUCDiff,color=TestType),
               width=0.7,position = position_dodge(0),size=0.5,outlier.size = 0) + 
  # 散点：
  geom_jitter(aes(x = TestType, y = AUCDiff,fill=TestType),
              width=0.2,shape=21,colour="black",size=1)+
  # 设置填充颜色
  scale_fill_manual(values=c("#e4933d","#8684b9"))+
  scale_color_manual(values=c("#e4933d","#8684b9"))+
  # 添加显著性检验：
  # stat_compare_means(aes(x = TestType, y = AUCDiff),
  #                   label="p.signif", method = "anova")+
  # 分面：
  facet_grid(.~Training, labeller = label_both)+
  # X轴y轴标题：
  labs(y="AUCDiff",x="Distribution")+
  # 设置主题：
  theme_minimal()+
  theme(strip.background = element_rect(fill="#eaeae0", color = "#dcddcf"),  # 分面颜色
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text.y = element_text(size=11,face="plain",color="black"),
        axis.text.x = element_blank(),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),
        # legend.position="none",
        # legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks.x=element_blank())
ggsave("boxplot_line_IBD_age_distribution.pdf", height = 6, width = 12)
# geom_jitter添加抖动效果，每次成图时点的水平位置会有随机性差异，但垂直位置保持一致

wilcox.test(AUCDiff~TestType,data = subset(rfPermutationTest_IBD_age,Training=="Child&Adolescent",paired=TRUE))
wilcox.test(AUCDiff~TestType,data = subset(rfPermutationTest_IBD_age,Training=="Young",paired=TRUE))
wilcox.test(AUCDiff~TestType,data = subset(rfPermutationTest_IBD_age,Training=="Middle-aged",paired=TRUE))
wilcox.test(AUCDiff~TestType,data = subset(rfPermutationTest_IBD_age,Training=="Elderly",paired=TRUE))

# threshold determination

rfThreshold <- function(data,name){
  threshold <- seq(0,0.99,0.01)
  rf_threshold <- as.data.frame(matrix(NA,ncol = 3))
  colnames(rf_threshold) <-c("threshold", paste0("mean_featureScore",name),paste0("bac_number",name))
  for (i in 1:length(threshold)){
    rfBacteria=names(which(rank_scale(colMeans(data$featureProfile))>=threshold[i]))
    rf_threshold[i,1] = threshold[i]
    rf_threshold[i,2] = mean(colMeans(data$featureProfile[,rfBacteria]))
    rf_threshold[i,3] = length(rfBacteria)
  }
  return(rf_threshold)
}

rfThreshold_IBD_child_adolescent <- rfThreshold(rfPermutedMultipleNull_IBD_child_adolescent,"child_adolescent")
rfThreshold_IBD_young <- rfThreshold(rfPermutedMultipleNull_IBD_young,"young")
rfThreshold_IBD_middle <- rfThreshold(rfPermutedMultipleNull_IBD_middle,"middle")
rfThreshold_IBD_old <- rfThreshold(rfPermutedMultipleNull_IBD_old,"old")

rfThreshold_IBD_age <- Reduce(function(x,y) merge(x,y,by="threshold",all=T),list(rfThreshold_IBD_child_adolescent,rfThreshold_IBD_young,rfThreshold_IBD_middle,rfThreshold_IBD_old))
rfThreshold_IBD_age_melt <- melt(rfThreshold_IBD_age[,c(1,2,4,6,8)],id.vars = "threshold",value.name = "mean",variable.name = "IBD_age")

# Figure S7a------------------------------------------------------------------
p <- ggplot(data = rfThreshold_IBD_age_melt, aes(x = threshold, y = mean, color = IBD_age)) +
  geom_line(data = rfThreshold_IBD_age_melt, aes(x = threshold, y = mean),linewidth=1.5) +
  scale_color_manual(values = c("#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="Percentile Rank of Mean Decrease in Gini ",y="The Mean of Mean Decrease in Gini of Remaining Species")+
  scale_x_continuous(breaks = seq(0,1,0.1),
                     label = c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))+
  geom_vline(xintercept = 0.85,linetype="dashed",color="black")+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12))
ggsave("rfThreshold_IBD_age_disease.pdf",width = 7,height = 6)

# species marker

rfSpeciesMarker_IBD_child_adolescent<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IBD_child_adolescent$featureProfile))>=0.85))
rfSpeciesMarker_IBD_young<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IBD_young$featureProfile))>=0.85))
rfSpeciesMarker_IBD_middle<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IBD_middle$featureProfile))>=0.85))
rfSpeciesMarker_IBD_old<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_IBD_old$featureProfile))>=0.85));rfSpeciesMarker_IBD_old[24] <- NA

rfSpeciesMarker_IBD_age <- data.frame(rfSpeciesMarker_IBD_child_adolescent,rfSpeciesMarker_IBD_young,rfSpeciesMarker_IBD_middle,rfSpeciesMarker_IBD_old)
write.csv(rfSpeciesMarker_IBD_age,"rfSpeciesMarker_IBD_age.csv",quote = F,row.names = F)

# species marker comparisons

rfSpeciesMarker_IBD_age_unique <- unique(c(rfSpeciesMarker_IBD_child_adolescent,rfSpeciesMarker_IBD_young,rfSpeciesMarker_IBD_middle,rfSpeciesMarker_IBD_old))
rfSpeciesMarker_IBD_age_unique <- rfSpeciesMarker_IBD_age_unique[!is.na(rfSpeciesMarker_IBD_age_unique)]
rfFeatureProfile_IBD_age <- rbind(rfPermutedMultipleNull_IBD_child_adolescent$featureProfile,rfPermutedMultipleNull_IBD_young$featureProfile,rfPermutedMultipleNull_IBD_middle$featureProfile,rfPermutedMultipleNull_IBD_old$featureProfile)
rfMeanMarkerScores_IBD_age <- cbind(colMeans(rfFeatureProfile_IBD_age[1:200,]),colMeans(rfFeatureProfile_IBD_age[201:400,]),colMeans(rfFeatureProfile_IBD_age[401:600,]),colMeans(rfFeatureProfile_IBD_age[601:800,]))
colnames(rfMeanMarkerScores_IBD_age) <- c("MeanFeatureProfile_IBD_child_adolescent","MeanFeatureProfile_IBD_young","MeanFeatureProfile_IBD_middle","MeanFeatureProfile_IBD_child_old")
rfComparison_IBD_age <- batch_dunns(rfFeatureProfile_IBD_age[,rfSpeciesMarker_IBD_age_unique],as.factor(c(rep("child_adolescent",200),rep("young",200),rep("middle",200),rep("old",200))))
rfComparisonSummary_IBD_age <- cbind(rfMeanMarkerScores_IBD_age[rfSpeciesMarker_IBD_age_unique,],rfComparison_IBD_age)
rfComparisonSummary_IBD_age <- as.data.frame(rfComparisonSummary_IBD_age)
rfComparisonSummary_IBD_age$IBD_child_adolescent <- ifelse(rownames(rfComparisonSummary_IBD_age) %in% rfSpeciesMarker_IBD_child_adolescent,1,0)
rfComparisonSummary_IBD_age$IBD_young <- ifelse(rownames(rfComparisonSummary_IBD_age) %in% rfSpeciesMarker_IBD_young,1,0)
rfComparisonSummary_IBD_age$IBD_middle <- ifelse(rownames(rfComparisonSummary_IBD_age) %in% rfSpeciesMarker_IBD_middle,1,0)
rfComparisonSummary_IBD_age$IBD_old <- ifelse(rownames(rfComparisonSummary_IBD_age) %in% rfSpeciesMarker_IBD_old,1,0)
write.csv(rfComparisonSummary_IBD_age,"rfComparisonSummary_IBD_age.csv",quote = F,row.names = T)

# heatmap of IBD age-specific species marker, Figure S7b-------------------------------------------
IBD_age_specific_markers <- read.csv("IBD_age_specific_markers.csv",header = T,row.names = 1)
library(pheatmap)
pdf("IBD_age_species_marker.pdf",height = 15,width = 5)
pheatmap(IBD_age_specific_markers,scale = "none",cluster_rows = F,cluster_cols = F,
         fontsize = 8,border_color = "black",
         cellwidth = 12,cellheight = 12,
         color = colorRampPalette(c("#c2d7e7","#9abcd7","#679ac2","#3378ae","#00579a"))(500))
dev.off()

# Bar plots, Figure 5c, d ---------------------------------------------------------------
IBD_proportion_1 <- read.xlsx("IBD_proportion.xlsx",sheet = 1)
IBD_proportion_2 <- read.xlsx("IBD_proportion.xlsx",sheet = 2)

IBD_proportion_1$group <- factor(IBD_proportion_1$group,levels = unique(IBD_proportion_1$group))
IBD_proportion_2$group <- factor(IBD_proportion_2$group,levels = unique(IBD_proportion_2$group))

IBD_proportion_1$category <- factor(IBD_proportion_1$category,levels = rev(unique(IBD_proportion_1$category)))
IBD_proportion_2$category <- factor(IBD_proportion_2$category,levels = unique(IBD_proportion_2$category))

p <- ggplot(IBD_proportion_1,aes(group, value, fill = category))+
  geom_bar(aes(group, value, fill = category),
           stat = "identity", position = "fill", color = "grey")+
  scale_fill_manual(values = c("#979797","#feda77"))+
  ylab("Proportion (%)")+
  xlab("")+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(color = "black",size=10),
        axis.text = element_text(color = "black",size=10))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","25%","50%","75%","100%"))+
  coord_flip()
ggsave("IBD_species_shared.pdf",width = 5,height = 3)

p <- ggplot(IBD_proportion_2,aes(group, value, fill = category))+
  geom_bar(aes(group, value, fill = category),
           stat = "identity", position = "fill", color = "grey")+
  scale_fill_manual(values = c("#a0c6e9","#eeeb65"))+
  ylab("Proportion (%)")+
  xlab("")+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(color = "black",size=10),
        axis.text = element_text(color = "black",size=10))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","25%","50%","75%","100%"))+
  coord_flip()
ggsave("age_group_species_shared.pdf",width = 5,height = 3)

# IBD Bray-Curtis Comparisons between disease status and age groups, Figure S8------------------------------------------------

dat_disease_2 <- dat_disease[IBDCountryCohort,]
bray_dis <- as.matrix(vegdist(dat_disease_2[,name_all[bac_species]],method = "bray",diag = TRUE,upper = TRUE))

pcoa <- cmdscale(bray_dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) <- c("x", "y", "z") 
eig <- pcoa$eig
points_2 <- cbind(points, dat_disease_2$age_group,dat_disease_2$disease)
colnames(points_2)[4] <- "age_group"
colnames(points_2)[5] <- "disease"

library(reshape2)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- detectCores()
registerDoParallel(cl)

# age groups
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ age_group+DNA_extraction_kit, data = dat_disease_2,permutations = 999)}
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.085 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=age_group)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=age_group), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("IBD_bray_age_group_after_selected.pdf", p, width = 8, height = 6)

# disease
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ disease+DNA_extraction_kit, data = dat_disease_2,permutations = 999)}
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.085 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=disease)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=age_group), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("IBD_bray_disease_after_selected.pdf", p, width = 7.5, height = 6)

# spearman distance change of IBD across age groups, Figure 5e, Figure S7c-e------------------------------------------------------

IBD_child_adolescent_featureProfile <- rfPermutedMultipleNull_IBD_child_adolescent$featureProfile
IBD_young_featureProfile <- rfPermutedMultipleNull_IBD_young$featureProfile
IBD_middle_featureProfile <- rfPermutedMultipleNull_IBD_middle$featureProfile
IBD_old_featureProfile <- rfPermutedMultipleNull_IBD_old$featureProfile

# child_adolescent, Figure 5e---------------------------------------------------------------
matrices <- list(IBD_child_adolescent_featureProfile[,rfSpeciesMarker_IBD_child_adolescent],
                 IBD_young_featureProfile[,rfSpeciesMarker_IBD_child_adolescent],
                 IBD_middle_featureProfile[,rfSpeciesMarker_IBD_child_adolescent],
                 IBD_old_featureProfile[,rfSpeciesMarker_IBD_child_adolescent])  # Add all the matrices you want to compare to this list
matrixA <- IBD_child_adolescent_featureProfile[,rfSpeciesMarker_IBD_child_adolescent]
distances_Child <- data.frame(matrix(NA,ncol = 4,nrow = 576))
for (i in 1:length(matrices)) {
  distances_Child[i] <- spearman_distance(matrixA, matrices[[i]])
}
colnames(distances_Child) <- c("Child & Adolescent","Young","Middle-aged","Elderly")
colMeans(distances_Child)

distances_Child <- melt(distances_Child,measure.vars=colnames(distances_Child),value.name="spearman_dis",variable.name="age_group")
# distances_Child <- subset(distances_Child,spearman_dis!=0)
p <- ggviolin(distances_Child,x="age_group",y="spearman_dis",add=c("boxplot","dotplot"),add.params = list(fill="white",binwidth=0.010,color="age_group",alpha=0.2,width=0.1,size=0.6),
              fill = "white",color = "age_group",width = 1,size = 1.2,position = "dodge",trim = T)+
  geom_pwc(aes(group=age_group),method = "wilcox.test",p.adjust.method="BH",label = "p.adj.signif")+
  # stat_compare_means(comparisons = boxplot_comparison,method = "wilcox.test",p.adjust.methods="BH")+
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = "Control") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(color="black",size=17),
        axis.text = element_text(color="black",size=15))+
  scale_color_manual(values = c("#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="",y="Spearman distance")
p
ggsave("child_spearman_dis_comparison2.pdf",width = 7.5,height = 5)

# young, Figure S7c---------------------------------------------------------------
matrices <- list(IBD_child_adolescent_featureProfile[,rfSpeciesMarker_IBD_young],
                 IBD_young_featureProfile[,rfSpeciesMarker_IBD_young],
                 IBD_middle_featureProfile[,rfSpeciesMarker_IBD_young],
                 IBD_old_featureProfile[,rfSpeciesMarker_IBD_young])  # Add all the matrices you want to compare to this list
matrixA <- IBD_young_featureProfile[,rfSpeciesMarker_IBD_young]
distances_Young <- data.frame(matrix(NA,ncol = 4,nrow = 576))
for (i in 1:length(matrices)) {
  distances_Young[i] <- spearman_distance(matrixA, matrices[[i]])
}
colnames(distances_Young) <- c("Child & Adolescent","Young","Middle-aged","Elderly")
colMeans(distances_Young)

distances_Young <- melt(distances_Young,measure.vars=colnames(distances_Young),value.name="spearman_dis",variable.name="age_group")
# distances_Young <- subset(distances_Young,spearman_dis!=0)
p <- ggviolin(distances_Young,x="age_group",y="spearman_dis",add=c("boxplot","dotplot"),add.params = list(fill="white",binwidth=0.010,color="age_group",alpha=0.2,width=0.1,size=0.6),
              fill = "white",color = "age_group",width = 1,size = 1.2,position = "dodge",trim = T)+
  geom_pwc(aes(group=age_group),method = "wilcox.test",p.adjust.method="BH",label = "p.adj.signif")+
  # stat_compare_means(comparisons = boxplot_comparison,method = "wilcox.test",p.adjust.methods="BH")+
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = "Control") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(color="black",size=17),
        axis.text = element_text(color="black",size=15))+
  scale_color_manual(values = c("#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="",y="Spearman distance")
p
ggsave("young_spearman_dis_comparison2.pdf",width = 7.5,height = 5)


# middle, Figure S7d---------------------------------------------------------------
matrices <- list(IBD_child_adolescent_featureProfile[,rfSpeciesMarker_IBD_middle],
                 IBD_young_featureProfile[,rfSpeciesMarker_IBD_middle],
                 IBD_middle_featureProfile[,rfSpeciesMarker_IBD_middle],
                 IBD_old_featureProfile[,rfSpeciesMarker_IBD_middle])  # Add all the matrices you want to compare to this list
matrixA <- IBD_middle_featureProfile[,rfSpeciesMarker_IBD_middle]
distances_Middle <- data.frame(matrix(NA,ncol = 4,nrow = 576))
for (i in 1:length(matrices)) {
  distances_Middle[i] <- spearman_distance(matrixA, matrices[[i]])
}
colnames(distances_Middle) <- c("Child & Adolescent","Young","Middle-aged","Elderly")
colMeans(distances_Middle)

distances_Middle <- melt(distances_Middle,measure.vars=colnames(distances_Middle),value.name="spearman_dis",variable.name="age_group")
# distances_Middle <- subset(distances_Middle,spearman_dis!=0)
p <- ggviolin(distances_Middle,x="age_group",y="spearman_dis",add=c("boxplot","dotplot"),add.params = list(fill="white",binwidth=0.010,color="age_group",alpha=0.2,width=0.1,size=0.6),
              fill = "white",color = "age_group",width = 1,size = 1.2,position = "dodge",trim = T)+
  geom_pwc(aes(group=age_group),method = "wilcox.test",p.adjust.method="BH",label = "p.adj.signif")+
  # stat_compare_means(comparisons = boxplot_comparison,method = "wilcox.test",p.adjust.methods="BH")+
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = "Control") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(color="black",size=17),
        axis.text = element_text(color="black",size=15))+
  scale_color_manual(values = c("#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="",y="Spearman distance")
p
ggsave("middle_spearman_dis_comparison2.pdf",width = 7.5,height = 5)

# old, Figure S7e---------------------------------------------------------------
rfSpeciesMarker_IBD_old <- rfSpeciesMarker_IBD_old[1:23]
matrices <- list(IBD_child_adolescent_featureProfile[,rfSpeciesMarker_IBD_old],
                 IBD_young_featureProfile[,rfSpeciesMarker_IBD_old],
                 IBD_middle_featureProfile[,rfSpeciesMarker_IBD_old],
                 IBD_old_featureProfile[,rfSpeciesMarker_IBD_old])  # Add all the matrices you want to compare to this list
matrixA <- IBD_old_featureProfile[,rfSpeciesMarker_IBD_old]
distances_old <- data.frame(matrix(NA,ncol = 4,nrow = 529))
for (i in 1:length(matrices)) {
  distances_old[i] <- spearman_distance(matrixA, matrices[[i]])
}
colnames(distances_old) <- c("Child & Adolescent","Young","Middle-aged","Elderly")
colMeans(distances_old)

distances_old <- melt(distances_old,measure.vars=colnames(distances_old),value.name="spearman_dis",variable.name="age_group")
# distances_old <- subset(distances_old,spearman_dis!=0)
p <- ggviolin(distances_old,x="age_group",y="spearman_dis",add=c("boxplot","dotplot"),add.params = list(fill="white",binwidth=0.010,color="age_group",alpha=0.2,width=0.1,size=0.6),
              fill = "white",color = "age_group",width = 1,size = 1.2,position = "dodge",trim = T)+
  geom_pwc(aes(group=age_group),method = "wilcox.test",p.adjust.method="BH",label = "p.adj.signif")+
  # stat_compare_means(comparisons = boxplot_comparison,method = "wilcox.test",p.adjust.methods="BH")+
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = "Control") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(color="black",size=17),
        axis.text = element_text(color="black",size=15))+
  scale_color_manual(values = c("#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="",y="Spearman distance")
p
ggsave("old_spearman_dis_comparison2.pdf",width = 7.5,height = 5)


# IBD species comparison, Figure 5f --------------------------------------------------

effect_size_calculator=function(x,y)
{
  library(effsize);
  species_list <- intersect(colnames(x),colnames(y));
  effect_size <- as.data.frame(matrix(NA,length(species_list),5))
  rownames(effect_size) <- species_list;
  for(i in 1:length(species_list))
  {
    print(species_list[i]);
    gg <- effsize::cohen.d(x[,species_list[i]],y[,species_list[i]])
    #print(gg);
    gt <- wilcox.test(x[,species_list[i]],y[,species_list[i]])
    effect_size[i,1] <- as.numeric(gg$estimate[1]);
    effect_size[i,2] <- as.factor(gg$magnitude[1]);
    effect_size[i,3] <- as.numeric(gt$p.value[1]);
    effect_size[i,4] <- mean(x[,species_list[i]]);
    effect_size[i,5] <- mean(y[,species_list[i]]);
  }
  colnames(effect_size) <- c("Estimate","Level","NominalP","MeanDisease","MedanControl");
  effect_size$p.adj <- p.adjust(effect_size$NominalP,method = "BH")
  effect_size <- apply(effect_size,2,function(x)(ifelse(is.nan(x),0,x)))
  #effect_size <- effect_size[names(which(!is.nan(rowSums(effect_size)))),];
  return(effect_size);
}

rfSpeciesMarker_IBD_age_intersect <- intersect(rfSpeciesMarker_IBD,rfSpeciesMarker_IBD_age_unique)
dat_disease_2 <- dat_disease[IBDCountryCohort,]

EffSize_IBD_child_adolescent <- effect_size_calculator(dat_disease_2[intersect(intersect(IBDIndividuals,IBDCountryCohort),c(child_adolescent)),rfSpeciesMarker_IBD_age_intersect],
                                                       dat_disease_2[intersect(intersect(SelectControls,IBDCountryCohort),c(child_adolescent)),rfSpeciesMarker_IBD_age_intersect])
EffSize_IBD_young <- effect_size_calculator(dat_disease_2[intersect(intersect(IBDIndividuals,IBDCountryCohort),c(young)),rfSpeciesMarker_IBD_age_intersect],
                                            dat_disease_2[intersect(intersect(SelectControls,IBDCountryCohort),c(young)),rfSpeciesMarker_IBD_age_intersect])
EffSize_IBD_middle <- effect_size_calculator(dat_disease_2[intersect(intersect(IBDIndividuals,IBDCountryCohort),c(middle)),rfSpeciesMarker_IBD_age_intersect],
                                             dat_disease_2[intersect(intersect(SelectControls,IBDCountryCohort),c(middle)),rfSpeciesMarker_IBD_age_intersect])
EffSize_IBD_old <- effect_size_calculator(dat_disease_2[intersect(intersect(IBDIndividuals,IBDCountryCohort),c(old)),rfSpeciesMarker_IBD_age_intersect],
                                          dat_disease_2[intersect(intersect(SelectControls,IBDCountryCohort),c(old)),rfSpeciesMarker_IBD_age_intersect])
EffSize_IBD_all <- data.frame(cbind(EffSize_IBD_child_adolescent[,1],EffSize_IBD_young[,1],EffSize_IBD_middle[,1],EffSize_IBD_old[,1]))
EffSize_IBD_all <- abs(EffSize_IBD_all)
colnames(EffSize_IBD_all) <- c("Child & Adolescent","Young","Middle-aged","Elderly")
EffSize_IBD_all_p <- cbind(EffSize_IBD_child_adolescent[,6],EffSize_IBD_young[,6],EffSize_IBD_middle[,6],EffSize_IBD_old[,6])
EffSize_IBD_all_p <- matrix(ifelse(
  EffSize_IBD_all_p < 0.05 & EffSize_IBD_all_p >= 0.01 ,"*",ifelse(
    EffSize_IBD_all_p < 0.01 & EffSize_IBD_all_p >= 0.001,"**",ifelse(
      EffSize_IBD_all_p < 0.001,"***",""
    ))),nrow(EffSize_IBD_all_p))

library(pheatmap)
pdf("IBD_speciesMarker_Cohend_ageGroup.pdf")
pheatmap(EffSize_IBD_all,scale = "none",cluster_rows = T,cluster_cols = F,legend_breaks = c(0,0.2,0.5,0.8,2.1),legend_labels = c("0","0.2","0.5","0.8","2.1"),
         fontsize = 8,border_color = "black",display_numbers=EffSize_IBD_all_p,
         cellwidth = 12,cellheight = 12,number_format = 8,number_color = "black",
         color = c(colorRampPalette(c("#c2d7e7"))(190),colorRampPalette(c("#679ac2"))(300),colorRampPalette(c("#3378ae"))(300),colorRampPalette(c("#00579a"))(1287)))
dev.off()













