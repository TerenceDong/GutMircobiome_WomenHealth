setwd("~/gut_health/health_code_recheck")
source("code_library.R")
# batch effect before adjusting------------------------------------------------------------
library(dplyr)
library(openxlsx)

d1_health <- read.xlsx("sampleMetadata.xlsx")
d2_health <- read.csv("RelativeAbundance_health.csv",row.names = 1)

d1_health$DNA_extraction_kit[is.na(d1_health$DNA_extraction_kit)==TRUE] <- "na" 
d1_health$sample_id <- gsub("\\.","-",d1_health$sample_id)

d2_reAbund_feature <- as.data.frame(matrix(NA,ncol = 3))
colnames(d2_reAbund_feature) <- c("species","prevalence","proportion")
Species <- colnames(d2_health)
for (i in 1:ncol(d2_health)){
  d2_reAbund_feature[i,1] = Species[i]
  d2_reAbund_feature[i,2] = 1-sum(d2_health[,Species[i]]==0)/nrow(d2_health)
  d2_reAbund_feature[i,3] = sum(d2_health[,Species[i]])/sum(d2_health)
}
d2_reAbund_select <- d2_reAbund_feature$species[d2_reAbund_feature$prevalence >= 0.1 & d2_reAbund_feature$proportion >= 0.00005]

d2_health_select <- d2_health[,d2_reAbund_select]
d2_health_select$sample_id <- row.names(d2_health_select)
d2_health_select$sample_id <- gsub("\\.","-",d2_health_select$sample_id)
d2_health_select <- merge(d2_health_select,d1_health,by="sample_id",all=FALSE)


# Figure S1-------------------------------------------------------------------- -----------------------------------------

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(d2_health_select[,d2_reAbund_select],method = "bray",diag = TRUE,upper = TRUE))
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

# sequencing_platform, Figure S1c -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select$sequencing_platform)
colnames(points_2)[4] <- "sequencing_platform"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ sequencing_platform, data = d2_health_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ sequencing_platform, data = d2_health_select,permutations = 999)
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

# DNA_extraction_kit, Figure S1a -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select$DNA_extraction_kit)
colnames(points_2)[4] <- "DNA_extraction_kit"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ DNA_extraction_kit, data = d2_health_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ DNA_extraction_kit, data = d2_health_select,permutations = 999)
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
  scale_color_manual(values = c("#91bde2", "#83c5b0", "#8fb79d","#eadfc3", "#eec16c","#d75032", "#ce9c83","#b9b7a2", "#f66f69","#289ebc", "#e8edf3","#8171b5"))
ggsave("bray_DNA_extraction_kit.pdf", p, width = 8, height = 6)

# median_read_length, Figure S1b -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select$median_read_length)
colnames(points_2)[4] <- "median_read_length"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ median_read_length, data = d2_health_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ median_read_length, data = d2_health_select,permutations = 999)
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

d2_health_select2 <- d2_health_select[,d2_reAbund_select]
rownames(d2_health_select2) <- d2_health_select$sample_id
d2_health_select2 <- t(d2_health_select2) # feature-by-sample
d2_health_select2 <- d2_health_select2/100 # divided by 100 to be compositional data

d2_health_select2_meta <- d2_health_select[,colnames(d1_health)] # sample-by-metadata
rownames(d2_health_select2_meta) <- d2_health_select$sample_id

d2_health_select2_meta$DNA_extraction_kit <- factor(d2_health_select2_meta$DNA_extraction_kit) # assume the variable adjusted to be factor variable
d2_health_select2_adjust1 <- adjust_batch(feature_abd = d2_health_select2,
                                          batch = c("DNA_extraction_kit"),
                                          # covariates = "study_condition",
                                          data = d2_health_select2_meta)$feature_abd_adj

d2_health_select2_meta$sequencing_platform <- factor(d2_health_select2_meta$sequencing_platform) # assume the variable adjusted to be factor variable
d2_health_select2_adjust2 <- adjust_batch(feature_abd = d2_health_select2_adjust1,
                                          batch = c("sequencing_platform"),
                                          # covariates = "study_condition",
                                          data = d2_health_select2_meta)$feature_abd_adj

d2_health_select2_meta$median_read_length <- factor(d2_health_select2_meta$median_read_length) # assume the variable adjusted to be factor variable
d2_health_select2_adjust3 <- adjust_batch(feature_abd = d2_health_select2_adjust2,
                                          batch = c("median_read_length"),
                                          # covariates = "study_condition",
                                          data = d2_health_select2_meta)$feature_abd_adj

# batch effect after adjusting, Figure S1------------------------------------------------------------

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(t(d2_health_select2_adjust3),method = "bray",diag = TRUE,upper = TRUE))

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

# sequencing_platform, Figure S1c -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select2_meta$sequencing_platform)
colnames(points_2)[4] <- "sequencing_platform"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ sequencing_platform, data = d2_health_select2_meta,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ sequencing_platform, data = d2_health_select2_meta,permutations = 999)
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

# DNA_extraction_kit, Figure S1a -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select2_meta$DNA_extraction_kit)
colnames(points_2)[4] <- "DNA_extraction_kit"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ DNA_extraction_kit, data = d2_health_select2_meta,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ DNA_extraction_kit, data = d2_health_select2_meta,permutations = 999)
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
  scale_color_manual(values = c("#91bde2", "#83c5b0", "#8fb79d","#eadfc3", "#eec16c","#d75032", "#ce9c83","#b9b7a2", "#f66f69","#289ebc", "#e8edf3","#8171b5"))
ggsave("bray_DNA_extraction_kit_after.pdf", p, width = 8, height = 6)

# median_read_length, Figure S1b -----------------------------------------------------------
points_2 <- cbind(points, d2_health_select$median_read_length)
colnames(points_2)[4] <- "median_read_length"

# PERMANOVA
set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ median_read_length, data = d2_health_select,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ median_read_length, data = d2_health_select,permutations = 999)
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

d2_health_select2_adjust3_t <- data.frame(t(d2_health_select2_adjust3))
d2_health_select2_adjust3_t$sample_id <- rownames(d2_health_select2_adjust3_t)
d2_health_select2_adjust3_t_merge <- merge(d2_health_select2_adjust3_t,d2_health_select2_meta,by="sample_id",all = T)
write.xlsx(d2_health_select2_adjust3_t_merge,"RelAbund_batchCorrect.xlsx")

# Figure 1B; microbial description of healthy women------------------------------------- ----------------------------------
dat_health <- read.xlsx("RelAbund_batchCorrect.xlsx")
continent <- data.frame(
  country=c("CAN","CHN","DEU","DNK","EST","FIN","FRA","GBR","IND","IRL","ISL","ITA","LUX","MDG","NLD","RUS","SWE","TZA","USA"),
  continent=c("North America","Asia","Europe","Europe","Europe","Europe","Europe","Europe","Asia","Europe","Europe","Europe","Europe","Africa","Europe","Europe","Europe","Africa","North America")
)
dat_health <- merge(dat_health,continent,by="country") # use merge function to simulate vlookup
bac_species <- c(3:172)

short_phylum <- read.xlsx("short_phylum.xlsx",sheet = 1)
phylum <- unique(short_phylum$phylum)
for (i in phylum){
  dat=subset(short_phylum,phylum==i)
  dat_health[,i]=tryCatch(rowSums(dat_health[,dat$short]),error=function(e){"error"})
}
dat_health$Euryarchaeota <- dat_health$Methanobrevibacter.smithii
dat_health$Verrucomicrobia <- dat_health$Akkermansia.muciniphila
bac_phylum <- c(308:313)

dat_health$age_group <- ifelse(dat_health$age <=1,"infant",
                               ifelse(dat_health$age >=2 & dat_health$age <=17,"child_adolescent",
                                      ifelse(dat_health$age >=18 & dat_health$age <=39,"young",
                                             ifelse(dat_health$age >=40 & dat_health$age <=59,"middle",
                                                    ifelse(dat_health$age >=60,"old",NA)))))
name_all <- colnames(dat_health)
rownames(dat_health) <- dat_health$sample_id
bac_species_name <- name_all[bac_species]
bac_phylum_name <- name_all[bac_phylum]

phylum_abundance_plot <- function(data,palette,name,height,width){
  results = as.data.frame(matrix(NA,ncol = 2))
  colnames(results) = c("phylum","mean")
  phylum = colnames(data)
  for (i in 1:length(phylum)) {
    results[i,1]=phylum[i]
    results[i,2]=mean(data[,phylum[i]])
  }
  results$phylum = results$phylum[order(results$mean,decreasing = T)]
  phylum_top6 = results$phylum[1:6]
  data_top6=data[,phylum_top6]
  data_top6$Others=1-rowSums(data_top6)
  # data_top6=data_top6[order(data_top6[,phylum_top6[1]],decreasing = T),]
  for (i in rev(1:length(phylum_top6))){
    data_top6=data_top6[order(data_top6[,phylum_top6[i]],decreasing = T),]
  }
  data_top6$ID=c(1:nrow(data_top6))
  library(reshape2)
  data_top6_mutate = melt(data_top6,id.vars="ID",variable.name="phylum",value.name = "abundance")
  data_top6_mutate$phylum=factor(data_top6_mutate$phylum,levels = rev(c(phylum_top6,"Others")))
  library(ggplot2)
  p <- ggplot(data_top6_mutate, aes(ID,  abundance, fill = phylum)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = palette ) +
    labs(x = 'Individual ID', y = 'Relative Abundance (%)', fill = "") +
    theme_classic()+
    guides(fill=guide_legend(nrow=7,reverse=TRUE))
  ggsave(paste0("RelativeAbundance_healthy_",name,".pdf"),height = height,width =width)
}
palette <- rev(c("#91bde2", "#83c5b0","#8fb79d","#eadfc3","#eec16c","#d75032","#e8edf3"))
phylum_abundance_plot(dat_health[,name_all[bac_phylum]],palette,"phylum.all",5,8)

mean(dat_health$Firmicutes)
mean(dat_health$Bacteroidetes)
mean(dat_health$Actinobacteria)
mean(dat_health$Proteobacteria)

# Table S1-------------------------------------------------- --------------------------------------------
# sample size
group_by(dat_health[,c("age","BMI","country")],country) %>% summarise_each(funs(length))

# mean and SD of age
dat_health_2 <- subset(dat_health,!is.na(age))
group_by(dat_health_2[,c("age","country")],country) %>% summarise_each(funs(mean,sd))

# mean and SD of BMI
dat_health_2 <- subset(dat_health,!is.na(BMI))
group_by(dat_health_2[,c("BMI","country")],country) %>% summarise_each(funs(mean,sd))
# Westernization和PMID信息在RelAbund_batchCorrect.xlsx对应non_westernized、PMID两列，统计分析手工进行

# Figure 1C; proportion of microbial variation among healthy women------------------------------ ----------------------------------

library(vegan)
library(ggplot2)
bray_dis <- as.matrix(vegdist(dat_health[,name_all[bac_species]],method = "bray",diag = TRUE,upper = TRUE))

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
points_2 <- cbind(points, dat_health$country)
colnames(points_2)[4] <- "country"

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ country+DNA_extraction_kit+sequencing_platform, data = dat_health,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ country, data = dat_health,permutations = 999)
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
ggsave("bray_country_after.pdf", p, width = 8, height = 6)

# continent -----------------------------------------------------------
points_2 <- cbind(points, dat_health$continent)
colnames(points_2)[4] <- "continent"

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ continent+DNA_extraction_kit+sequencing_platform, data = dat_health,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ continent, data = dat_health,permutations = 999)
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
ggsave("bray_continent_after.pdf", p, width = 8, height = 6)

# age -----------------------------------------------------------
points_2 <- cbind(points, dat_health$age)
colnames(points_2)[4] <- "age"
points_2 <- subset(points_2,!is.na(age))
idx <- rownames(bray_dis) %in% rownames(points_2)

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis[idx,idx] ~ age+DNA_extraction_kit+sequencing_platform, data = dat_health[idx,],permutations = 999)}
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
ggsave("bray_age_after.pdf", p, width = 8, height = 6)

# non_westernized -----------------------------------------------------------
points_2 <- cbind(points, dat_health$non_westernized)
colnames(points_2)[4] <- "non_westernized"

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis ~ non_westernized+DNA_extraction_kit+sequencing_platform, data = dat_health,permutations = 999)}
# group_permanova <- adonis(bray_dis ~ non_westernized, data = dat_health,permutations = 999)
a <- round(group_permanova$aov.tab$R2[1],digits = 3)
R2 <- paste("adonis R2: ",a, sep = "")
b <- group_permanova$aov.tab$`Pr(>F)`[1]
p_v <- paste("p: ",b, sep = "")
title <- paste(R2," ",p_v, sep = "");title # adonis R2: 0.007 p: 0.001

p <- ggplot(points_2, aes(x=x, y=y, color=non_westernized)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  # stat_ellipse(level=0.95) +
  # geom_encircle(aes(fill=non_westernized), alpha = 0.1, show.legend = F) +
  theme_bw()
ggsave("bray_non_westernized_after.pdf", p, width = 8, height = 6)

# BMI -----------------------------------------------------------
points_2 <- cbind(points, dat_health$BMI)
colnames(points_2)[4] <- "BMI"
points_2 <- subset(points_2,!is.na(BMI))
idx <- rownames(bray_dis) %in% rownames(points_2)

set.seed(2023)
group_permanova <- foreach(
  n=1,.combine=cbind,.packages = 'vegan') %dopar% {
    adonis(bray_dis[idx,idx] ~ BMI+DNA_extraction_kit+sequencing_platform, data = dat_health[idx,],permutations = 999)}
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
ggsave("bray_BMI_after.pdf", p, width = 8, height = 6)

# radarchart 
plotdata1 <- read.xlsx("plotdata.xlsx",sheet = 1,rowNames = T)
library(fmsb)
pdf("radarplot.pdf",height = 5,width = 5)
radarchart(
  plotdata1, axistype = 1,
  # Customize the polygon
  pcol = "#ff5a5f", plwd = 2, plty = "solid",
  # Customize the grid
  cglcol = "grey", cglty = "dashed", cglwd = 1,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels
  vlcex = 1.1, vlabels = colnames(plotdata1),
  caxislabels = c("0", "2.5%", "5.0%", "7.5%", "10.0%"))
dev.off()


# Figure 1D;bootstrapped_PERMANOVA-------------------------------------------------- -------

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
        adonis(data_sample ~ variable_select+DNA_extraction_kit+sequencing_platform, data = metadata_sample,permutations = 999)
      }
    results[i,1]=i
    results[i,2]=round(group_permanova$aov.tab$R2[1],digits = 3)
    results[i,3]=group_permanova$aov.tab$`Pr(>F)`[1]
  }
  return(results)
}

result_country <- bootstrap_PERMANOVA(bray_dis,dat_health,"country",100,c(0.2,0.8))
write.xlsx(result_country,"result_country.xlsx")

result_continent <- bootstrap_PERMANOVA(bray_dis,dat_health,"continent",100,c(0.2,0.8))
write.xlsx(result_continent,"result_continent.xlsx")

points_2 <- cbind(points, dat_health$age)
colnames(points_2)[4] <- "age"
points_2 <- subset(points_2,!is.na(age))
idx <- rownames(bray_dis) %in% rownames(points_2)
result_age <- bootstrap_PERMANOVA(bray_dis[idx,idx],dat_health[idx,],"age",100,c(0.2,0.8))
write.xlsx(result_age,"result_age.xlsx")

result_non_westernized <- bootstrap_PERMANOVA(bray_dis,dat_health,"non_westernized",100,c(0.2,0.8))
write.xlsx(result_non_westernized,"result_non_westernized.xlsx")

points_2 <- cbind(points, dat_health$BMI)
colnames(points_2)[4] <- "BMI"
points_2 <- subset(points_2,!is.na(BMI))
idx <- rownames(bray_dis) %in% rownames(points_2)
result_BMI <- bootstrap_PERMANOVA(bray_dis[idx,idx],dat_health[idx,],"BMI",100,c(0.2,0.8))
write.xlsx(result_BMI,"result_BMI.xlsx")

bootstrapped_permanova <- read.xlsx("bootstrapped_permanova.xlsx",sheet = 1)
bootstrapped_permanova_reshape <- melt(bootstrapped_permanova,id.vars="ID",variable.name = "group",value.name = "r2")
bootstrapped_permanova_reshape$group <- factor(bootstrapped_permanova_reshape$group,levels = c("Country","Age","Continent","Non_westernization","BMI"),labels = c("Country","Age","Continent","Westernization","BMI"))

p <- ggplot(bootstrapped_permanova_reshape,aes(group,r2))+
  stat_boxplot(geom = "errorbar",aes(color=group),width=0.2)+
  geom_boxplot(aes(color=group),outlier.alpha = 0,width=0.6)+
  geom_jitter(aes(fill=group),width = 0.2,size=1.4,shape=21,color="black",alpha=0.6)+
  theme_bw()+
  scale_color_manual(values = c("#336ba2","#f13333","#33adc3","#68c466","#a87eb2"))+
  scale_fill_manual(values = c("#336ba2","#f13333","#33adc3","#68c466","#a87eb2"))+
  theme(axis.text = element_text(color = "black",size = 8))+
  labs(x="",y="R2")+
  theme(axis.title = element_text(color="black",size=15),
        axis.text.x = element_text(color = "black",size = 12,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_text(color = "black",size = 12))
p
ggsave("bootstrapped_permanova.pdf",width=7,height=6)

# Figure 1E,Table S3; contribution of age and region to individual microbiota-------------------------- -----------------

species_contribution <- function(data,variable){
  results=data.frame(matrix(NA,ncol = 3))
  colnames(results)=c("species",paste0("adj.R2",variable),"p")
  for (i in 1:length(bac_species_name)){
    model=lm(data[,bac_species_name[i]]~data[,variable],data = data)
    model_summary=summary(model)
    results[i,1]=bac_species_name[i]
    results[i,2]=model_summary$adj.r.squared
    results[i,3]=pf(model_summary$fstatistic[1],model_summary$fstatistic[2],model_summary$fstatistic[3],lower.tail = F)
  }
  results$p.adj=p.adjust(results$p)
  results[,paste0("adj.R2",variable)]=ifelse(results$p.adj<0.05,results[,paste0("adj.R2",variable)],0)
  return(results)
}

library(compositions)
# age
dat_health_2 <- subset(dat_health,!is.na(age))
dat_health_2[,bac_species_name] <- clr(dat_health_2[,bac_species_name])
species_contribution_age <- species_contribution(dat_health_2,"age")

# country
dat_health_2 <- dat_health
dat_health_2[,bac_species_name] <- clr(dat_health_2[,bac_species_name])
species_contribution_country <- species_contribution(dat_health_2,"country")

# continent
dat_health_2 <- dat_health
dat_health_2[,bac_species_name] <- clr(dat_health_2[,bac_species_name])
species_contribution_continent <- species_contribution(dat_health_2,"continent")

# BMI
dat_health_2 <- subset(dat_health,!is.na(BMI))
dat_health_2[,bac_species_name] <- clr(dat_health_2[,bac_species_name])
species_contribution_BMI <- species_contribution(dat_health_2,"BMI")

# non_westernized
dat_health_2 <- dat_health
dat_health_2$non_westernized <- factor(dat_health_2$non_westernized,levels = c("yes","no"),labels = c(0,1))
dat_health_2[,bac_species] <- clr(dat_health_2[,bac_species])
species_contribution_non_westernized <- species_contribution(dat_health_2,"non_westernized")

species_contribution_all <- Reduce(function(x,y) merge(x,y,by="species",all=T),list(species_contribution_age[,c(1,2)],species_contribution_continent[,c(1,2)],species_contribution_country[,c(1,2)],species_contribution_BMI[,c(1,2)],species_contribution_non_westernized[,c(1,2)]))
species_contribution_all$adj.R2.sum <- rowSums(species_contribution_all[,2:6])
species_contribution_all <- subset(species_contribution_all,adj.R2.sum>0.05) # species with R2 >5% remained
colnames(species_contribution_all)[2:6] <- c("Age","Continent","Country","BMI","Westernization")
species_contribution_all$category <- colnames(species_contribution_all[,2:6])[max.col(species_contribution_all[,2:6])]
table(species_contribution_all$category)
write.csv(species_contribution_all,"species_contribution_all.csv",quote = F,row.names = F) # refering to Table S3

# barplot
species_contribution_all <- read.csv("species_contribution_all.csv")
species_contribution_all <- arrange(species_contribution_all,category,desc(adj.R2.sum))
species_contribution_all$ID <- 1:dim(species_contribution_all)[1]
species_contribution_all_mutate <- melt(species_contribution_all,measure.vars = c("Age","Continent","Country","BMI","Westernization"),variable.name = "group",value.name = "r2")
species_contribution_all_mutate$group <- factor(species_contribution_all_mutate$group,levels = rev(c("Age","Country","Continent","Westernization","BMI")))

p <- ggplot(species_contribution_all_mutate, aes(ID,  r2, fill = group)) +
  geom_col(position = "stack", width = 1,color="white",linewidth=0.2) +
  scale_fill_manual(values = rev(c("#f13333","#336ba2","#33adc3","#68c466","#a87eb2"))) +
  labs(x = '78 species (adjusted R2 >5%)', y = 'Adjusted R2', fill = "") +
  theme_classic()+
  theme(axis.title = element_text(size=12,color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(linewidth = 0),
        axis.text.y = element_text(size=12,color = "black"))+
  guides(fill=guide_legend(nrow=5,reverse=TRUE))
p
ggsave("species_contribution_all.pdf",width = 15,height = 6)

species_contribution_all_mutate_country <- subset(species_contribution_all_mutate,category=="Country")
species_contribution_all_mutate_country$group <- factor(species_contribution_all_mutate_country$group,levels = rev(c("Country","Age","Continent","Westernization","BMI")))
p <- ggplot(species_contribution_all_mutate_country, aes(ID,  r2, fill = group)) +
  geom_col(position = "stack", width = 1,color="white",linewidth=0.2) +
  scale_fill_manual(values = rev(c("#336ba2","#f13333","#33adc3","#68c466","#a87eb2"))) +
  labs(x = '78 species (adjusted R2 >5%)', y = 'Adjusted R2', fill = "") +
  theme_classic()+
  theme(axis.title = element_text(size=12,color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(linewidth = 0),
        axis.text.y = element_text(size=12,color = "black"))+
  guides(fill=guide_legend(nrow=5,reverse=TRUE))
p
ggsave("species_contribution_country.pdf",width = 9,height = 6)

# phylum barplot
species_contribution_all <- merge(species_contribution_all,short_phylum,by.x="species",by.y="short",all=F)
species_contribution_all$phylum <- factor(species_contribution_all$phylum,levels = c("Actinobacteria","Bacteroidetes","Euryarchaeota","Firmicutes","Proteobacteria","Verrucomicrobia"))
species_contribution_all$bar <- rep(1,dim(species_contribution_all)[1])

p <- ggplot(species_contribution_all, aes(ID,  bar, fill = phylum)) +
  geom_col(position = "stack", width = 1,color="white",linewidth=0.2) +
  scale_fill_manual(values = c("#8fb79d","#83c5b0","#d75032","#91bde2","#eadfc3","#eec16c")) +
  labs(x = '78 species (adjusted R2 >5%)', y = 'Adjusted R2', fill = "") +
  theme_classic()+
  theme(axis.title = element_text(size=12,color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(linewidth = 0),
        axis.text.y = element_text(size=12,color = "black"))+
  guides(fill=guide_legend(nrow=5,reverse=TRUE))
p
ggsave("species_contribution_bar.pdf",width = 15,height = 0.6)


# include country > 0.05*total number of subject (0.05*3644=182.2)
GBR <- rownames(dat_health)[dat_health$country == "GBR"]
CHN <- rownames(dat_health)[dat_health$country == "CHN"]
NLD <- rownames(dat_health)[dat_health$country == "NLD"]
SWE <- rownames(dat_health)[dat_health$country == "SWE"]
USA <- rownames(dat_health)[dat_health$country == "USA"]

Africa <- rownames(dat_health)[dat_health$continent == "Africa"]
Asia <- rownames(dat_health)[dat_health$continent == "Asia"]
Europe <- rownames(dat_health)[dat_health$continent == "Europe"]
North_America <- rownames(dat_health)[dat_health$continent == "North America"]

infant <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age <=1]
child_adolescent <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=2 & dat_health$age <=17]
young <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=18 & dat_health$age <=39]
middle <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=40 & dat_health$age <=59]
old <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=60]

# Table S4; Comparison of relative abundance of species across age groups ---------------------------------------------------------------
dat_health_2 <- subset(dat_health,!is.na(age_group))
AgeGroupsComparisons <- as.data.frame(matrix(NA,ncol = 2))
colnames(AgeGroupsComparisons) <- c("species","kruskal_p")
for (i in 1:length(bac_species_name)){
  comparison=kruskal.test(dat_health_2[,bac_species_name[i]]~age_group,data = dat_health_2)
  AgeGroupsComparisons[i,1]=bac_species_name[i]
  AgeGroupsComparisons[i,2]=comparison$p.value
}
AgeGroupsComparisons$p.adjust <- p.adjust(AgeGroupsComparisons$kruskal_p,method = "BH")
write.csv(AgeGroupsComparisons,"AgeGroupsComparisons.csv",row.names = F,quote = F)

# include country > 0.05*total number of subject (0.05*3644=182.2)
GBR <- rownames(dat_health)[dat_health$country == "GBR"]
CHN <- rownames(dat_health)[dat_health$country == "CHN"]
NLD <- rownames(dat_health)[dat_health$country == "NLD"]
SWE <- rownames(dat_health)[dat_health$country == "SWE"]
USA <- rownames(dat_health)[dat_health$country == "USA"]

Africa <- rownames(dat_health)[dat_health$continent == "Africa"]
Asia <- rownames(dat_health)[dat_health$continent == "Asia"]
Europe <- rownames(dat_health)[dat_health$continent == "Europe"]
North_America <- rownames(dat_health)[dat_health$continent == "North America"]

infant <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age <=1]
child_adolescent <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=2 & dat_health$age <=17]
young <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=18 & dat_health$age <=39]
middle <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=40 & dat_health$age <=59]
old <- rownames(dat_health)[!is.na(dat_health$age) & dat_health$age >=60]

# Figure 2a; pie chart of phylum proportion across age groups and regions------------------------ ------------

pie_chart <- function(data,region,age){
  summary_data=data.frame(matrix(NA,ncol = 2))
  colnames(summary_data)=c("phylum","mean")
  col_data=colnames(data)
  for (i in 1:length(col_data)){
    summary_data[i,1]=col_data[i]
    summary_data[i,2]=mean(data[,col_data[i]])
  }
  library(ggplot2)
  p=ggplot(summary_data,aes(x="",y=mean,fill=phylum))+
    geom_bar(stat = "identity",width = 1)+
    coord_polar("y",start = 0)+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          panel.grid = element_blank())+
    scale_fill_manual(values = c("#8fb79d","#83c5b0","#d75032","#91bde2","#eadfc3","#eec16c"))
  ggsave(paste0(region,"_",age,"_piechart.pdf"),width = 7.5,height = 5)
}

# Asia
dat_health_2 <- dat_health[intersect(infant,Asia),bac_phylum_name]
pie_chart(dat_health_2,"Asia","infant")

dat_health_2 <- dat_health[intersect(child_adolescent,Asia),bac_phylum_name]
pie_chart(dat_health_2,"Asia","child_adolescent")

dat_health_2 <- dat_health[intersect(young,Asia),bac_phylum_name]
pie_chart(dat_health_2,"Asia","young")

dat_health_2 <- dat_health[intersect(middle,Asia),bac_phylum_name]
pie_chart(dat_health_2,"Asia","middle")

dat_health_2 <- dat_health[intersect(old,Asia),bac_phylum_name]
pie_chart(dat_health_2,"Asia","old")

# Africa
dat_health_2 <- dat_health[intersect(young,Africa),bac_phylum_name]
pie_chart(dat_health_2,"Africa","young")

dat_health_2 <- dat_health[intersect(middle,Africa),bac_phylum_name]
pie_chart(dat_health_2,"Africa","middle")

# Europe
dat_health_2 <- dat_health[intersect(infant,Europe),bac_phylum_name]
pie_chart(dat_health_2,"Europe","infant")

dat_health_2 <- dat_health[intersect(child_adolescent,Europe),bac_phylum_name]
pie_chart(dat_health_2,"Europe","child_adolescent")

dat_health_2 <- dat_health[intersect(young,Europe),bac_phylum_name]
pie_chart(dat_health_2,"Europe","young")

dat_health_2 <- dat_health[intersect(middle,Europe),bac_phylum_name]
pie_chart(dat_health_2,"Europe","middle")

dat_health_2 <- dat_health[intersect(old,Europe),bac_phylum_name]
pie_chart(dat_health_2,"Europe","old")

# North_America
dat_health_2 <- dat_health[intersect(child_adolescent,North_America),bac_phylum_name]
pie_chart(dat_health_2,"North_America","child_adolescent")

dat_health_2 <- dat_health[intersect(young,North_America),bac_phylum_name]
pie_chart(dat_health_2,"North_America","young")

dat_health_2 <- dat_health[intersect(middle,North_America),bac_phylum_name]
pie_chart(dat_health_2,"North_America","middle")

dat_health_2 <- dat_health[intersect(old,North_America),bac_phylum_name]
pie_chart(dat_health_2,"North_America","old")

# SWE
dat_health_2 <- dat_health[intersect(infant,SWE),bac_phylum_name]
pie_chart(dat_health_2,"SWE","infant")

dat_health_2 <- dat_health[intersect(young,SWE),bac_phylum_name]
pie_chart(dat_health_2,"SWE","young")

dat_health_2 <- dat_health[intersect(old,SWE),bac_phylum_name]
pie_chart(dat_health_2,"SWE","old")

# CHN
dat_health_2 <- dat_health[intersect(young,CHN),bac_phylum_name]
pie_chart(dat_health_2,"CHN","young")

dat_health_2 <- dat_health[intersect(middle,CHN),bac_phylum_name]
pie_chart(dat_health_2,"CHN","middle")

dat_health_2 <- dat_health[intersect(old,CHN),bac_phylum_name]
pie_chart(dat_health_2,"CHN","old")

# NLD
dat_health_2 <- dat_health[intersect(young,NLD),bac_phylum_name]
pie_chart(dat_health_2,"NLD","young")

dat_health_2 <- dat_health[intersect(middle,NLD),bac_phylum_name]
pie_chart(dat_health_2,"NLD","middle")

dat_health_2 <- dat_health[intersect(old,NLD),bac_phylum_name]
pie_chart(dat_health_2,"NLD","old")

# GBR
dat_health_2 <- dat_health[intersect(infant,GBR),bac_phylum_name]
pie_chart(dat_health_2,"GBR","infant")

dat_health_2 <- dat_health[intersect(young,GBR),bac_phylum_name]
pie_chart(dat_health_2,"GBR","young")

dat_health_2 <- dat_health[intersect(middle,GBR),bac_phylum_name]
pie_chart(dat_health_2,"GBR","middle")

dat_health_2 <- dat_health[intersect(old,GBR),bac_phylum_name]
pie_chart(dat_health_2,"GBR","old")

# USA
dat_health_2 <- dat_health[intersect(child_adolescent,USA),bac_phylum_name]
pie_chart(dat_health_2,"USA","child_adolescent")

dat_health_2 <- dat_health[intersect(young,USA),bac_phylum_name]
pie_chart(dat_health_2,"USA","young")

dat_health_2 <- dat_health[intersect(middle,USA),bac_phylum_name]
pie_chart(dat_health_2,"USA","middle")

dat_health_2 <- dat_health[intersect(old,USA),bac_phylum_name]
pie_chart(dat_health_2,"USA","old")

# Figure 2b; alpha diversity comparisons among group----------------------------------------------- ----------------------------------

shannon_index <- data.frame(vegan::diversity(dat_health[,bac_species_name],index = "shannon"))
alpha_index <- shannon_index
colnames(alpha_index) <- c("shannon")
alpha_index$sample_id <- rownames(alpha_index)
alpha_index <- merge(alpha_index,dat_health[,c("sample_id","age","age_group","country","continent")],by="sample_id",all=F)
alpha_index_select <- subset(alpha_index,!is.na(age_group))
alpha_index_select$age_group <- factor(alpha_index_select$age_group,levels = c("infant","child_adolescent","young","middle","old"),labels = c("Infant","Child & Adolescent","Young","Middle-aged","Elderly"))

# shannon 
p <- ggviolin(alpha_index_select,x="age_group",y="shannon",add=c("boxplot","dotplot"),add.params = list(fill="white",binwidth=0.025,color="age_group",alpha=0.2,width=0.1,size=0.6),
              fill = "white",color = "age_group",width = 1,size = 1.2,position = "dodge",trim = T)+
  geom_pwc(aes(group=age_group),method = "wilcox.test",p.adjust.method="BH",label = "p.adj.signif")+
  # stat_compare_means(comparisons = boxplot_comparison,method = "wilcox.test",p.adjust.methods="BH")+
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = "Control") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(color="black",size=17),
        axis.text = element_text(color="black",size=15))+
  scale_color_manual(values = c("#cdbe6a","#00a087","#86a4cf","#e64b35","#82491e"))+
  labs(x="",y="Shannon index")
p
ggsave("age_group_shannon_comparison2.pdf",width = 7.5,height = 5)

# Figure 2c;topological structure comparison of different age --------------------------------------

library(igraph)
library(psych)

# network parameter
# infant
dat_health_2 <- subset(dat_health,age_group=="infant")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0 
occor.igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
occor.bad.vs = V(occor.igraph)[igraph::degree(occor.igraph) == 0]
occor.igraph = delete.vertices(occor.igraph, occor.bad.vs)
network.occor <- as.data.frame(matrix(NA,ncol = 10))
names(network.occor) <- c("group","num.edges","num.vertices","connectance","average.degree","edge.connectivity","Trajectorying.coefficient","notTrajectorying.coefficient","centralization.betweenness","centralization.degree")
network.occor[1,1] = "infant"
network.occor[1,2] = length(E(occor.igraph))
network.occor[1,3] = length(V(occor.igraph))
network.occor[1,4] = edge_density(occor.igraph,loops=FALSE)
network.occor[1,5] = mean(igraph::degree(occor.igraph))
network.occor[1,6] = edge_connectivity(occor.igraph)
network.occor[1,7] = transitivity(occor.igraph)
network.occor[1,8] = no.clusters(occor.igraph)
network.occor[1,9] = centralization.betweenness(occor.igraph)$centralization 
network.occor[1,10] = centralization.degree(occor.igraph)$centralization

# child_adolescent
dat_health_2 <- subset(dat_health,age_group=="child_adolescent")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0 
occor.igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
occor.bad.vs = V(occor.igraph)[igraph::degree(occor.igraph) == 0]
occor.igraph = delete.vertices(occor.igraph, occor.bad.vs)
network.occor[2,1] = "child_adolescent"
network.occor[2,2] = length(E(occor.igraph))
network.occor[2,3] = length(V(occor.igraph))
network.occor[2,4] = edge_density(occor.igraph,loops=FALSE)
network.occor[2,5] = mean(igraph::degree(occor.igraph))
network.occor[2,6] = edge_connectivity(occor.igraph)
network.occor[2,7] = transitivity(occor.igraph)
network.occor[2,8] = no.clusters(occor.igraph)
network.occor[2,9] = centralization.betweenness(occor.igraph)$centralization 
network.occor[2,10] = centralization.degree(occor.igraph)$centralization

# young
dat_health_2 <- subset(dat_health,age_group=="young")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0 
occor.igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
occor.bad.vs = V(occor.igraph)[igraph::degree(occor.igraph) == 0]
occor.igraph = delete.vertices(occor.igraph, occor.bad.vs)
network.occor[3,1] = "young"
network.occor[3,2] = length(E(occor.igraph))
network.occor[3,3] = length(V(occor.igraph))
network.occor[3,4] = edge_density(occor.igraph,loops=FALSE)
network.occor[3,5] = mean(igraph::degree(occor.igraph))
network.occor[3,6] = edge_connectivity(occor.igraph)
network.occor[3,7] = transitivity(occor.igraph)
network.occor[3,8] = no.clusters(occor.igraph)
network.occor[3,9] = centralization.betweenness(occor.igraph)$centralization 
network.occor[3,10] = centralization.degree(occor.igraph)$centralization

# middle
dat_health_2 <- subset(dat_health,age_group=="middle")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0 
occor.igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
occor.bad.vs = V(occor.igraph)[igraph::degree(occor.igraph) == 0]
occor.igraph = delete.vertices(occor.igraph, occor.bad.vs)
network.occor[4,1] = "middle"
network.occor[4,2] = length(E(occor.igraph))
network.occor[4,3] = length(V(occor.igraph))
network.occor[4,4] = edge_density(occor.igraph,loops=FALSE)
network.occor[4,5] = mean(igraph::degree(occor.igraph))
network.occor[4,6] = edge_connectivity(occor.igraph)
network.occor[4,7] = transitivity(occor.igraph)
network.occor[4,8] = no.clusters(occor.igraph)
network.occor[4,9] = centralization.betweenness(occor.igraph)$centralization 
network.occor[4,10] = centralization.degree(occor.igraph)$centralization

# old
dat_health_2 <- subset(dat_health,age_group=="old")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.4] = 0 
occor.igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
occor.bad.vs = V(occor.igraph)[igraph::degree(occor.igraph) == 0]
occor.igraph = delete.vertices(occor.igraph, occor.bad.vs)
network.occor[5,1] = "old"
network.occor[5,2] = length(E(occor.igraph))
network.occor[5,3] = length(V(occor.igraph))
network.occor[5,4] = edge_density(occor.igraph,loops=FALSE)
network.occor[5,5] = mean(igraph::degree(occor.igraph))
network.occor[5,6] = edge_connectivity(occor.igraph)
network.occor[5,7] = transitivity(occor.igraph)
network.occor[5,8] = no.clusters(occor.igraph)
network.occor[5,9] = centralization.betweenness(occor.igraph)$centralization 
network.occor[5,10] = centralization.degree(occor.igraph)$centralization
write.csv(network.occor,"network.age_group.csv",quote = F,row.names = F)

# visualization
library(igraph)
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

# infant
dat_health_2 <- subset(dat_health,age_group=="infant")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_df <- CorrDF(occor$r , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.4),] 
cor_df <- cor_df[which(cor_df$p < 0.05),]

igraph <- graph_from_data_frame(cor_df, direct=F)
V(igraph)$size <- 4
V(igraph)$color <- "grey" 
E(igraph)$color[E(igraph)$cor >= 0.4] <- "red" 
E(igraph)$color[E(igraph)$cor <= -0.4] <- "lightblue" 
E(igraph)$width <- 1
coords <- layout_with_fr(igraph, niter=999,grid="nogrid") 
pdf("infant_igraph.pdf", height = 5, width = 5)
plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) 
dev.off()
# 每次成图时，关联的排列位置会有变化，但关联的组成情况不变

# child_adolescent
dat_health_2 <- subset(dat_health,age_group=="child_adolescent")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_df <- CorrDF(occor$r , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.4),] 
cor_df <- cor_df[which(cor_df$p < 0.05),]

igraph <- graph_from_data_frame(cor_df, direct=F)
V(igraph)$size <- 4
V(igraph)$color <- "grey" 
E(igraph)$color[E(igraph)$cor >= 0.4] <- "red" 
E(igraph)$color[E(igraph)$cor <= -0.4] <- "lightblue" 
E(igraph)$width <- 1
coords <- layout_with_fr(igraph, niter=999,grid="nogrid") 
pdf("child_adolescent_igraph.pdf", height = 5, width = 5)
plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) 
dev.off()
# 每次成图时，关联的排列位置会有变化，但关联的组成情况不变

# young
dat_health_2 <- subset(dat_health,age_group=="young")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_df <- CorrDF(occor$r , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.4),] 
cor_df <- cor_df[which(cor_df$p < 0.05),]

igraph <- graph_from_data_frame(cor_df, direct=F)
V(igraph)$size <- 4
V(igraph)$color <- "grey" 
E(igraph)$color[E(igraph)$cor >= 0.4] <- "red" 
E(igraph)$color[E(igraph)$cor <= -0.4] <- "lightblue" 
E(igraph)$width <- 1
coords <- layout_with_fr(igraph, niter=999,grid="nogrid") 
pdf("young_igraph.pdf", height = 5, width = 5)
plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) 
dev.off()
# 每次成图时，关联的排列位置会有变化，但关联的组成情况不变

# middle
dat_health_2 <- subset(dat_health,age_group=="middle")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_df <- CorrDF(occor$r , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.4),] 
cor_df <- cor_df[which(cor_df$p < 0.05),]

igraph <- graph_from_data_frame(cor_df, direct=F)
V(igraph)$size <- 4
V(igraph)$color <- "grey" 
E(igraph)$color[E(igraph)$cor >= 0.4] <- "red" 
E(igraph)$color[E(igraph)$cor <= -0.4] <- "lightblue" 
E(igraph)$width <- 1
coords <- layout_with_fr(igraph, niter=999,grid="nogrid") 
pdf("middle_igraph.pdf", height = 5, width = 5)
plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) 
dev.off()
# 每次成图时，关联的排列位置会有变化，但关联的组成情况不变

# old
dat_health_2 <- subset(dat_health,age_group=="old")[,bac_species_name]
dat_health_2 <- data.frame(clr(dat_health_2))
occor <- corr.test(dat_health_2,use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_df <- CorrDF(occor$r , occor$p) 
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.4),] 
cor_df <- cor_df[which(cor_df$p < 0.05),]

igraph <- graph_from_data_frame(cor_df, direct=F)
V(igraph)$size <- 4
V(igraph)$color <- "grey" 
E(igraph)$color[E(igraph)$cor >= 0.4] <- "red" 
E(igraph)$color[E(igraph)$cor <= -0.4] <- "lightblue" 
E(igraph)$width <- 1
coords <- layout_with_fr(igraph, niter=999,grid="nogrid") 
pdf("old_igraph.pdf", height = 5, width = 5)
plot(igraph, layout=coords, vertex.label = NA, vertex.frame.color=NA) 
dev.off()
# 每次成图时，关联的排列位置会有变化，但关联的组成情况不变

# Figure 2d;species distribution by age ----------------------- ----------------------------------------

results<-as.data.frame(matrix(NA,ncol = 2))
colnames(results)<-c("species","mean")
for (i in 1:length(bac_species_name)) {
  results[i,1]=bac_species_name[i]
  results[i,2]=mean(dat_health[,bac_species_name[i]])
}
results<-results[order(results$mean,decreasing = T),]
species_top10 <- results$species[1:10]

# age_group
dat_health_2 <- subset(dat_health,!is.na(age_group))
dat_health_2 <- dat_health_2[,c(species_top10,"age_group")]
dat_health_2 <- melt(dat_health_2,measure.vars = species_top10,variable.name = "species",value.name = "abundance")
dat_health_2$species <- factor(dat_health_2$species,levels = species_top10)
dat_health_2$age_group <- factor(dat_health_2$age_group,levels = rev(c("infant","child_adolescent","young","middle","old")),labels = rev(c("Infant","Child & Adolescent","Young adult","Middle-aged adult","Elderly")))

p <- ggplot(dat_health_2,aes(x = abundance,y = age_group ,fill=age_group))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10))+
  geom_violin(scale = "width",draw_quantiles = c(0.25,0.5,0.75),size=0.1)+
  # geom_dotplot(aes(size = abundance))+
  scale_fill_manual(values = rev(c("#cdbe6a","#00a087","#86a4cf","#e64b35","#82491e")))+
  facet_grid(cols = vars(species))+
  theme(legend.position = "bottom",
        panel.grid = element_blank(), #移除背景网格线
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 10,color="black"),
        strip.background = element_blank(), #移除分面外围边框
        strip.text.x = element_text(size = 10, angle = 60,vjust = 0,hjust = 0) #分面文本标签倾斜60°
  ) +
  labs(x="Relative abundance (%)",y="")+
  scale_x_continuous(limits = c(0,1),breaks = c(0,0.5,1),labels = c("0","50%","100%"))
p
ggsave("species_age_group_distribution.pdf",width = 16,height = 8)

# Figure S2a,species distribution by continent and age ---------------------------------------------------------------
dat_health_2 <- subset(dat_health,!is.na(dat_health$age))
dat_health_2 <- dat_health_2[,c(species_top10,"continent","age")]
dat_health_2 <- melt(dat_health_2,measure.vars = species_top10,variable.name = "species",value.name = "abundance")
dat_health_2$species <- factor(dat_health_2$species,levels = rev(species_top10))
dat_health_2$continent <- factor(dat_health_2$continent,levels = c("Africa","Asia","Europe","North America"))

p <- ggplot(dat_health_2,aes(x = age,y = species,color=continent))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10))+
  #geom_point(size = (trim$abundance*10)^1.5/8,alpha = .3)+
  geom_point(aes(size = abundance),shape = 1,alpha=0.5)+
  scale_color_manual(values = c("#7ea3b5","#ffe09c","#f05ea2","#fec7b2"))+
  # scale_size(range=c(1,5))+
  facet_grid(.~continent,scales = "free")+
  theme(legend.position = "top")+
  labs(x="Age (year)",y="")
p
ggsave("species_continent_age_distribution.pdf",width = 10,height = 8)

# Figure S2b,species distribution by country and age ---------------------------------------------------------------
dat_health_2 <- subset(dat_health,!is.na(dat_health$age))
dat_health_2 <- subset(dat_health_2,dat_health_2$country=="GBR"|dat_health_2$country=="CHN"|dat_health_2$country=="NLD"|dat_health_2$country=="SWE"|dat_health_2$country=="USA")
dat_health_2 <- dat_health_2[,c(species_top10,"country","age")]
dat_health_2 <- melt(dat_health_2,measure.vars = species_top10,variable.name = "species",value.name = "abundance")
dat_health_2$species <- factor(dat_health_2$species,levels = rev(species_top10))
dat_health_2$country <- factor(dat_health_2$country,levels = c("CHN","NLD","SWE","GBR","USA"),labels = c("China","Netherlands","Sweden","the UK","the USA"))

p <- ggplot(dat_health_2,aes(x = age,y = species,color=country))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10))+
  #geom_point(size = (trim$abundance*10)^1.5/8,alpha = .3)+
  geom_point(aes(size = abundance),shape = 1,alpha=0.5)+
  scale_color_manual(values = c("#1a80bd","#379160","#8684b9","#c3503f","#e4933d"))+
  # scale_size(range=c(1,5))+
  facet_grid(.~country,scales = "free")+
  theme(legend.position = "top")+
  labs(x="Age (year)",y="")
p
ggsave("species_country_age_distribution.pdf",width = 12,height = 8)

# Figure 3a, c; iterative RF-------------------------------------------------------- -------------------

AgeCohort <- c(infant,child_adolescent,young,middle,old)
set.seed(2023)
idx <- sample(c(0,1),length(AgeCohort),replace = T,prob = c(0.5,0.5))
AgeCohortTrain <- AgeCohort[idx==0]
AgeCohortTest <- AgeCohort[idx==1]

rfPermutedMultipleNull_infant <- iterative_rf_multiple_with_null_distribution(clr(dat_health[,bac_species_name]),
                                                                              AgeCohortTrain,AgeCohortTest,infant,c(child_adolescent,young,middle,old),200,10,round(0.3*length(infant),digits = 0),round(0.3*length(infant),digits = 0))
# trainsize和testsize保持一致，且为训练集或验证集中样本量最小集合中关注结局（此处为Disease）样本量的0.6

rfPermutedMultipleNull_child_adolescent <- iterative_rf_multiple_with_null_distribution(clr(dat_health[,bac_species_name]),
                                                                                        AgeCohortTrain,AgeCohortTest,child_adolescent,c(infant,young,middle,old),200,10,round(0.3*length(child_adolescent),digits = 0),round(0.3*length(child_adolescent),digits = 0))

rfPermutedMultipleNull_young <- iterative_rf_multiple_with_null_distribution(clr(dat_health[,bac_species_name]),
                                                                             AgeCohortTrain,AgeCohortTest,young,c(child_adolescent,infant,middle,old),200,10,round(0.3*length(young),digits = 0),round(0.3*length(young),digits = 0))

rfPermutedMultipleNull_middle <- iterative_rf_multiple_with_null_distribution(clr(dat_health[,bac_species_name]),
                                                                              AgeCohortTrain,AgeCohortTest,middle,c(child_adolescent,young,infant,old),200,10,round(0.3*length(middle),digits = 0),round(0.3*length(middle),digits = 0))

rfPermutedMultipleNull_old <- iterative_rf_multiple_with_null_distribution(clr(dat_health[,bac_species_name]),
                                                                           AgeCohortTrain,AgeCohortTest,old,c(child_adolescent,young,middle,infant),200,10,round(0.3*length(old),digits = 0),round(0.3*length(old),digits = 0))

# infant
rfTestAUCs_infant <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_infant) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_infant$AUC <- c(apply(rfPermutedMultipleNull_infant$AUCSameAge,1,median),apply(rfPermutedMultipleNull_infant$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_infant$AUCNull1,1,median))
rfTestAUCs_infant$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_infant$Group <- c(rep("infant",600))
rfTestAUCs_infant$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# child_adolescent
rfTestAUCs_child_adolescent <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_child_adolescent) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_child_adolescent$AUC <- c(apply(rfPermutedMultipleNull_child_adolescent$AUCSameAge,1,median),apply(rfPermutedMultipleNull_child_adolescent$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_child_adolescent$AUCNull1,1,median))
rfTestAUCs_child_adolescent$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_child_adolescent$Group <- c(rep("child_adolescent",600))
rfTestAUCs_child_adolescent$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# young
rfTestAUCs_young <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_young) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_young$AUC <- c(apply(rfPermutedMultipleNull_young$AUCSameAge,1,median),apply(rfPermutedMultipleNull_young$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_young$AUCNull1,1,median))
rfTestAUCs_young$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_young$Group <- c(rep("young",600))
rfTestAUCs_young$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# middle
rfTestAUCs_middle <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_middle) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_middle$AUC <- c(apply(rfPermutedMultipleNull_middle$AUCSameAge,1,median),apply(rfPermutedMultipleNull_middle$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_middle$AUCNull1,1,median))
rfTestAUCs_middle$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_middle$Group <- c(rep("middle",600))
rfTestAUCs_middle$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# old
rfTestAUCs_old <- as.data.frame(matrix(NA,600,4))
colnames(rfTestAUCs_old) <- c("AUC","ClassificationGroup","Group","Iter")
rfTestAUCs_old$AUC <- c(apply(rfPermutedMultipleNull_old$AUCSameAge,1,median),apply(rfPermutedMultipleNull_old$AUCDiffAge,1,median),apply(rfPermutedMultipleNull_old$AUCNull1,1,median))
rfTestAUCs_old$ClassificationGroup <- factor(c(rep("TrainSet",200),rep("TestSet",200),rep("MergedSet",200)),levels=c("TrainSet","TestSet","MergedSet"))
rfTestAUCs_old$Group <- c(rep("old",600))
rfTestAUCs_old$Iter <- c(paste0("Iter",c(1:200)),paste0("Iter",c(1:200)),paste0("Iter",c(1:200)))

# merge
rfTestAUCs_age <- rbind(rfTestAUCs_infant,rfTestAUCs_child_adolescent,rfTestAUCs_young,rfTestAUCs_middle,rfTestAUCs_old)
rfTestAUCs_age$Group <- factor(rfTestAUCs_age$Group,levels = c("infant","child_adolescent","young","middle","old",labels=c("Infant","Child&Adolescent","Young adult","Middle-aged adult","Elderly")))
write.csv(rfTestAUCs_age,"rfTestAUCs_age.csv",quote = F,row.names = F)

# mean and SD of AUC across groups
rfTestAUCs_age$type <- paste(rfTestAUCs_age$Group,rfTestAUCs_age$ClassificationGroup,sep = "_")
group_by(rfTestAUCs_age[,c("AUC","type")],type) %>% summarise_each(funs(mean,sd))

# plot Figure 3a---------------------------------------------------------------------
library(reshape2)
library(ggforce)
library(dplyr)
library(ggpubr)

p <- ggplot(rfTestAUCs_age)+
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
ggsave("boxplot_line_age.pdf", height = 6, width = 14)

# threshold determination Figure 3c--------------------------------------------------------

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

rfThreshold_infant <- rfThreshold(rfPermutedMultipleNull_infant,"infant")
rfThreshold_child_adolescent <- rfThreshold(rfPermutedMultipleNull_child_adolescent,"child_adolescent")
rfThreshold_young <- rfThreshold(rfPermutedMultipleNull_young,"young")
rfThreshold_middle <- rfThreshold(rfPermutedMultipleNull_middle,"middle")
rfThreshold_old <- rfThreshold(rfPermutedMultipleNull_old,"old")

rfThreshold_age <- Reduce(function(x,y) merge(x,y,by="threshold",all=T),list(rfThreshold_infant,rfThreshold_child_adolescent,rfThreshold_young,rfThreshold_middle,rfThreshold_old))
rfThreshold_age_melt <- melt(rfThreshold_age[,c(1,2,4,6,8,10)],id.vars = "threshold",value.name = "mean",variable.name = "age")

sample_select <- rownames(rfThreshold_age_melt) %in% rownames(rfThreshold_age_melt[rfThreshold_age_melt$age=="mean_featureScorechild_adolescent"|rfThreshold_age_melt$age=="mean_featureScoreold",])
p1 <- ggplot(data = rfThreshold_age_melt[sample_select,], aes(x = threshold, y = mean, color = age)) +
  geom_line(data = rfThreshold_age_melt[sample_select,], aes(x = threshold, y = mean),linewidth=1.5) +
  scale_color_manual(values = c("#00a087","#82491e"))+
  labs(x="Percentile Rank of Mean Decrease in Gini ",y="The Mean of Mean Decrease in Gini of Remaining Species")+
  scale_x_continuous(breaks = seq(0,1,0.1),
                     label = c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))+
  geom_vline(xintercept = 0.8,linetype="dashed",color="black")+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12))
ggsave("rfThreshold_age1.pdf",width = 6,height = 6)
p2 <- ggplot(data = rfThreshold_age_melt[!sample_select,], aes(x = threshold, y = mean, color = age)) +
  geom_line(data = rfThreshold_age_melt[!sample_select,], aes(x = threshold, y = mean),linewidth=1.5) +
  scale_color_manual(values = c("#cdbe6a","#86a4cf","#e64b35"))+
  labs(x="Percentile Rank of Mean Decrease in Gini ",y="The Mean of Mean Decrease in Gini of Remaining Species")+
  scale_x_continuous(breaks = seq(0,1,0.1),
                     label = c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))+
  geom_vline(xintercept = 0.8,linetype="dashed",color="black")+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(color = "black",size = 12))
ggsave("rfThreshold_age2.pdf",width = 6,height = 6)

# species marker

rfSpeciesMarker_infant<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_infant$featureProfile))>=0.8))
rfSpeciesMarker_child_adolescent<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_child_adolescent$featureProfile))>=0.8))
rfSpeciesMarker_young<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_young$featureProfile))>=0.8))
rfSpeciesMarker_middle<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_middle$featureProfile))>=0.8))
rfSpeciesMarker_old<-names(which(rank_scale(colMeans(rfPermutedMultipleNull_old$featureProfile))>=0.8))

rfSpeciesMarker_age <- data.frame(rfSpeciesMarker_infant,rfSpeciesMarker_child_adolescent,rfSpeciesMarker_young,rfSpeciesMarker_middle,rfSpeciesMarker_old)
write.csv(rfSpeciesMarker_age,"rfSpeciesMarker_age.csv",quote = F,row.names = F)

# Figure 3d; heatmap of age-specific species marker -----------------------
age_species_marker <- t(read.csv("age_species_marker.csv",header = T,row.names = 1))
library(pheatmap)
pdf("age_species_marker.pdf",height = 5,width = 15)
pheatmap(age_species_marker,scale = "none",cluster_rows = F,cluster_cols = F,
         fontsize = 8,border_color = "white",
         cellwidth = 12,cellheight = 12,
         color = colorRampPalette(c("#c2d7e7","#9abcd7","#679ac2","#3378ae","#00579a"))(500))
dev.off()

# species marker comparisons, Table S5-----------------------------------------------------

rfSpeciesMarker_age_unique <- unique(c(rfSpeciesMarker_infant,rfSpeciesMarker_child_adolescent,rfSpeciesMarker_young,rfSpeciesMarker_middle,rfSpeciesMarker_old))
rfFeatureProfile_age <- rbind(rfPermutedMultipleNull_infant$featureProfile,rfPermutedMultipleNull_child_adolescent$featureProfile,rfPermutedMultipleNull_young$featureProfile,rfPermutedMultipleNull_middle$featureProfile,rfPermutedMultipleNull_old$featureProfile)
rfMeanMarkerScores_age <- cbind(colMeans(rfFeatureProfile_age[1:200,]),colMeans(rfFeatureProfile_age[201:400,]),colMeans(rfFeatureProfile_age[401:600,]),colMeans(rfFeatureProfile_age[601:800,]),colMeans(rfFeatureProfile_age[801:1000,]))
colnames(rfMeanMarkerScores_age) <- c("MeanFeatureProfile_infant","MeanFeatureProfile_child_adolescent","MeanFeatureProfile_young","MeanFeatureProfile_middle","MeanFeatureProfile_old")
rfComparison_age <- batch_dunns(rfFeatureProfile_age[,rfSpeciesMarker_age_unique],as.factor(c(rep("infant",200),rep("child_adolescent",200),rep("young",200),rep("middle",200),rep("old",200))))
rfComparisonSummary_age <- cbind(rfMeanMarkerScores_age[rfSpeciesMarker_age_unique,],rfComparison_age)
rfComparisonSummary_age <- as.data.frame(rfComparisonSummary_age)
rfComparisonSummary_age$infant <- ifelse(rownames(rfComparisonSummary_age) %in% rfSpeciesMarker_infant,1,0)
rfComparisonSummary_age$child_adolescent <- ifelse(rownames(rfComparisonSummary_age) %in% rfSpeciesMarker_child_adolescent,1,0)
rfComparisonSummary_age$young <- ifelse(rownames(rfComparisonSummary_age) %in% rfSpeciesMarker_young,1,0)
rfComparisonSummary_age$middle <- ifelse(rownames(rfComparisonSummary_age) %in% rfSpeciesMarker_middle,1,0)
rfComparisonSummary_age$old <- ifelse(rownames(rfComparisonSummary_age) %in% rfSpeciesMarker_old,1,0)
write.csv(rfComparisonSummary_age,"rfComparisonSummary_age.csv",quote = F,row.names = T)

# Figure 3b, Figure S3; microbiota age-------------------------------------------- ----------------------------

library(reshape2)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- detectCores()
registerDoParallel(cl)

# all age groups; Figure 3b-----------------------------------------------------------------------------
dat_health_2 <- subset(dat_health,!is.na(age))
dat_health_2[,bac_species] <- clr(dat_health_2[,bac_species])

results <- as.data.frame(matrix(NA,nrow = nrow(dat_health_2),ncol = length(bac_species)))
rownames(results) <- rownames(dat_health_2)
colnames(results) <- bac_species_name
for (i in 1:length(bac_species_name)){
  model=glm(dat_health_2[,bac_species_name[i]]~as.factor(country)+as.factor(continent)+as.factor(sequencing_platform)+as.factor(DNA_extraction_kit)+non_westernized,data = dat_health_2,family = "gaussian")
  results[,i]=model$residuals
}

dat_health_2 <- cbind(results,dat_health_2$age)
colnames(dat_health_2)[dim(dat_health_2)[2]] <- "age"

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2, dat_health_2$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2[, c(optimal_num_features, ncol(dat_health_2))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2$age)
write.table(microbiota.age,"microbiotaAge_microbiota.txt",quote = F,row.names = F)

plot(microbiota.age$MicrobiotaAge,microbiota.age$ChronologicalAge)

pdf("ChrMicro_Age_SmoothScatter.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess.pdf",width = 6,height = 6)

# infant; Figure S3a-----------------------------------------------------------------------------
dat_health_2.new <- dat_health_2[infant,]

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2.new[, -ncol(dat_health_2.new)], dat_health_2.new$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2.new[, c(optimal_num_features, ncol(dat_health_2.new))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2.new)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2.new$age)
write.table(microbiota.age,"microbiotaAge_microbiota_infant.txt",quote = F,row.names = F)

pdf("ChrMicro_Age_SmoothScatter_infant.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  # geom_line(linewidth = 2,color="#4677b7",method)+
  geom_smooth(method = "lm", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess_infant.pdf",width = 6,height = 6)

# child_adolescent; Figure S3b-----------------------------------------------------------------------------
dat_health_2.new <- dat_health_2[child_adolescent,]

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2.new[, -ncol(dat_health_2.new)], dat_health_2.new$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2.new[, c(optimal_num_features, ncol(dat_health_2.new))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2.new)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2.new$age)
write.table(microbiota.age,"microbiotaAge_microbiota_child_adolescent.txt",quote = F,row.names = F)

pdf("ChrMicro_Age_SmoothScatter_child_adolescent.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  # geom_line(linewidth = 2,color="#4677b7",method)+
  geom_smooth(method = "loess", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess_child_adolescent.pdf",width = 6,height = 6)

# young; Figure S3c-----------------------------------------------------------------------------
dat_health_2.new <- dat_health_2[young,]

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2.new[, -ncol(dat_health_2.new)], dat_health_2.new$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2.new[, c(optimal_num_features, ncol(dat_health_2.new))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2.new)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2.new$age)
write.table(microbiota.age,"microbiotaAge_microbiota_young.txt",quote = F,row.names = F)

pdf("ChrMicro_Age_SmoothScatter_young.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  # geom_line(linewidth = 2,color="#4677b7",method)+
  geom_smooth(method = "loess", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess_young.pdf",width = 6,height = 6)

# middle; Figure S3d-----------------------------------------------------------------------------
dat_health_2.new <- dat_health_2[middle,]

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2.new[, -ncol(dat_health_2.new)], dat_health_2.new$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2.new[, c(optimal_num_features, ncol(dat_health_2.new))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2.new)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2.new$age)
write.table(microbiota.age,"microbiotaAge_microbiota_middle.txt",quote = F,row.names = F)

pdf("ChrMicro_Age_SmoothScatter_middle.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  # geom_line(linewidth = 2,color="#4677b7",method)+
  geom_smooth(method = "loess", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess_middle.pdf",width = 6,height = 6)

# old; Figure S3e-----------------------------------------------------------------------------
dat_health_2.new <- dat_health_2[old,]

rfcv_result <- foreach(
  n=1,.combine=cbind,.packages = 'randomForest') %dopar% {
    rfcv(dat_health_2.new[, -ncol(dat_health_2.new)], dat_health_2.new$age, cv.fold = 100)}
optimal_num_features <- rfcv_result$n.var[rfcv_result$error.cv == min(rfcv_result$error.cv)]

microbiota.rf <- randomForest(age ~ ., data = dat_health_2.new[, c(optimal_num_features, ncol(dat_health_2.new))], ntree=10000, nPerm=100, mtry=optimal_num_features/3, proximity=TRUE, importance=TRUE) 
microbiota.taxa <- microbiota.rf$importance
microbiota.pred <- predict(microbiota.rf, dat_health_2.new)
microbiota.age <- data.frame(sample = names(microbiota.pred) , MicrobiotaAge = microbiota.pred, ChronologicalAge = dat_health_2.new$age)
write.table(microbiota.age,"microbiotaAge_microbiota_old.txt",quote = F,row.names = F)

pdf("ChrMicro_Age_SmoothScatter_old.pdf",width = 6,height = 6)
with(microbiota.age, smoothScatter(ChronologicalAge,MicrobiotaAge, main =""))
dev.off()

p <- ggplot(microbiota.age, aes(x = ChronologicalAge, y = MicrobiotaAge))+
  #geom_point(alpha = 0.3)+
  # geom_line(linewidth = 2,color="#4677b7",method)+
  geom_smooth(method = "loess", linewidth = 2,color="#4677b7")+
  labs(x = "Chronological Age", y = "Microbiota Age")+
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.title = element_text(size = 20,color = "black"))+
  theme_classic()
ggsave("ChrMicro_Age_Loess_old.pdf",width = 6,height = 6)

