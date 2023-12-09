
# batch effect identification ------------------------------------------------------------
library(dplyr)
library(openxlsx)
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
cl <- detectCores()
registerDoParallel(cl)

library(patchwork)
library(ape) 
library(vegan) 
library(ggsci)
library(ggplot2) 
library(openxlsx)
library(ggpubr)

# sequencing_platform (as an example) -----------------------------------------------------------
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

# batch effect after adjusting -----------------------------------------------------------

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
cl <- detectCores()
registerDoParallel(cl)

# sequencing_platform (as an example) -----------------------------------------------------------
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

# proportion of microbial variation among healthy women------------------------------ ----------------------------------

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
cl <- detectCores()
registerDoParallel(cl)

# country (as ana example)-----------------------------------------------------------
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

# contribution of age and region to individual microbiota-------------------------- -----------------

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

# pie chart of phylum proportion across age groups and regions------------------------ ------------

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

# alpha diversity comparisons among group----------------------------------------------- ----------------------------------

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

# topological structure comparison of different age --------------------------------------

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


# species distribution by age ----------------------- ----------------------------------------

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

# heatmap of age-specific species marker -----------------------
age_species_marker <- t(read.csv("age_species_marker.csv",header = T,row.names = 1))
library(pheatmap)
pdf("age_species_marker.pdf",height = 5,width = 15)
pheatmap(age_species_marker,scale = "none",cluster_rows = F,cluster_cols = F,
         fontsize = 8,border_color = "white",
         cellwidth = 12,cellheight = 12,
         color = colorRampPalette(c("#c2d7e7","#9abcd7","#679ac2","#3378ae","#00579a"))(500))
dev.off()

# microbiota age-------------------------------------------- ----------------------------

library(reshape2)
library(foreach)
library(doParallel)
cl <- detectCores()
registerDoParallel(cl)

# all age groups -----------------------------------------------------------------------------
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
