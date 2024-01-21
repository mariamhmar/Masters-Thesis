library(tidyverse)
library(vegan)
library(indicspecies)
library(devtools)
library(lme4)
library(patchwork)
library(lmerTest)
##Import OTU counts table and Metadata
FullSRERDataset<-read.csv("FinalSRERdataOTUsandMetaData-ForAnalysis-27June2023.csv", header=TRUE, fileEncoding="UTF-8-BOM")
##Remove samples where richness=0 
filtered_FullData <- subset(FullSRERDataset, Richness != 0)
##Seperate data to generate an OTU table for distance matrix
##MetaDataTable
MetaSRER<-filtered_FullData[,1:12]
  column_to_rownames(var = "Sample")
  rownames(MetaSRER) <- MetaSRER$Sample
###Convert distance to a factor since it is being used as a metric to measure disturbance levels
MetaSRER$Distance <- as.factor(MetaSRER$Distance)
###OTU Table
OTU2<-filtered_FullData[,c(1,13:350)]
rownames(OTU2) <- OTU2$Sample
###Filtering samples with low read counts  
sumtable<-function(x){
  if(is.numeric(x)){
    sum(x) > 10
  } else {
    TRUE
  }
}
CH2OTUfilter<-OTU2[, sapply(OTU2,  sumtable)]
##nothing filtered

OTU2$Sample <- NULL
shannon<-diversity(OTU2)
  column_to_rownames(var="Sample")
###Generate Distance matrix and NMDSplot
dm.method <- 'bray'
dm <- vegdist(OTU2, method=dm.method)
otu.nmds <- metaMDS(dm,
                    k = 2,
                    maxit = 999,
                    trymax = 500,
                    wascores = TRUE)
##Check model
stressplot(otu.nmds)
###NMDS all Pastures
plot(otu.nmds)

nmds.scores <- as.data.frame(scores(otu.nmds, display = 'sites'))
nmds.scores$site <- rownames(nmds.scores) 
nmds.scores$Sample <- nmds.scores$site
nmds.scores$site <- NULL
nmds.scores <- left_join(nmds.scores, MetaSRER, by = 'Sample')
###Graphing all colored by Distance
nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(size=3, aes(color=Distance, shape=Pasture)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", shape="Pasture", y = "NMDS2") 
nmds_plot
###Graphing all colored by pasture
nmds_plot_PG <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(size=3, aes(color=Pasture, shape=Distance)) +
  stat_ellipse(aes(color=Pasture)) +
  labs(x = "NMDS1", colour = "Pasture", shape="Distance", y = "NMDS2") 
nmds_plot_PG


##Re-Running NMDS but individually for each pasture
Meta5s<-subset(MetaSRER, Pasture == "5S")
OTU5s<- rownames_to_column(OTU2, var = "Sample")
OTU5sonly<-merge(Meta5s, OTU5s, by.x="Sample", all=FALSE)
OTU5sonly<-OTU5sonly[,-c(2:12)] %>%
column_to_rownames(var="Sample")
###Distance matrix and stressplot
dm.method <- 'bray'
dm5s <- vegdist(OTU5sonly, method=dm.method)
otu5s.nmds <- metaMDS(dm5s,
                    k = 2,
                    maxit = 999,
                    trymax = 500,
                    wascores = TRUE)
stressplot(otu5s.nmds)
otu5s.nmds
nmds5s.scores <- as.data.frame(scores(otu5s.nmds, display = 'sites'))
nmds5s.scores$site <- rownames(nmds5s.scores) 
nmds5s.scores$Sample <- nmds5s.scores$site
nmds5s.scores$site <- NULL
nmds5s.scores <- left_join(nmds5s.scores, MetaSRER, by = 'Sample')
nmds_plot5s <- ggplot(nmds5s.scores,
                     aes(x = NMDS1,
                         y = NMDS2)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 5s") 
nmds_plot5s

####Visualizing Pasture 5S richness
ggplot(Meta5s, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
##Getting rid of dataset outliers based on richness values
list_quantiles5s <- tapply(Meta5s$Richness, Meta5s$Distance, quantile)
Q1s5s <- sapply(1:3, function(i) list_quantiles5s[[i]][2])
Q3s5s <- sapply(1:3, function(i) list_quantiles5s[[i]][4])
IQRs5s <- tapply(Meta5s$Richness, Meta5s$Distance, IQR)
Lowers5s <- Q1s5s - 1.5*IQRs5s
Uppers5s <- Q3s5s + 1.5*IQRs5s
datas5s <- split(Meta5s, Meta5s$Distance)
##Re-Visualizing Richness 5s no outliers
data_no_outlier5s <- NULL
for (i in 1:3){
  out5s <- subset(datas5s[[i]], datas5s[[i]]$Richness > Lowers5s[i] & datas5s[[i]]$Richness < Uppers5s[i])
  data_no_outlier5s <- rbind(data_no_outlier5s, out5s)
}

dim(data_no_outlier5s)
###Graph it
ggplot(data_no_outlier5s, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
###Re-Trying NMDS with no_outliers dataset
##Get rid of outliers from OTU table
different_rows5s <- anti_join(Meta5s, data_no_outlier5s, by = "Sample")
print(rownames(different_rows5s))
delete5s<-c("BSC5562_004", "DOW5562_031", "BSC5562_016", "BSC5562_018")
otu5s_nooutliers<-OTU5sonly[!(rownames(OTU5sonly) %in% delete5s), ]

dm.method <- 'bray'
dmnoout5s <- vegdist(otu5s_nooutliers, method=dm.method)
otunoout5s.nmds <- metaMDS(dmnoout5s,
                      k = 2,
                      maxit = 999,
                      trymax = 500,
                      wascores = TRUE)
stressplot(otunoout5s.nmds)
otunoout5s.nmds
nmds5snoout.scores <- as.data.frame(scores(otunoout5s.nmds, display = 'sites'))
nmds5snoout.scores$site <- rownames(nmds5snoout.scores) 
nmds5snoout.scores$Sample <- nmds5snoout.scores$site
nmds5snoout.scores$site <- NULL
nmds5snoout.scores <- left_join(nmds5snoout.scores, MetaSRER, by = 'Sample')
nmds_plot5snoout <- ggplot(nmds5snoout.scores,
                      aes(x = NMDS1,
                          y = NMDS2)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 5S") + theme(panel.grid.major = element_blank(), 
                                                                                  panel.grid.minor = element_blank())
nmds_plot5snoout
##5s no outliers check statistics if significant differences by distance
Meta5s_nooutliers<-Meta5s[!(rownames(Meta5s) %in% delete5s), ]
adonis2(vegdist(otu5s_nooutliers, method='bray')~Distance, data=Meta5s_nooutliers, permutations = 9999)

###Pasture 5N NMDS
Meta5N<-subset(MetaSRER, Pasture == "5N")
OTU5N<- rownames_to_column(OTU2, var = "Sample")
OTU5Nonly<-merge(Meta5N, OTU5N, by.x="Sample", all=FALSE)
OTU5Nonly<-OTU5Nonly[,-c(2:12)] %>%
  column_to_rownames(var="Sample")
dm.method <- 'bray'
dm5N <- vegdist(OTU5Nonly, method=dm.method)
otu5N.nmds <- metaMDS(dm5N,
                      k = 2,
                      maxit = 999,
                      trymax = 500,
                      wascores = TRUE)
stressplot(otu5N.nmds)
nmds5N.scores <- as.data.frame(scores(otu5N.nmds, display = 'sites'))
nmds5N.scores$site <- rownames(nmds5N.scores) 
nmds5N.scores$Sample <- nmds5N.scores$site
nmds5N.scores$site <- NULL
nmds5N.scores <- left_join(nmds5N.scores, MetaSRER, by = 'Sample')
nmds_plot5N <- ggplot(nmds5N.scores,
                      aes(x = NMDS1,
                          y = NMDS2)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 5N") 
nmds_plot5N
###Richness 5n
ggplot(Meta5N, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
##Attempt to get rid of outliers by identifying upper and lower IQR for 5N pasture 
quartiles5N<- quantile(Meta5N$Richness, probs=c(.25, .75), na.rm= FALSE)
IQR5N<- IQR(Meta5N$Richness)
Lower5N<- quartiles5N[1]-1.5*IQR5N
upper5N<- quartiles5N[2]+1.5*IQR5N
Richness5Nnooutliers<- subset(Meta5N, Meta5N$Richness> Lower5N & Meta5N$Richness <upper5N)
##Plotting with no outliers
ggplot(Richness5Nnooutliers, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
###Try getting rid of outliers by richness and distance
list_quantiles5n <- tapply(Meta5N$Richness, Meta5N$Distance, quantile)
Q1s5N <- sapply(1:3, function(i) list_quantiles5n[[i]][2])
Q3s5N <- sapply(1:3, function(i) list_quantiles5n[[i]][4])
IQRs5ng <- tapply(Meta5N$Richness, Meta5N$Distance, IQR)
Lowers5ng <- Q1s5N - 1.5*IQRs5ng
Uppers5ng <- Q3s5N + 1.5*IQRs5ng
datas5n <- split(Meta5N, Meta5N$Distance)

data_no_outlier5n <- NULL
for (i in 1:3){
  out5n <- subset(datas5n[[i]], datas5n[[i]]$Richness > Lowers5ng[i] & datas5n[[i]]$Richness < Uppers5ng[i])
  data_no_outlier5n <- rbind(data_no_outlier5n, out5n)
}

dim(data_no_outlier5n)
###Graph it
ggplot(data_no_outlier5n, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
##Trying the NMDS after removing outliers 5N
different_rows5n <- anti_join(Meta5N, data_no_outlier5n, by = "Sample")
print(rownames(different_rows5n))
delete5n<-c("BSC5562_028", "BSC5562_031", "DOW5562_048")
otu5n_nooutliers<-OTU5Nonly[!(rownames(OTU5Nonly) %in% delete5n), ]

dm.method <- 'bray'
dmnoout5n <- vegdist(otu5n_nooutliers, method=dm.method)
otunoout5n.nmds <- metaMDS(dmnoout5n,
                           k = 2,
                           maxit = 999,
                           trymax = 500,
                           wascores = TRUE)
stressplot(otunoout5n.nmds)
otunoout5n.nmds
nmds5nnoout.scores <- as.data.frame(scores(otunoout5n.nmds, display = 'sites'))
nmds5nnoout.scores$site <- rownames(nmds5nnoout.scores) 
nmds5nnoout.scores$Sample <- nmds5nnoout.scores$site
nmds5nnoout.scores$site <- NULL
nmds5nnoout.scores <- left_join(nmds5nnoout.scores, MetaSRER, by = 'Sample')
nmds_plot5nnoout <- ggplot(nmds5nnoout.scores,
                           aes(x = NMDS1,
                               y = NMDS2)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 5N") + theme(panel.grid.major = element_blank(), 
                                                                                  panel.grid.minor = element_blank())
nmds_plot5nnoout


###6B Only
Meta6B<-subset(MetaSRER, Pasture == "6B")
OTU6B<- rownames_to_column(OTU2, var = "Sample")
OTU6Bonly<-merge(Meta6B, OTU6B, by.x="Sample", all=FALSE)
OTU6Bonly<-OTU6Bonly[,-c(2:12)] %>%
  column_to_rownames(var="Sample")
dm.method <- 'bray'
dm6B <- vegdist(OTU6Bonly, method=dm.method)
otu6B.nmds <- metaMDS(dm6B,
                      k = 2,
                      maxit = 999,
                      trymax = 500,
                      wascores = TRUE)
stressplot(otu6B.nmds)
nmds6B.scores <- as.data.frame(scores(otu6B.nmds, display = 'sites'))
nmds6B.scores$site <- rownames(nmds6B.scores) 
nmds6B.scores$Sample <- nmds6B.scores$site
nmds6B.scores$site <- NULL
nmds6B.scores <- left_join(nmds6B.scores, MetaSRER, by = 'Sample')
nmds_plot6B <- ggplot(nmds6B.scores,
                      aes(x = NMDS1,
                          y = NMDS2,)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance))+
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 6B") 
nmds_plot6B
ordihull(nmds6B.scores,groups=Distance,draw="polygon",col="grey90",label=F)
Most_Disturbance6B<- nmds.scores[nmds.scores$Distance == "0", ][chull(nmds6B.scores[nmds.scores$Distance == 
                                                                                  "0", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
Middle_Disturbance6B<- nmds.scores[nmds6B.scores$Distance == "500", ][chull(nmds.scores[nmds.scores$Distance == 
                                                                                      "500", c("NMDS1", "NMDS2")]), ]  
Least_Disturbance6B<- nmds.scores[nmds6B.scores$Distance == "1000", ][chull(nmds.scores[nmds.scores$Distance == 
                                                                                      "1000", c("NMDS1", "NMDS2")]), ]
hull.data6B <- rbind(Most_Disturbance6B, Middle_Disturbance6B, Least_Disturbance6B)
###Richness6B
ggplot(Meta6B, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))
##getting rid of outliers
list_quantiles6b <- tapply(Meta6B$Richness, Meta6B$Distance, quantile)
Q1s6b <- sapply(1:3, function(i) list_quantiles6b[[i]][2])
Q3s6b <- sapply(1:3, function(i) list_quantiles6b[[i]][4])
IQRs6b <- tapply(Meta6B$Richness, Meta6B$Distance, IQR)
Lowers6b <- Q1s6b - 1.5*IQRs6b
Uppers6b <- Q3s6b + 1.5*IQRs6b
datas6b <- split(Meta6B, Meta6B$Distance)

data_no_outlier6b <- NULL
for (i in 1:3){
  out6b <- subset(datas6b[[i]], datas6b[[i]]$Richness > Lowers6b[i] & datas6b[[i]]$Richness < Uppers6b[i])
  data_no_outlier6b <- rbind(data_no_outlier6b, out6b)
}

dim(data_no_outlier6b)
###Graph it
ggplot(data_no_outlier6b, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance") + 
  ylab("Ln Richness") + theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))

###NMDS for 6B
different_rows6B <- anti_join(Meta6B, data_no_outlier6b, by = "Sample")
print(rownames(different_rows6B))
delete6b<-c("BSC5562_052", "BSC5562_064", "DOW5562_079", "DOW5562_085", "BSC5562_069", "DOW5562_087")
otu6b_nooutliers<-OTU6Bonly[!(rownames(OTU6Bonly) %in% delete6b), ]

dm.method <- 'bray'
dmnoout6b <- vegdist(otu6b_nooutliers, method=dm.method)
otunoout6b.nmds <- metaMDS(dmnoout6b,
                           k = 2,
                           maxit = 999,
                           trymax = 500,
                           wascores = TRUE)
stressplot(otunoout6b.nmds)
otunoout6b.nmds
nmds6bnoout.scores <- as.data.frame(scores(otunoout6b.nmds, display = 'sites'))
nmds6bnoout.scores$site <- rownames(nmds6bnoout.scores) 
nmds6bnoout.scores$Sample <- nmds6bnoout.scores$site
nmds6bnoout.scores$site <- NULL
nmds6bnoout.scores <- left_join(nmds6bnoout.scores, MetaSRER, by = 'Sample')
nmds_plot6bnoout <- ggplot(nmds6bnoout.scores,
                           aes(x = NMDS1,
                               y = NMDS2)) +
  geom_point(size=3, aes(color=Distance)) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="Pasture 6B") + theme(panel.grid.major = element_blank(), 
                                                                                  panel.grid.minor = element_blank())
nmds_plot6bnoout


###NMDS all no outliers in one plot
deleteall<-c("BSC5562_004", "DOW5562_031", "BSC5562_016", "BSC5562_018","BSC5562_028", "BSC5562_031", "DOW5562_048","BSC5562_052", "BSC5562_064", "DOW5562_079", "DOW5562_085", "BSC5562_069", "DOW5562_087")
####Meta no outliers

Meta_no_outliers<-MetaSRER[!(rownames(MetaSRER) %in% deleteall), ]
richnesspast<-ggplot(Meta_no_outliers, aes(y = Richness, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance in Meters") + 
  ylab("Richness") + ggtitle(label="Species Richness by Distance and Pasture")+ theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))+facet_wrap(~Pasture)+theme(panel.grid.major = element_blank(), 
                                                                                                             panel.grid.minor = element_blank())
richnesspast
dm.method <- 'bray'
dm <- vegdist(all_no_outliers, method=dm.method)
otunoout.nmds <- metaMDS(dm,
                         k = 2,
                         maxit = 999,
                         trymax = 500,
                         wascores = TRUE)
otunoout.nmds

nmdsnoout.scores <- as.data.frame(scores(otunoout.nmds, display = 'sites'))
nmdsnoout.scores$site <- rownames(nmdsnoout.scores) 
nmdsnoout.scores$Sample <- nmdsnoout.scores$site
nmdsnoout.scores$site <- NULL
nmdsnoout.scores <- left_join(nmdsnoout.scores, Meta_no_outliers, by = 'Sample')
nmds_noout_all <- ggplot(nmdsnoout.scores,
                         aes(x = NMDS1,
                             y = NMDS2)) +
  geom_point(size=3, aes(color=Distance, shape=Pasture)) +theme(panel.grid.major = element_blank(), 
                                                                panel.grid.minor = element_blank()) +
  stat_ellipse(aes(color=Distance)) +
  labs(x = "NMDS1", colour = "Distance", y = "NMDS2", title="NMDS All Pastures") 

all_no_outliers<-OTU2[!(rownames(OTU2) %in% deleteall), ]
dm.method <- 'bray'
dm2 <- vegdist(all_no_outliers, method=dm.method)
otu2.nmds <- metaMDS(dm2,
                     k = 2,
                     maxit = 999,
                     trymax = 500,
                     wascores = TRUE)
stressplot(otu2.nmds)
##Graph
nmds_noout_all

###Are there significant differences in composition by distance
adonis2(vegdist(all_no_outliers, method="bray") ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers, permutations=9999)

###Repeating by trying to filter data set further
##Fulldataset Significance
RichnessandPastureanova <- adonis2(all_no_outliers ~ Distance + Pasture, data = Meta_no_outliers, permutations = 999)
print(RichnessandPastureanova)


##Export datasets
write_excel_csv(all_no_outliers, "Otutablenooutliers.xls")
write_excel_csv(Meta_no_outliers, "Metadatanooutliers.xls")

##Indic species
##Add all no outliers to site Info
##group by pasture and transect 
###Full dataset
spec_for_indic<-read.csv("All_no_outliers.ind.csv", header=T)
#indspec_all <- multipatt(spec_for_indic, groups, 
                   # control = how(nperm=999)) 
#write.csv(Meta_no_outliers, "meta_No_outliers.csv")
groups<-read.csv("Groupsforindicator.csv",header=T)
write.csv(all_no_outliers, "All_no_outliers.csv")
spec_for_indic<-read.csv("All_no_outliers.ind.csv", header=T)
abund<- spec_for_indic[,3:ncol(spec_for_indic)]
group<-spec_for_indic$Pasture.Distance
indspec_all <- multipatt(abund, group, func = "IndVal.g",
                         control = how(nperm=999)) 
summary(indspec_all)
indspec_all

### By Pasture
Ind.pst<-read.csv("All_no_outliers.csv", header = T)
abund2<- Ind.pst[,3:ncol(spec_for_indic)]
group2<-Ind.pst$Pasture
indspec_past <- multipatt(abund2, group2, func = "IndVal.g", duleg = TRUE, 
                         control = how(nperm=999)) 
summary(indspec_past)
###By Distance
spec_for_indic<-read.csv("All_no_outliers.ind2.csv", header=T)
abund2<- spec_for_indic[,3:ncol(spec_for_indic)]
group2<-spec_for_indic$Distance
indspec_all2 <- multipatt(abund2, group2, func = "IndVal.g",
                         control = how(nperm=999)) 
summary(indspec_all2)


###Shannon_Diversity
Shannon_div_nooutliers<-diversity(all_no_outliers, index="shannon")
Shannon_div <- as.data.frame(scores(Shannon_div_nooutliers, display = 'sites'))
Shannon_div$site <- rownames(Shannon_div) 
Shannon_div$Sample <- Shannon_div$site
Shannon_div$site <- NULL
Meta_with_shannon <- left_join(Shannon_div, Meta_no_outliers, by = 'Sample')
Shannondivplot<-ggplot(Meta_with_shannon, aes(y = Dim1, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance in Meters") + 
  ylab("Shannon Diversity") + ggtitle(label="Shannon Diversity by Distance and Pasture")+ theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))+facet_wrap(~Pasture)+ theme(panel.grid.major = element_blank(),
                                                                                                            panel.grid.minor = element_blank())
Shannondivplot
##Evenness plot
data_evenness <- diversity(all_no_outliers) / log(specnumber(all_no_outliers)) 
evenness <- as.data.frame(scores(data_evenness, display = 'sites'))
evenness$site <- rownames(evenness) 
evenness$Sample <- evenness$site
evenness$site <- NULL
Meta_with_even<- left_join(evenness, Meta_no_outliers, by = 'Sample')
evennessplot<-ggplot(Meta_with_even, aes(y = Dim1, x = Distance, fill = Distance, group= Distance)) + geom_boxplot(outlier.shape = NA, alpha=0.8)+
  geom_jitter(aes(color = Distance))+ 
  xlab("Distance in Meters") + 
  ylab("Evenness") + ggtitle(label="Evenness by Distance and Pasture")+ theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90))+facet_wrap(~Pasture)
evennessplot
###Making a list of OTUs used in analysis to pull sequences from for TBAS identification 
tocopyotu<-t(all_no_outliers)
write.csv(tocopyotu, "Otususednooutlier.csv")
OTUList<-read.csv("OTUListonly.csv", header=TRUE)
View(OTUList)
OTUSequencemaster<-read.csv("OTUs_and_seq.csv")
Sequences_used<-merge(OTUList, OTUSequencemaster, by='ID', all=FALSE)
write.csv(Sequences_used, "SRERforTBAS.csv")

##ln Richness to normalize data for statistical analysis
# Assuming your dataset is already loaded as "Meta_no_outliers"
# Take the natural logarithm of the "Richness" column and store it in a new column
Meta_no_outliers$ln_Richness<- log(Meta_no_outliers$Richness)
MetaSRER$ln_Richness<-log(MetaSRER$Richness)
# View the updated dataset
head(Meta_no_outliers)
##Stats
Richnessstat<-aov(Meta_no_outliers$ln_Richness ~ Distance+ Error(Pasture), data=Meta_no_outliers)
summary(Richnessstat)
###Patching Figures together 
patchwork + plot_annotation(tag_levels = 'A')
patchwork<-(nmds_noout_all)/(nmds_plot5nnoout |nmds_plot5snoout | nmds_plot6bnoout )
Richnessinfo<-(richnesspast)/(Shannondivplot)
Richnessinfo + plot_annotation(tag_levels = 'A')
###Statistics
adonis2(Meta_no_outliers$ln_Richness ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers)
adonis2(Meta_no_outliers$ln_Richness ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers)
###Retrying above as ANOVA not PERMANOVA

model = lmer(Meta_no_outliers$ln_Richness ~ Distance +(1|Pasture), data=Meta_no_outliers,
             REML=TRUE)
anova(model)
rand(model)
  # Example: One-way ANOVA
anova_model_Richness <- aov(Meta_no_outliers$ln_Richness ~ Distance +Error(Pasture), data=Meta_no_outliers)
summary(anova_model_Richness)
Meta_no_outliers<- mutate(Meta_no_outliers, Pasture=as.factor(Pasture), ordered(Pasture))
anova_model_Richness
adonis2(MetaSRER$ln_Richness ~ Distance, strata=MetaSRER$Pasture, data=MetaSRER, permutations=9999)
adonis2(MetaSRER$ln_Richness ~ Distance, strata=MetaSRER$Pasture, data=MetaSRER)
adonis2(diversity(all_no_outliers, index="shannon") ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers, permutations=9999)
adonis2(diversity(all_no_outliers, index="shannon") ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers)

adonis2(diversity(OTU2, index="shannon") ~ Distance, strata=MetaSRER$Pasture, data=MetaSRER, permutations=9999)
adonis2(vegdist(all_no_outliers) ~ Distance, strata=Meta_no_outliers$Pasture, data=Meta_no_outliers, permutations=9999)


###Making tables of OTUs based on their distance 
Zero_Dist<-subset(Meta_no_outliers, Distance == "0")
Zero_dist_OTUS<-Zero_Dist[,-c(2:14)]
rownames(Zero_dist_OTUS) <- NULL
# Convert the matrix to a data frame
Zero_dist_OTUS <- as.data.frame(Zero_dist_OTUS)
# Apply column_to_rownames
Zero_dist_OTUS <- column_to_rownames(Zero_dist_OTUS, var = "Row.names")
Highest_Dist<-subset(Meta_no_outliers, Distance=="1000")
Highest_Dist<-Highest_Dist[,-c(2:14)]
rownames(Highest_Dist) <- NULL
Highest_Dist_OTUS <- as.data.frame(Highest_Dist)
# Apply column_to_rownames
Highest_Dist_OTUS <- column_to_rownames(Highest_Dist_OTUS, var = "Row.names")

# Calculate column sums
col_sums <- colSums(Highest_Dist_OTUS)

# Get the column names where sums are greater than 0
selected_column_names_1000 <- names(col_sums[col_sums > 0])

# Calculate column sums
col_sums_0 <- colSums(Zero_dist_OTUS )

# Get the column names where sums are greater than 0
selected_column_names_0 <- names(col_sums_0[col_sums_0 > 0])
maxlength = max(length(selected_column_names_0), length(selected_column_names_1000))
Otus_for_Tbas_2<- data.frame(
  Distance_0 = c(selected_column_names_0, rep(NA, maxlength - length(selected_column_names_0))),
  Distance_1000 = c(selected_column_names_1000, rep(NA, maxlength -length(selected_column_names_1000)))
)
OTU_TABLE_1000_OTUS<-subset(Highest_Dist_OTUS, select=c(selected_column_names_1000))
OTU_TABLE_0_OTUS<-subset(Zero_dist_OTUS, select=c(selected_column_names_0))
###Convert to presense absense

write.csv(OTU_TABLE_0_OTUS, "OTU_TABLE_0.csv")
OTU_TABLE_0_OTUS<-read.csv("OTU_TABLE_0.csv", header=T)
write.csv(OTU_TABLE_1000_OTUS, "OTU_TABLE_1000.csv")
OTU_TABLE_1000_OTUS<-read.csv("OTU_TABLE_1000.csv", header=T)
library(relabund)
install.packages("relabund")
PA_Table_1000<- ifelse(OTU_TABLE_1000_OTUS != 0, 1, 0)
PA_Table_0<- ifelse(OTU_TABLE_0_OTUS != 0, 1, 0)
library(tidyverse)
sums_1000<-colSums(PA_Table_1000)
sums_1000<- as.data.frame(sums_1000)
colnames(sums_1000) <- "Sums"
sums_1000<-rownames_to_column(sums_1000, var = "Query.sequence")

sums_0<-colSums(PA_Table_0)
sums_0<- as.data.frame(sums_0)
colnames(sums_0) <- "Sums"
sums_0<-rownames_to_column(sums_0, var = "Query.sequence")

Tax_Sums_1000<-merge(sums_1000, All_TBAS_TAX, by='Query.sequence', all=FALSE)

Tax_Sums_0<-merge(sums_0, All_TBAS_TAX, by='Query.sequence', all=FALSE)
###Aggregating by class and making a pie chart

aggregated_1000 <- Tax_Sums_1000 %>%
  group_by(Most.common.class.level.assignment) 
aggregated_1000<- aggregated_1000%>%
  summarize(SumSamples = sum(Sums), by=)
aggregated_0<- Tax_Sums_0%>%
  group_by(Most.common.class.level.assignment) %>%
  summarize(SumSamples = sum(Sums))


scale_fill_manual(values = c("red","#FF9900", "yellow","#CC9900", "#CC9999", "#CC99FF", "#99FF33", "#33FFFF", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#66CC33", "#FF99ff"))    
aggregated_1000 <- aggregated_1000[-c(4, 5, 6, 8, 9, 10), ]                           
Pie_1000<- ggplot(aggregated_1000, aes(x="", y=SumSamples, fill=Most.common.class.level.assignment))+
  labs(title = "Most Disturbance 1000M", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values = c("red","#FF9900", "yellow","#FF6600", "#CC9999", "#CC99FF", "#99FF33", "#33FFFF","#FFCCCC", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#CCCCCC", "#003300", "#FF99ff")) +
  theme_minimal()
Pie_1000
aggregated_0 <- aggregated_0[-c(4, 5, 6, 8, 9), ] 
Pie_0<- ggplot(aggregated_0, aes(x="", y=SumSamples, fill=Most.common.class.level.assignment))+ 
  labs(title = "Least Disturbance 0M", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values = c("red","#FF9900", "yellow","#FF6600", "#CC9999", "#CC99FF", "#99FF33", "#33FFFF", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#003300", "#FF99ff"))+ 
  theme_minimal()
Pie_0
write.csv(aggregated_0, "agg_0.csv")
write.csv(aggregated_1000, "agg_1000.csv")
###Putting into one chart
agg_0_1000<-read.csv("agg_1000_0.csv", header=T)
Pie_0_1000<- ggplot(agg_0_1000, aes(x="", y=Proportion, fill=Class))+ facet_wrap(~Distance)+
  labs(title = "Sample Distribution of Ascomycota by Class", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = .5, stat = "identity")+ 
  scale_fill_manual(values = c("red","#FF9900", "yellow","#FF6600", "#CC9999", "#99FF33", "#33FFFF", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#003300", "#FF99ff"))+ 
  theme_minimal()+theme(aspect.ratio = 2/1, panel.spacing = unit(-.08,"cm"))
Pie_0_1000
unitecombo<-read.csv("unite_1000_0.csv", header=T)
Unite_com<- ggplot(unitecombo, aes(x="", y=Proportion, fill=Class))+ facet_wrap(~Distance)+
  labs(title = "Sample Distribution of Basidiomycota by Class", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = .5, stat = "identity")+ 
  scale_fill_manual(values = c("#FF9900", "yellow", "#CC9999", "#99FF33", "#33FFFF", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#003300", "#FF99ff"))+ 
  theme_minimal()+theme(aspect.ratio = 2/1, panel.spacing = unit(-.08,"cm"))
Unite_com
Class +plot_annotation(
  title = "Ascomycota Sample Distribution by Class",
  theme = theme(plot.title = element_text(hjust = 0.5),)
)
Class
# Arrange the plots and title using patchwork
combined_plot <- (plot1 | plot2) / title
##Trying with Unite results to see if I get the same results
Tax_UNITE<-read.csv("SRERIDS_V2.csv", header=T)
Tax_UNITE_1000<-merge(sums_1000, Tax_UNITE, by='Query.sequence', all=FALSE)

Tax_UNITE_0<-merge(sums_0, Tax_UNITE, by='Query.sequence', all=FALSE)
###Aggregating by class and making a pie chart

aggregated_1000_Unite <- Tax_UNITE_1000 %>%
  group_by(Class) %>%
  summarize(SumSamples = sum(Sums))
aggregated_0_Unite<- Tax_UNITE_0%>%
  group_by(Class) %>%
  summarize(SumSamples = sum(Sums))
##Graphing Unite Results
write.csv(aggregated_1000_Unite, "unite_1000.csv")
write.csv(aggregated_0_Unite, "unite_0.csv")
aggregated_1000_Unite <- aggregated_1000_Unite[-12, ]
Pie_1000_Unite<- ggplot(aggregated_1000_Unite, aes(x="", y=SumSamples, fill=Class))+
  labs(title = "Basidiomycota Sample Distibution by Class 1000M", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values = c("red","#FF9900", "yellow","#FF6600", "#CC9999", "#99FF33", "#33FFFF","#FFCCCC", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#CCCCCC", "#003300", "#FF99ff")) +
  theme_minimal()
Pie_1000_Unite
aggregated_0_Unite <- aggregated_0_Unite[-13, ]
Pie_0_Unite<- ggplot(aggregated_0_Unite, aes(x="", y=SumSamples, fill=Class))+ 
  labs(title = "Basidiomycota Sample Distibution by Class 0M", groups="Taxonomic Class", fill="Taxonomic class", y="", x="") +
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values = c("red","#FF9900", "yellow","#FF6600", "#CC9999", "#99FF33", "#33FFFF","#FFCCCC", "#CC00FF", "#3399FF", "#660066", "#CCC000", "#CCCCCC", "#003300", "#FF99ff"))+ 
  theme_minimal()
Pie_0_Unite
Basidio<-(Pie_0_Unite)|(Pie_1000_Unite)
Basidio 
# Create a color palette
extended_palette <- rep(colors, length.out = 16)

# Create the pie chart with labels and a title
pie(
  aggregated_1000$SumSamples,
  labels = aggregated_1000$Most.common.class.level.assignment,
  col = extended_palette,
  main = "Sample Distribution by Class Level 1000 Distance",
  cex = 0.7  # Adjust label size
)
install.packages("viridis")
library(viridis)
# Create a color palette with 16 distinct colors from the 'viridis' palette
my_palette <- rainbow(16)



##Importing full assignments for merge
All_TBAS_TAX<-read.csv("assignments_TBAS_allSRER.csv", header=TRUE)


write.csv(Otus_for_Tbas_2, "OTUs_by_distance.csv")
###Merge with the Sequences used datafile 
Zero_OTU_only <-Otus_for_Tbas_2[,1] 
write.csv(Otus_for_Tbas_2, "OTUs_by_distance.csv")
thousand<-read.csv("1000_Dist_otus.csv", header=TRUE)
Zero<-read.csv("Zero_Dist_otus.csv")
ZeroForTBAS<-merge(Sequences_used, Zero, by='ID', all=FALSE)
ThousandForTBAS<-merge(Sequences_used, thousand, by='ID', all=FALSE)
write.csv(ThousandForTBAS, 'ThousandforTbas.csv')
write.csv(ZeroForTBAS, 'ZeroforTbas.csv')
#### Guilds for Indic Spec
All_Guilds<-read.csv("otus_tableL3H7L5KJ.guilds_all.csv", header=TRUE)
Zero_spec<-read.csv("Zero_Spec_indic_otus.csv", header=TRUE,  fileEncoding="UTF-8-BOM")
Zero_indic_guilds<-merge(Zero_spec, All_Guilds, by='strain', all=FALSE)
Thousand_Indic<-read.csv("thousand_Spec_indic_otus.csv", header=TRUE,  fileEncoding="UTF-8-BOM")
thousand_indic_guilds<-merge(Thousand_Indic, All_Guilds, by='strain', all=FALSE)
write.csv(thousand_indic_guilds, "thousand_indic.csv")
write.csv(Zero_indic_guilds, "zero_guilds_indic.csv")
###Combining Guild Info and Ta
##Indic Spec plus 
write.csv(Tax_Sums_0, "Zero_Dis_Weighted.csv")
write.csv(Tax_Sums_1000, "Thousand_Dis_weighted.csv")
zero_comb<-read.csv("0_comb.csv", header=T)
Indic_Zero_counts<-merge(zero_comb, Tax_Sums_0, by="Query.sequence", all=FALSE)
thousand_comb<-read.csv("1000_comb.csv", header=T)
Indic_1000_counts<-merge(thousand_comb, Tax_Sums_1000, by="Query.sequence", all=FALSE)
write.csv(Indic_1000_counts, "Indic_1000_weighted.csv")
write.csv(Indic_Zero_counts, "Indic_0_weighted.csv")
####Comparing to the UNITE IDS
comparing_unite<-read.csv("UniteIDs_For_Comparison.csv", header=TRUE, fileEncoding="UTF-8-BOM")
Unite_TBAS_1000<-merge(comparing_unite, Tax_Sums_1000, by="Query.sequence", all=FALSE)
Unite_TBAS_0<-merge(comparing_unite, Tax_Sums_0, by="Query.sequence", all=FALSE)
write.csv(Unite_TBAS_0, "Zero_Unite_TBAS_combinedoutputs.csv")
write.csv(Unite_TBAS_1000, "1000_Unite_TBAS_combinedoutputs.csv")
###Making Pie Chart
Pie_Data<-read.csv("Weighted_Pie_Chart_Data.csv", header=T)
pie<-ggplot(data=Pie_Data, aes(x=" ", y=count, group=Guild, colour=Guild, fill=Guild)) + 
  coord_polar("y")  +theme_void() 
library(ggplot2)

Pie_Data <- read.csv("Weighted_Pie_Chart_Data.csv", header = TRUE)
##To put labels outside of the plot
library(ggrepel)
# Create the pie chart
##get positions for labels

#Pie_Data$label_pos <- cumsum(Pie_Data$Percentage) - 0.5 * Pie_Data$Percentage

pie_chart <- ggplot(data = Pie_Data2, aes(x = "", y = Percentage, fill = Guild, label = Labels)) +
  geom_bar(stat = "identity") +
  coord_polar("y") +
  theme_void() + geom_col(color = "black") +
  guides(fill = guide_legend(title = "Guild")) +
  facet_wrap(~ Distance)+
  theme(text = element_text(size = 16, face = "bold"),  
    plot.title = element_text(size = 26, face = "bold")) + geom_text(aes(label = Labels,
                                                                         x = 1.2), position = position_stack(vjust = 0.5), size=5) +
   labs(fill = "Guild")

pie_chart
Propofclass<-read.csv("Proportionsofclasses.csv", header=T)
ggplot(Propofclass, aes(x = Distance, y = Proportion, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#07AA81","#F2872D", "#6AE458","#41EF89", "#DD00A7", "#26938E", "#AD7811", "#1C8924","#C83453","#9C07AA","yellow","red" ))+
  labs(title = "Fungal Class Distribution by Disturbance level",face = "bold", x = NULL, y = NULL) +
  theme_minimal()+theme(axis.text = element_text(size = 12, color="black"),
                        plot.title = element_text(hjust = 0.5))


weights<-read.csv("Weighted_Guilds.csv", header=T)
indicator<-read.csv("indicatorbydistance2.csv", header=T)
supp<-merge(weights, indicator, by='OTU')
write.csv(supp, "Supp3.csv")
Suppsimp<-supp[,c(1:4, 6:15)]
identity<-read.csv('UnitaandTBAS_indicator.csv', header=T)
#identity<-read.csv('Comparing_Unite_and_TBAS_Outputs.csv', header=T)
IDNodup <- identity[!duplicated(identity$OTU), ]
Finalsupp<-merge(Suppsimp, IDNodup, by='OTU', all=FALSE)
write.csv(Finalsupp, "SupplementaryFigINDICFINAL2.csv")
