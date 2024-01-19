library(tidyverse)
library(vegan)
library(plyr)
AllGoLIFEFULL<-read.csv("ALL-GOLIFE-TBAS.csv", header = TRUE)
ALLTBASOTU<-AllGoLIFEFULL[,c("otu", "Query.sequence")]
OTUtableAll<-ALLTBASOTU%>%
  count(otu, Query.sequence)%>%
  pivot_wider(names_from=otu, values_from=n, values_fill=0)

###Filter out the singletons
sumtable<-function(x){
  if(is.numeric(x)){
    sum(x) > 1
  } else {
    TRUE
  }
}

AllOTUfilter<-OTUtableAll[, sapply(OTUtableAll,  sumtable)]
All_Golife_OTU_filtered <- AllOTUfilter %>%
  filter(rowSums(select(., starts_with("OTU")), na.rm = TRUE) != 0)

write.csv(All_Golife_OTU_filtered, "AllGoLifeOTUNov_24.csv")

MetaAll<-read.csv("GoLifeAllMetaCorrected.csv", header=TRUE)#%>%
#rename(Query.sequence = Code)
MetaAllFilter<-merge(MetaAll, All_Golife_OTU_filtered, by.x="Query.sequence", all=FALSE)


META_no_duplicates <- MetaAllFilter[!duplicated(MetaAllFilter$Query.sequence), ]
write.csv(META_no_duplicates, "AllGoLifeOTUandMeta_Nov24.csv")

JustMetaAllFilter<-META_no_duplicates[,1:12]


###First running an NMDS to see how they cluster by region
AllOTUTableNMDS<-column_to_rownames(All_Golife_OTU_filtered, var="Query.sequence")
# AbundanceGoLife<-decostand(AllOTUTableNMDS, method="total")
# write.csv(AllOTUTableNMDS, "allcheckingtable.csv")


dm.method <- 'bray'
dmall <- vegdist(AllOTUTableNMDS, method=dm.method)
Allotu.nmds <- metaMDS(dmall,
                       k = 2,
                       maxit = 999,
                       trymax = 500,
                       wascores = TRUE)
stressplot(Allotu.nmds)


###Couldnt analyze sample to sample trying to analyze by site 
#MetaPlusOTU<-merge(MetaAll, OTUtableAll, by.x="Query.sequence")
#MetaPlusOTU<-MetaPlusOTU[!duplicated(MetaPlusOTU$Query.sequence), ]
#write.csv(MetaPlusOTU, "Combiningsitesall.csv")
###excel simplify site names with transect into
# MetaPlusOTU2<-read.csv("Combiningsitesall.csv", header=T)
# OTU_Nove<-MetaPlusOTU2[,c(2,15:178)]
OTUbySite<-META_no_duplicates[,c(2,15:178)] %>%
  group_by(Site)%>%
  summarize(across(starts_with("OTU"), sum, na.rm = TRUE))
library(dplyr)

##Get Rid of the singletons
sumtable<-function(x){
  if(is.numeric(x)){
    sum(x) > 1
  } else {
    TRUE
  }
}


AllOTUfilter<-OTUbySite[, sapply(OTUbySite,  sumtable)]


###get rid of empty rows
filtered_OTU <- AllOTUfilter %>%
  filter(rowSums(select(., starts_with("OTU")), na.rm = TRUE) != 0)
AllOTUNMDS<-column_to_rownames(filtered_OTU, var="Site")
ALLOTUAbund<-decostand(AllOTUNMDS, method='total')
##Richness
##dont load plyr until now or summarize across function will not work 
library(plyr)
Richness_site<-apply(filtered_OTU[,-1]>0,1,sum)
Richness<-ddply(filtered_OTU,~Site,function(x) {
  data.frame(RICHNESS=sum(x[-1]>0))
})
###Now ln transform the values 
Richness$ln_richness <- log(Richness$RICHNESS)
###Getting rid of sites where ln Richness is 0
Richness_no_zero <- Richness[!Richness$ln_richness == "0", ]
##Merge with Metadata
Site_info<-read.csv("GoLife_Meta_All_nov.csv", header=T)
Site_info_Rich<-merge(Site_info, Richness_no_zero, by="Site", all=FALSE)
Go_climate<-read.csv("GoLifeAllMetaCorrected_nov.csv", header=TRUE)
All_Meta_Rich<-merge(Site_info_Rich, Go_climate, by="Strain", all=FALSE)
All_Meta_Rich <- All_Meta_Rich[!duplicated(All_Meta_Rich$Strain), ]
###getting rid of the two highest richness sites because they are skewing the data left
All_Meta_Rich <- All_Meta_Rich[!All_Meta_Rich$Site == "Tucson Mountains", ]
All_Meta_Rich <- All_Meta_Rich[!All_Meta_Rich$Site == "Tucson Mountains (Saguaro)", ]
###To Add to your appendix
write.csv(All_Meta_Rich_NoOutlier, "Final_Meta_For_appendix2.csv")
anova_Richness <- aov(All_Meta_Rich$ln_richness ~ Location, data=All_Meta_Rich)
summary(anova_Richness)
##Location Random Factor
library(nlme)
anova_temp <- aov(All_Meta_Rich$ln_richness ~ Mean_annual_temperature + Error(Location), data=All_Meta_Rich)
summary(anova_temp)
### Mean_precipitation_driest_quarter
anova_Mean_precipitation_driest_quarter <- aov(All_Meta_Rich$ln_richness ~ Mean_precipitation_driest_quarter + Error(Location), data=All_Meta_Rich)
summary(anova_Mean_precipitation_driest_quarter)
####Mean_precipitation_wettest_quarter
anova_Mean_precipitation_wettest_quarter <- aov(All_Meta_Rich$ln_richness ~ Mean_precipitation_wettest_quarter + Error(Location), data=All_Meta_Rich)
summary(anova_Mean_precipitation_wettest_quarter)
###Mean_annual_precipitation
anova_Mean_annual_precipitation <- aov(All_Meta_Rich$ln_richness ~ Mean_annual_precipitation + Error(Location), data=All_Meta_Rich)
summary(anova_Mean_annual_precipitation)

### Look for pairwise correlations
###Making the cor_mat
##simplify to numerical variables 
Cor_Max_df<- All_Meta_Rich[,-c(1:8, 11)]
Cor_Max_df$Location[Cor_Max_df$Location == "Southern_Africa"] <-1
Cor_Max_df$Location[Cor_Max_df$Location == "Arizona"] <-2
Cor_Max_df$Location[Cor_Max_df$Location == "Chile"] <-3
Cor_Max_df$Location<- as.numeric(Cor_Max_df$Location)
###Check variable types
cor_matrix <- cor(Cor_Max_df, method = "pearson")  # You can use "spearman" or "kendall" for other correlation methods
##Confidence Intervals
print(cor_matrix)
###
install.packages("psych")
library(psych)
correlation_mat<-corr.test(Cor_Max_df, method = "pearson")
print(correlation_mat, short=FALSE)

write.csv(correlation_mat$r, "Correlation_matrix-Golife.csv")
write.csv(correlation_mat$p, "Pvals_Correlation_matrix-Golife.csv")
write.csv(correlation_mat$ci, "CI_Correlation_matrix-Golife.csv")
##Make it pretty 
lowerCor(Cor_Max_df, method = "pearson")
install.packages("apaTables")
library(apaTables)
apa.cor.table(Cor_Max_df, "Correlation_Table_Golife.doc")
install.packages("metan")
library(metan)
All<-corr_coef(Cor_Max_df)
All
plot(All)
##Structural equation model 
install.packages("lavaan", dependencies=TRUE)
library(lavaan)
cov(Cor_Max_df)
mvnout <- mardia(Cor_Max_df)
mvnout$ln_richness
print(mvnout)
## Saw that ln Richness skews to the left so I got ris of outliers on the high end. 
###PCA Time since our variables are correlated to cluster them into variables 
Climate_only<- select(Cor_Max_df, c(5:8))
###Need to filter Datasets
OTU_Site_NoOu<- AllOTUfilter[!AllOTUfilter$Site == "Tucson Mountains", ]
OTU_Site_NoOu<- AllOTUfilter[!AllOTUfilter$Site == "Tucson Mountains (Saguaro)", ]
OTU_Site_NoOu<-OTU_Site_NoOu %>%
  column_to_rownames(var = "Site")
zero_sum_columns <- colSums(OTU_Site_NoOu) == 0
OTU_Site_NoOu <- OTU_Site_NoOu[, !zero_sum_columns]
#PCA<-rda(decostand(OTU_Site_NoOu, method="hellinger", scale=T, Center=T))
#PCA
#pcaData <- as.data.frame(PCA$x[, 1:2])
ClimdatPCA<-All_Meta_Rich[,c(2,10, 14:17)]
pca <- prcomp(All_Meta_Rich[,14:17], scale = TRUE)
pca2<-prcomp(ClimdatPCA[,2:5], scale = TRUE)
autoplot(pca)
autoplot(pca)
str(pca)
pcaData <- as.data.frame(pca$x[, 1:2])
pcaData <- cbind(pcaData, All_Meta_Rich$Site, All_Meta_Rich$Location, All_Meta_Rich$ln_richness)
PCA_clim_site1<-merge(PCA_clim_site, Richness, by='Site', all=FALSE)

colnames(pcaData) <- c("PC1", "PC2", "Site_Name", "Location", "ln_Richness")
##Statistics on our PCA data
anova_PCA1 <- aov(pcaData$ln_Richness ~ PC1, data=pcaData)
summary(anova_PCA1)
anova_PCA2 <- aov(pcaData$ln_Richness ~ PC2, data=pcaData)
summary(anova_PCA2)
###Generalized linear model
Temp<-lm(ln_richness ~ Mean_annual_temperature, data=PCA_clim_site1)
water<-lm(ln_richness ~ Mean_annual_precipitation, data=PCA_clim_site1)
summary(Temp)
summary(water)

colnames(pcaData) <- c("PC1", "PC2", "Site_Name", "Location")
ggplot(pcaData) +
  aes(PC1, PC2, color = Location, shape = Location) + 
  geom_point(size = 3) 
coord_fixed() # fixing coordinates
library(ggfortify)
library(ggrepel)
result <- ClimdatPCA %>%
  group_by(Site, Location) 
summarize(
  annual_temperature = mean(Mean_annual_temperature),
  precipitation_driest_quarter = mean(Mean_precipitation_driest_quarter), 
  precipitation_wettest_quarter = mean(Mean_precipitation_wettest_quarter), 
  precipitation_wettest_quarter = mean(Mean_annual_precipitation))

PCA_clim_site<-result[!duplicated(result$Site), ]
###FOR PCA 
pca_sites <- prcomp(PCA_clim_site[,3:6], scale = TRUE)

autoplot(pca, data = PCA_clim_site, colour = 'Location', size=3,
         loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, loadings.label.size = 3, loadings.label.color='black')+
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2)

pca2<-autoplot(pca, data = All_Meta_Rich, 
               loadings = TRUE, loadings.colour = 'blue', 
               loadings.label = TRUE, loadings.label.size = 3, loadings.label.color='black')+
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2)

pca2

unique_variables <- unique(All_Meta_Rich$Plant.or.Lichen.host)
print(unique_variables)

All_Meta_Rich$Plant.or.Lichen.host<-  trimws(All_Meta_Rich$Plant.or.Lichen.host)
All_Meta_Rich$Plant.or.Lichen.host <- gsub("Lichen-saxicolous", "Lichen", All_Meta_Rich$Plant.or.Lichen.host, ignore.case = TRUE)



#+ theme(panel.grid.major = element_line(linetype = "blank"))




library(ggplot2)
library(plotly)
library(ggfortify)
install.packages('ggfortify')

p <- autoplot(PCA, data = All_Meta_Rich, colour = 'Location')

##NMDS
dm.method <- 'bray'
dmall <- vegdist(AllOTUNMDS, method=dm.method)
Allotu.nmds <- metaMDS(dmall,
                       k = 2,
                       maxit = 999,
                       trymax = 500,
                       wascores = TRUE)
stressplot(Allotu.nmds)
##Meta simplifies
Meta.NMDS<-merge(MetaAll, filtered_OTU, by.x="Site")
Meta.NMDS<-Meta.NMDS[!duplicated(Meta.NMDS$Site), ]
plot(Allotu.nmds)

nmds.scores <- as.data.frame(scores(Allotu.nmds, display = 'sites'))
nmds.scores$site <- rownames(nmds.scores) 
nmds.scores$Site <- nmds.scores$site
nmds.scores$site <- NULL
nmds.scores <- left_join(nmds.scores, Meta.NMDS, by = 'Site')
nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(size=3, aes(color=Mean_annual_precipitation, shape=Location)) +
  labs(x = "NMDS1", colour = "Mean Annual Percipitation", shape="Location") 
nmds_plot
nmds_plot2 <- ggplot(nmds.scores,
                     aes(x = NMDS1,
                         y = NMDS2)) +
  geom_point(size=3, aes(color=Host.clade, shape=Location)) +
  labs(x = "NMDS1", colour = "Host Type", shape="Location") 
nmds_plot2
nmds_plot3 <- ggplot(nmds.scores,
                     aes(x = NMDS1,
                         y = NMDS2)) +
  geom_point(aes(color=Mean_annual_precipitation, shape=Location, size=Mean_annual_temperature)) +
  labs(x = "NMDS1", colour = "Mean_annual_precipitation", shape="Location", size= "Mean annual temperature") +scale_colour_gradientn(colours=rainbow(4))
nmds_plot3
nmds_plot4 <- ggplot(nmds.scores,
                     aes(x = NMDS1,
                         y = NMDS2)) +
  geom_point(aes(color=Mean_precipitation_wettest_quarter, shape=Location, size=Mean_annual_precipitation)) +
  labs(x = "NMDS1", colour = "Precipiation wettest quarter", size= "Mean annual precipiation", shape="Location") + ggtitle("All GoLife sites Ordination", )+ theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_gradientn(colours=rainbow(4))
nmds_plot4 <- ggplot(nmds.scores,
                     aes(x = NMDS1,
                         y = NMDS2)) +
  geom_point(aes(color = Mean_precipitation_wettest_quarter, shape = Location, size = Mean_annual_precipitation)) +
  labs(x = "NMDS1", color = "Precipitation wettest quarter", size = "Mean annual precipitation", shape = "Location") +
  ggtitle("All GoLife sites Ordination") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_gradientn(colours = rainbow(4))

print(nmds_plot4)

nmds_plot4
adonis2(vegdist(AllOTUNMDS, method="bray")~ Mean_precipitation_wettest_quarter + Location, data=Meta.NMDS, permutations=9999)
adonis2(vegdist(AllOTUNMDS, method="bray") ~ Mean_annual_precipitation + region, data=Meta.NMDS, permutations=9999)
adonis2(vegdist(AllOTUNMDS, method="bray") ~ Host.clade + Location, data=Meta.NMDS, permutations=9999)
adonis2(vegdist(AllOTUNMDS, method="bray") ~ Mean_annual_temperature + Location, data=Meta.NMDS, permutations=9999)
write.csv(AllOTUTableNMDS, "GoLifeOTUsbysite.csv")

###Combining the Chile sites
MetaPlusOTU3<-read.csv("Combiningsitesallupdated8_27_23.csv", header=T)
OTUbySite2<-MetaPlusOTU3[,c(2,13:176)] %>%
  group_by(Site) %>%
  summarize(across(starts_with("OTU"), sum, na.rm = TRUE))
##Get Rid of the singletons
sumtable<-function(x){
  if(is.numeric(x)){
    sum(x) > 1
  } else {
    TRUE
  }
}


AllOTUfilter2<-OTUbySite2[, sapply(OTUbySite2,  sumtable)]


###get rid of empty rows
filtered_OTU2 <- AllOTUfilter2 %>%
  filter(rowSums(select(., starts_with("OTU")), na.rm = TRUE) != 0)
AllOTUNMDS2<-column_to_rownames(filtered_OTU2, var="Site")
dm.method <- 'bray'
dmall2 <- vegdist(AllOTUNMDS2, method=dm.method)
Allotu.nmds2 <- metaMDS(dmall2,
                        k = 2,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)
stressplot(Allotu.nmds2)
meta2<-MetaPlusOTU3[,1:13]
metabysite<-meta2[,c(2,3,9,10,11,12)]%>%
  group_by(Site)
metabysite_nodup<-metabysite[!duplicated(metabysite$Site), ]
nmds.scores2 <- as.data.frame(scores(Allotu.nmds2, display = 'sites'))
nmds.scores2$site <- rownames(nmds.scores2) 
nmds.scores2$Site <- nmds.scores2$site
nmds.scores2$site <- NULL
nmds.scores2 <- left_join(nmds.scores2, meta2, by = 'Site')
nmds_plot_d2 <- ggplot(nmds.scores2,
                       aes(x = NMDS1,
                           y = NMDS2)) +
  geom_point(size=3, aes(color=Mean_annual_precipitation, shape=Location)) +
  labs(x = "NMDS1", colour = "Mean Annual Percipitation", shape="Location") 
nmds_plot_d2
nmds_plot_d3 <- ggplot(nmds.scores2,
                       aes(x = NMDS1,
                           y = NMDS2)) +
  geom_point(size=3, aes(color=Host.clade, shape=Location)) +
  labs(x = "NMDS1", colour = "Host Type", shape="Location") 
nmds_plot_d3
nmds_plot_d4 <- ggplot(nmds.scores2,
                       aes(x = NMDS1,
                           y = NMDS2)) +
  geom_point(aes(color=Mean_annual_precipitation, shape=Location, size=Mean_annual_temperature)) +
  labs(x = "NMDS1", colour = "Mean_annual_precipitation", shape="Location", size= "Mean annual temperature") +scale_colour_gradientn(colours=rainbow(4))
nmds_plot_d4
nmds_plot_d5 <- ggplot(nmds.scores2,
                       aes(x = NMDS1,
                           y = NMDS2)) +
  geom_point(aes(color=Mean_precipitation_wettest_quarter, shape=Location, size=Mean_annual_precipitation)) +
  labs(x = "NMDS1", colour = "Precipiation wettest quarter", size= "Mean annual precipiation", shape="Location") + ggtitle("All GoLife sites Ordination", )+ theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_gradientn(colours=rainbow(4))
nmds_plot_d5 
nmds_plot_d6<- ggplot(nmds.scores2,
                      aes(x = NMDS1,
                          y = NMDS2)) +
  geom_point(aes(color = Mean_precipitation_wettest_quarter, shape = Location, size = Mean_annual_precipitation)) +
  labs(x = "NMDS1", color = "Precipitation wettest quarter", size = "Mean annual precipitation", shape = "Location") +
  ggtitle("All GoLife sites Ordination") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_gradientn(colours = rainbow(4))
nmds_plot_d6
###NMDS AGAIN FOR SITES
dm.method <- 'bray'
dmall <- vegdist(AllOTUNMDS, method=dm.method)
Allotu.nmds <- metaMDS(dmall,
                       k = 2,
                       maxit = 999,
                       trymax = 500,
                       wascores = TRUE)
stressplot(Allotu.nmds)
##Meta simplifies
Meta.NMDS<-merge(MetaAll, filtered_OTU, by.x="Site")
Meta.NMDS<-Meta.NMDS[!duplicated(Meta.NMDS$Site), ]
plot(Allotu.nmds)

nmds.scores <- as.data.frame(scores(Allotu.nmds, display = 'sites'))
nmds.scores$site <- rownames(nmds.scores) 
nmds.scores$Site <- nmds.scores$site
nmds.scores$site <- NULL
nmds.scores <- left_join(nmds.scores, Meta.NMDS, by = 'Site')
nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(size=3, aes(color=Mean_annual_precipitation, shape=Location)) +
  labs(x = "NMDS1", colour = "Mean Annual Percipitation", shape="Location") 
nmds_plot

library(plyr)
Richness_site<-apply(filtered_OTU2[,-1]>0,1,sum)
Richness<-ddply(filtered_OTU2,~Site,function(x) {
  data.frame(RICHNESS=sum(x[-1]>0))
})
Richnessbysite<-merge(Richness, metabysite_nodup, by.x="Site")
write.csv(Richnessbysite, "richness_by_site30_8_23.csv")
###Rename col in OTU by site
OTU_Nove<-OTUtableAll%>%
  rename(Strain=Query.sequence)
OTUplusMeta2<-merge( All_Meta_Rich, OTUtableAll, by.x='Strain')###DATASET GOT SMALLER CUZ WE TOOK OUT 2 OUTLIERS IN RICHNESS ANALYSIS
###LETS MAKE AN NMDS THAT SHOWS PLANT INFO
MetabyHost.Type<-OTUplusMeta2[,c(2:3, 6:10, 14:17)] %>%
  group_by(Site.Host)
MetabyHost.Type<-MetabyHost.Type[!duplicated(MetabyHost.Type$Site.Host), ]

#%>%
OTUbyHost.Type<-OTUplusMeta2[,c(7, 17:181)] %>%
  group_by(Site.Host)%>%
  summarize(across(starts_with("OTU"), sum, na.rm = TRUE))%>%
  column_to_rownames(var="Site.Host")

zero_sum_columns2 <- colSums(OTUbyHost.Type) == 0
OTUbyHost.Type <- OTUbyHost.Type[, !zero_sum_columns2]

detach("package:plyr", unload = TRUE)

library(dplyr)

Meta.NMDS<-merge(MetaAll, filtered_OTU, by.x="Site")
Meta.NMDS<-Meta.NMDS[!duplicated(Meta.NMDS$Site), ]
plot(Allotu.nmds)

dm.method <- 'jaccard'
dmall5 <- vegdist(OTUbyHost.Type, method=dm.method)
Allotu.nmds5 <- metaMDS(dmall5,
                        k = 2,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)
stressplot(Allotu.nmds5)
###Trying PCOA
install.packages("ecodist")
bray_curtis_dist <- vegdist(OTUbyHost.Type, method = "bray")
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
bray_curtis_plot <- ggplot(data = Bray_Meta_pcoa_df, aes(x=pcoa1, y=pcoa2, color= Host, shape=Location)) +
  geom_point(size=3) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) 
bray_curtis_plot

MetabyHost.Type<-MetabyHost.Type %>%
  column_to_rownames(var="Site.Host")
###Lets see if statistically significant
Permanova_HostType<- adonis2(vegdist(OTUbyHost.Type, method = "bray")~ Plant.or.Lichen.host, data=MetabyHost.Type, strata=MetabyHost.Type$Location, permutations = 9999)

summary(Permanova_HostType)
Permanova_HostType

install.packages('DESeq2')
library(DESeq2)
Bray_Meta_pcoa_df <- cbind(bray_curtis_pcoa_df,
                           Host = MetabyHost.Type$Plant.or.Lichen.host, 
                           Location= MetabyHost.Type$Location)

Permanova_HostType2<- adonis2(vegdist(OTUbyHost.Type, method = "bray")~ Plant.or.Lichen.host*Location, data=MetabyHost.Type, permutations = 9999)

Permanova_HostType2



nmds.scores <- as.data.frame(scores(OTUbyHost.Type, display = 'sites'))
nmds.scores$site <- rownames(nmds.scores) 
nmds.scores$Site <- nmds.scores$site
nmds.scores$site <- NULL
nmds.scores <- left_join(nmds.scores, MetabyHost.Type, by = 'Site')
nmds_plot <- ggplot(nmds.scores,
                    aes(x = NMDS1,
                        y = NMDS2)) +
  geom_point(size=3, aes(color=Mean_annual_precipitation, shape=Location)) +
  labs(x = "NMDS1", colour = "Mean Annual Percipitation", shape="Location") 
nmds_plot
