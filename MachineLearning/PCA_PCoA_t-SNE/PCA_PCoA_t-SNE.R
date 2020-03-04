rm(list = ls())
mydata = iris


# PCA
pca = prcomp(mydata[,-5],
            center = T,
            scale. = T)
pca$center
pca$scale

library(ggbiplot)
pca.plot = ggbiplot(pca,
             obs.scale = 1,
             var.scale = 1,
             groups = mydata$Species,
             ellipse = T,
             circle = F,
             ellipse.prob = 0.68)
pca.plot = pca.plot + 
  scale_color_discrete(name = '') + 
  theme(legend.direction = 'horizontal',legend.position = 'top')+
  labs(title = 'PCA')

print(pca.plot)


# PCoA 
library(vegan)
pcoa = cmdscale(vegdist(mydata[,-5], method = 'bray'), k = nrow(mydata) - 1, eig = T)
pcoa_eig <- (pcoa$eig)[1:2]/sum(pcoa$eig)
library(ggplot2)

pcoa.plot = ggplot(data = as.data.frame(pcoa$points),mapping = aes(x = V1, y = V2, col = mydata$Species))+
  geom_point()+
  labs(x = paste('PCoA 1:', round(pcoa_eig[1],4)*100,'%'),
       y = paste('PCoA 2:', round(pcoa_eig[2],4)*100,'%'),
       title = 'PCoA')


# t-SNE
library(Rtsne)
mydata.tsne = unique(iris)
tsne = Rtsne(as.matrix(unique(mydata.tsne[,-5])), dims = 2,
             perplexity = floor((nrow(mydata)-1)/3.1),
             max_iter = 999)

tsne.plot = ggplot(as.data.frame(tsne[["Y"]]), mapping = aes(x = V1, y = V2, 
                                                 col = mydata.tsne$Species))+
  geom_point()+
  labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2', title = 't-SNE')

library(patchwork)

p = pca.plot | pcoa.plot | tsne.plot
p = p & theme(legend.position = 'none')
p


# PERMANOVA
# all
(perma = adonis2(mydata[,-5]~Species, data = mydata, method  = 'bray',
                 permutations = 999))
perma$sig = ifelse(perma[5]<0.01,'**', 
                   ifelse(perma[5]<0.05 & perma[5]>0.01, '*','N.S'))
# one by one
group_name <- as.character(unique(mydata$Species))
adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    data_ij <- subset(mydata, mydata$Species %in% c(group_name[i], group_name[j]))
    adonis_result_otu_ij <- adonis2(data_ij[,-5]~Species, data_ij, permutations = 999, method  = 'bray')
    adonis_result_otu_ij$cat = paste(group_name[i],'/',group_name[j])
    adonis_result_two <- rbind(adonis_result_two, adonis_result_otu_ij[1,])
    if (i == (length(group_name) - 1) & j == length(group_name)) {
      rownames(adonis_result_two) = adonis_result_two$cat
      adonis_result_two$sig = ifelse(adonis_result_two[5]<0.01,'**', 
                                     ifelse(adonis_result_two[5]<0.05 & adonis_result_two[5]>0.01, '*','N.S'))
      adonis_result_two = adonis_result_two[,-(ncol(adonis_result_two)-1)]
      adonis_result_two <- as.data.frame(adonis_result_two)
    }
  }
}
adonis_result_all = rbind(perma[1,], adonis_result_two)
