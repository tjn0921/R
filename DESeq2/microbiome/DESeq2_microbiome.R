rm(list = ls())

mydata = read.table('data/otutab.txt', row.names = 1, header = T) %>%
  dplyr::select(seq(1:60))

coldata = data.frame(row.names = colnames(mydata),
                     group_list = rep(c('CK','Treatment'),each = 30))


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mydata, colData = coldata, design= ~group_list)
dds2 <- DESeq(dds) #标准化
res <- results(dds2, contrast=c("group_list",'CK','Treatment')) #差异分析结果

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)

# 火山图
resdata$g=ifelse(resdata$pvalue>0.5,'stable',
                 ifelse(resdata$log2FoldChange > 1,'up',
                        ifelse(resdata$log2FoldChange < -1,'down','stable')))
table(resdata$g)
resdata$v = -log10(resdata$pvalue)

library(ggplot2)
ggplot(data = resdata, 
       aes(x = log2FoldChange, 
           y = log2(baseMean), 
           color = g))+
  geom_point()+
  coord_flip()+
  theme_classic()+
  theme(panel.grid.major = element_line(),
        legend.title = element_blank())
