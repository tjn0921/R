rm(list = ls())

mydata = read.table('data/otutab.txt', row.names = 1, header = T) %>%
  dplyr::select(seq(1:20))

coldata = data.frame(row.names = colnames(mydata),
                     group_list = rep(c('CK','Treatment'),each = 10))


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mydata, colData = coldata, design= ~group_list)
dds2 <- DESeq(dds) #标准化
res <- results(dds2, contrast=c("group_list",'CK','Treatment')) #差异分析结果

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)

# 火山图



