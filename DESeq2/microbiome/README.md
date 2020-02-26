**DESeq2**常常用于**RNA-Seq**数据分析,也有用于微生物组数据分析的$^{[1]}$。

![文献图](https://github.com/lixiang117423/R/raw/master/DESeq2/microbiome/2.png)

**DESeq**可以直接处理**FPKM**、**FTPM**、**TPM**及**CPM**数据，在处理微生物数据时用**OTU**数据。

下面我们就用一篇文献的数据进行简单的测试。（测试文件下载地址：https://github.com/lixiang117423/R/blob/master/DESeq2/microbiome/data/otutab.txt）

```R
rm(list = ls())

# 导入数据并选择数据
mydata = read.table('data/otutab.txt', row.names = 1, header = T) %>%
  dplyr::select(seq(1:60))

# 构建分组矩阵
coldata = data.frame(row.names = colnames(mydata),
                     group_list = rep(c('CK','Treatment'),each = 30))


library(DESeq2)

# 构建DESeq2用的举证
dds <- DESeqDataSetFromMatrix(countData = mydata, colData = coldata, design= ~group_list)

# 标准化数据
dds2 <- DESeq(dds)

# 提取差异分析结果
res <- results(dds2, contrast=c("group_list",'CK','Treatment'))

# 合并差异分析结果和表达矩阵
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)

# 火山图
# pvalue选择0.5是为了有更多的点作图便于展示，实际研究按需选择
resdata$g=ifelse(resdata$pvalue>0.5,'stable',
                 ifelse(resdata$log2FoldChange > 1,'up',
                        ifelse(resdata$log2FoldChange < -1,'down','stable')))
table(resdata$g)

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

```

最终得到如下的图：

![差异微生物示意图](https://github.com/lixiang117423/R/raw/master/DESeq2/microbiome/1.png)



---

**参考文献：**

[1]	EDWARDS J, JOHNSON C, SANTOS-MEDELLíN C et al. Structure, variation, and assembly of the root-associated microbiomes of rice[J]. ***Proceedings of the National Academy of Sciences***, 2015: 112: E911-E920.

---
**致谢：**感谢**StatQuest**视频栏目。

