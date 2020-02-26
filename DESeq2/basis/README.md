**DESeq2**常常用于**RNA-Seq**数据分析,也有用于微生物组数据分析的$^{[1]}$。**DESeq**可以直接处理**FPKM**、**FTPM**、**TPM**及**CPM**数据，在处理微生物数据时用**OTU**数据。
本贴就简单说说**DESeq2**是如何校正数据的。

---
假设有2个样本，6个基因：

![样本信息](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/1.png)

可以发现每个基因在样本2中的表达量都是样本1中的2倍。这个很可能不是由于样本本身的生物学因素引起的，而是由于测序的深度等影响的。需要对样本基因的表达量进行校正。

---
# Step1：取对数值
假设有如下的基因表达情况：

![样本信息](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/2.png)

先对每个表达量取以*e*为底数的对数：

![取对数](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/3.png)

# Step2：求每个基因表达量的均值
对每个基因（没行）求均值：

![每个基因表达量求均值](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/4.png)

# Step3：筛选表达量非零的基因
基因A在取对数值后出现**负无穷大**，所以将这个基因暂时剔除。

![筛选基因](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/5.png)

# Step4：求表达量的中位数
求出每个样本基因表达的中位数：

![求中位数](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/6.png)

# Step5：求出scaleing factor
根据公式$scaleing factor = e^{median}$求出每个样本基因表达量的scaleing factor。

![求scaleing factor](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/7.png)

# Step6：每个表达量的值除以scaleing factor
将每个样品中每个基因的表达量除以样品对应的scaleing factor即可。

![校正后的表达矩阵](https://github.com/lixiang117423/R/raw/master/DESeq2/basis/figures/8.png)

----
# 参考文献
[1]	EDWARDS J, JOHNSON C, SANTOS-MEDELLíN C et al. Structure, variation, and assembly of the root-associated microbiomes of rice[J]. ***Proceedings of the National Academy of Sciences***, 2015: 112: E911-E920.

---
**致谢：**感谢**StatQuest**视频栏目。