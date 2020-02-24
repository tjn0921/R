**朴素贝叶斯分类**是**贝叶斯分类**中最常见最简单的分类方法。相关的解释可以参照知乎的帖子：https://zhuanlan.zhihu.com/p/26262151

本教程重点是朴素贝叶斯分类在R中的实现。

----

# 数据导入

```R
# Libraries
library(naivebayes)
library(dplyr)
library(ggplot2)
library(psych)

# Data
data <- read.csv('binary.csv', header = T)
str(data)
xtabs(~admit+rank, data = data)
data$rank <- as.factor(data$rank)
data$admit <- as.factor(data$admit)
```

---

# 绘图展示几个变量之间的关系

```R
pairs.panels(data[-1])
```

![变量之间的关系](https://github.com/lixiang117423/R/raw/master/MachineLearning/NaiveBayesClassification/figures/1.png)

---

# 箱线图展示gpa和admitd的关系

```R
data %>%
         ggplot(aes(x=admit, y=gpa, fill = admit)) +
         geom_boxplot() +
         ggtitle("Box Plot")
```

![变量之间的关系](https://github.com/lixiang117423/R/raw/master/MachineLearning/NaiveBayesClassification/figures/2.png)

可以看出的是“1”的学生的GRE成绩的平均值高于“0”的学生。也可以绘制其他变量的箱线图。

---

# 展示密度分布

```R
data %>% ggplot(aes(x=gpa, fill = admit)) +
         geom_density(alpha=0.8, color= 'black') +
         ggtitle("Density Plot")
```

![密度分布](https://github.com/lixiang117423/R/raw/master/MachineLearning/NaiveBayesClassification/figures/3.png)

---

# 数据预处理

将数据划分成训练集和测试集：

```R
# Data Partition
set.seed(1234)
ind <- sample(2, nrow(data), replace = T, prob = c(0.8, 0.2))
train <- data[ind == 1,]
test <- data[ind == 2,]
```

---

# 构建模型

```R
# Naive Bayes Model
model <- naive_bayes(admit ~ ., data = train, usekernel = T)
model
```

---

# 模型可视化

```R
plot(model)
```

![模型可视化](https://github.com/lixiang117423/R/raw/master/MachineLearning/NaiveBayesClassification/figures/4.png)

---

# 用模型进行预测

## 预测训练集

````R
# Predict
p <- predict(model, train, type = 'prob')
head(cbind(p, train))

# Confusion Matrix - train data
p1 <- predict(model, train)
(tab1 <- table(p1, train$admit))
1 - sum(diag(tab1)) / sum(tab1)
````

```R
> head(cbind(p, train)) # 展示预测的结果的可能性
          0         1 admit gre  gpa rank
1 0.8528794 0.1471206     0 380 3.61    3
2 0.5621460 0.4378540     1 660 3.67    3
3 0.2233490 0.7766510     1 800 4.00    1
4 0.8643901 0.1356099     1 640 3.19    4
6 0.6263274 0.3736726     1 760 3.00    2
7 0.5933791 0.4066209     1 560 2.98    1

> (tab1 <- table(p1, train$admit))
   
p1    0   1
  0 203  69
  1  20  33
> 1 - sum(diag(tab1)) / sum(tab1)
[1] 0.2738462
```

预测训练集时，错误率在27.4%左右。

## 预测测试集

```R
# Confusion Matrix - test data
p2 <- predict(model, test)
(tab2 <- table(p2, test$admit))
1 - sum(diag(tab2)) / sum(tab2)
```

```R
> (tab2 <- table(p2, test$admit))
   
p2   0  1
  0 47 20
  1  3  5
> 1 - sum(diag(tab2)) / sum(tab2)
[1] 0.3066667
```

预测测试集时，错误率在31%左右。

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！