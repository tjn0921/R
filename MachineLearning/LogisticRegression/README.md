**逻辑回归**和**线性回归**类似。线性回归是使用连续变量对另外的变量进行预测。逻辑回归可以使用连续变量，也可以使用分类变量对别的变量进行预测。预测的结果通常是**YES**或**NO**。

![线性回归与逻辑回归](https://github.com/lixiang117423/R/raw/master/MachineLearning/LogisticRegression/figure/1.jpg)

---

**本教程使用GRE、GPA和RANK三个变量使用逻辑回归预测学生的申请是否被批准。**数据下载地址：https://github.com/lixiang117423/R/blob/master/MachineLearning/LogisticRegression/binary.csv。

----

# 数据导入处理

```R
# Read data file
mydata <- read.csv('binary.csv', header = T)
str(mydata)
mydata$admit <- as.factor(mydata$admit) # 设置成因子
mydata$rank <- as.factor(mydata$rank) # 设置成因子

# Two-way table of factor variables
xtabs(~admit + rank, data = mydata) # 查看数据

# 将数据划分成训练集和测试集
# Partition data - train (80%) & test (20%)
set.seed(1234)
ind <- sample(2, nrow(mydata), replace = T, prob = c(0.8, 0.2))
train <- mydata[ind==1,]
test <- mydata[ind==2,]
```

---

# 构建模型

```R
mymodel <- glm(admit ~ gre + gpa + rank, data = train, family = 'binomial')
summary(mymodel)
```

```R
> summary(mymodel)

Call:
glm(formula = admit ~ gre + gpa + rank, family = "binomial", 
    data = train)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.5873  -0.8679  -0.6181   1.1301   2.1178  

Coefficients:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept) -5.009514   1.316514  -3.805 0.000142 ***
gre          0.001631   0.001217   1.340 0.180180    
gpa          1.166408   0.388899   2.999 0.002706 ** 
rank2       -0.570976   0.358273  -1.594 0.111005    
rank3       -1.125341   0.383372  -2.935 0.003331 ** 
rank4       -1.532942   0.477377  -3.211 0.001322 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 404.39  on 324  degrees of freedom
Residual deviance: 369.99  on 319  degrees of freedom
AIC: 381.99

Number of Fisher Scoring iterations: 4
```

从运行的结果中可以看出，**gre**在模型中不显著（P=0.180180）。所以将**gre**剔除，重新构建模型。

---

# 构建新模型

```R
mymodel <- glm(admit ~ gpa + rank, data = train, family = 'binomial')
summary(mymodel)
```

```R
> summary(mymodel)

Call:
glm(formula = admit ~ gpa + rank, family = "binomial", data = train)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.5156  -0.8880  -0.6318   1.1091   2.1688  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -4.7270     1.2918  -3.659 0.000253 ***
gpa           1.3735     0.3590   3.826 0.000130 ***
rank2        -0.5712     0.3564  -1.603 0.108976    
rank3        -1.1645     0.3804  -3.061 0.002203 ** 
rank4        -1.5642     0.4756  -3.289 0.001005 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 404.39  on 324  degrees of freedom
Residual deviance: 371.81  on 320  degrees of freedom
AIC: 381.81

Number of Fisher Scoring iterations: 4
```

可见剔除**gre**变量后，整个模型更好了。

---

# 使用模型进行预测

## 使用模型预测训练集

```R
# Prediction
p1 <- predict(mymodel, train, type = 'response')
head(p1)
head(train)

# Misclassification error - train data
pred1 <- ifelse(p1>0.5, 1, 0)
tab1 <- table(Predicted = pred1, Actual = train$admit)
tab1
1 - sum(diag(tab1))/sum(tab1)
```

```R
> tab1
         Actual
Predicted   0   1
        0 208  73
        1  15  29
> 1 - sum(diag(tab1))/sum(tab1)
[1] 0.2707692
```

错误率是**27.1%**。

## 使用模型预测测试集

```R
# Misclassification error - test data
p2 <- predict(mymodel, test, type = 'response')
pred2 <- ifelse(p2>0.5, 1, 0)
tab2 <- table(Predicted = pred2, Actual = test$admit)
tab2
1 - sum(diag(tab2))/sum(tab2)
```

```R
> tab2
         Actual
Predicted  0  1
        0 48 20
        1  2  5
> 1 - sum(diag(tab2))/sum(tab2)
[1] 0.2933333
```

错误率为**29.3&**。

---

# 评估模型性能

这一步大概是评估模型的性能，我也没太看懂。

```R
# Goodness-of-fit test
with(mymodel, pchisq(null.deviance - deviance, df.null-df.residual, lower.tail = F))
[1] 1.450537e-06
```

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！