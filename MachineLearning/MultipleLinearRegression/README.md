**多元线性回归**（Multiple Linear Regression）是**简单线性回归** 的变形。下面这张图可以看出两者的一些区别：

![一元线性回归和多元线性回归 区别](https://github.com/lixiang117423/R/raw/master/MachineLearning/MultipleLinearRegression/figures/1.png)

---

本教程数据下载地址：https://github.com/lixiang117423/R/blob/master/MachineLearning/MultipleLinearRegression/vehicle.csv

---

# 导入数据建立模型

```R
rm(list = ls())

data = read.csv('vehicle.csv')

# multiple linear regression
res = lm(lc ~ Mileage + lh, data = data) # 选择其中两个变量进行建模
summary(res)
```

```R
> summary(res)

Call:
lm(formula = lc ~ Mileage + lh, data = data)

Residuals:
    Min      1Q  Median      3Q     Max 
-672.79  -14.73   -0.62   12.89  741.05 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.375e+00  2.218e+00    0.62    0.535    
Mileage     -8.475e-05  6.622e-05   -1.28    0.201    
lh           7.355e+01  4.155e-01  177.01   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 48.62 on 1621 degrees of freedom
Multiple R-squared:  0.951,	Adjusted R-squared:  0.951 
F-statistic: 1.574e+04 on 2 and 1621 DF,  p-value: < 2.2e-16
```

根据这个模型，我们可以得到的回归方程是：

$$lc = 1.375 - 0.00008475×Mileage + 0.7335×lh$$

而且，可以看到，变量**Mileage**在模型中的“作用”并不显著。那是不是可以省去这个变量呢？

---

# 一元和多元的比较

我们可以把变量**Mileage**不加入到模型中，那开始的模型就变成一元线性回归了。

```R
reduced = lm(lc~lh, data = data)
summary(reduced)
```

```R
> summary(reduced)

Call:
lm(formula = lc ~ lh, data = data)

Residuals:
    Min      1Q  Median      3Q     Max 
-670.09  -14.72   -0.32   12.96  742.70 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  -0.2359     1.8262  -0.129    0.897    
lh           73.5088     0.4144 177.387   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 48.63 on 1622 degrees of freedom
Multiple R-squared:  0.951,	Adjusted R-squared:  0.9509 
F-statistic: 3.147e+04 on 1 and 1622 DF,  p-value: < 2.2e-16
```

此时得到的就是一元线性回归。

**比较二元线性回归和一元线性回归是否有差异**：

```R
anova(res, reduced)
> anova(res, reduced)
Analysis of Variance Table

Model 1: lc ~ Mileage + lh
Model 2: lc ~ lh
  Res.Df     RSS Df Sum of Sq      F Pr(>F)
1   1621 3831889                           
2   1622 3835760 -1     -3871 1.6376 0.2008
```

可以看到，*p=0.2008 > 0.05*。说明一元回归和二元回归是没有差异的。

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！