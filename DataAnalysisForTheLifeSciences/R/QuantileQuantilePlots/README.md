# 为什么要进行正态性检验

进行*t*检验、方差分析、相关性分析等数据分析时，都要求数据服从正态分布或者近似正态分布。但是这个条件往往被忽略。为了保证数据满足上述方法的使用条件，对数据进行正态性检验是十分重要的。本文就简单介绍如何在R语言中对数据进行正态性检验。
本文的测试数据为R内嵌数据集**iris**。

# shapiro.test函数检验

R基础函数**shapiro.test()**可以对数据进行正态性检验。根据输出的*p-value*和*W*进行判断数据是否服从正态分布。*p-value*<0.05说明数据不符合正态分布；*W*越接近1说明数据越接近正态分布。

```{r}
suppressMessages(T)
shapiro.test(iris$Sepal.Length)
```

```R
> shapiro.test(iris$Sepal.Length)

	Shapiro-Wilk normality test

data:  iris$Sepal.Length
W = 0.97609, p-value = 0.01018
```



根据输出的*p-value*=0.01018<0.06,说明数据不符合正态分布。

# 直方图判断

直接根据直方图判断数据是否服从正态分布。

先看一个服从正态分布的例子

```{r}
set.seed(123)
hist(rnorm(100,mean = 5, sd = 1))
```

![服从正态分布](https://github.com/lixiang117423/PLANTOMIX/raw/master/20200309%E3%80%90R%E6%8A%80%E8%83%BD%E3%80%91R%E8%AF%AD%E8%A8%80%E4%B8%AD%E5%AF%B9%E6%95%B0%E6%8D%AE%E8%BF%9B%E8%A1%8C%E6%AD%A3%E6%80%81%E6%80%A7%E6%A3%80%E9%AA%8C/figures/1.png)

再看一个不服从正态分布的例子

```{r}
set.seed(123)
hist(runif(100, min = 2, max = 4))
```

![不服从正态分布](https://github.com/lixiang117423/PLANTOMIX/raw/master/20200309%E3%80%90R%E6%8A%80%E8%83%BD%E3%80%91R%E8%AF%AD%E8%A8%80%E4%B8%AD%E5%AF%B9%E6%95%B0%E6%8D%AE%E8%BF%9B%E8%A1%8C%E6%AD%A3%E6%80%81%E6%80%A7%E6%A3%80%E9%AA%8C/figures/2.png)

再看我们的测试数据

```{r}
hist(iris$Sepal.Length)
```

![测试数据](https://github.com/lixiang117423/PLANTOMIX/raw/master/20200309%E3%80%90R%E6%8A%80%E8%83%BD%E3%80%91R%E8%AF%AD%E8%A8%80%E4%B8%AD%E5%AF%B9%E6%95%B0%E6%8D%AE%E8%BF%9B%E8%A1%8C%E6%AD%A3%E6%80%81%E6%80%A7%E6%A3%80%E9%AA%8C/figures/3.png)

# Q-Q图判断

Q-Q图是根据数据的分位情况进行判断数据是否服从正态分布。使用R基础函数**qqnorm()**和**qqline()**即可绘制Q-Q图。x轴是理论分位数，Y轴是数据的分位数。如果数据符合正态分布，那数据点应该很好地拟合给出的直线。
```{r}
qqnorm(iris$Sepal.Length)
qqline(iris$Sepal.Length)
```

![Q-Q图](https://github.com/lixiang117423/PLANTOMIX/raw/master/20200309%E3%80%90R%E6%8A%80%E8%83%BD%E3%80%91R%E8%AF%AD%E8%A8%80%E4%B8%AD%E5%AF%B9%E6%95%B0%E6%8D%AE%E8%BF%9B%E8%A1%8C%E6%AD%A3%E6%80%81%E6%80%A7%E6%A3%80%E9%AA%8C/figures/4.png)

# 小节

数据正态性检验很重要。数据正态性检验的方法很多，掌握一两种即可。