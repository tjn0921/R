**支持向量机**（Support Vector Machine），简称**SVM**，是一种二分类模型。先不探究算法的数学推理，直接看如何在R语言中实现。

----

# 数据探索

使用的数据集是R内嵌的数据集**iris**。先以**Petal.Width**为横坐标，**Petal.Length**为纵坐标，查看数据：

```R
# Data
data("iris")

library(ggplot2)
qplot(Petal.Length,Petal.Width,data = iris,color = Species)
```

![不同种之间分离状况](https://github.com/lixiang117423/R/raw/master/MachineLearning/SupportVectorMachine/figures/1.png)

可以看出，**setosa**和其他两个完全分开，**versicolor**和**virginica**有点交叉。

---

# 构建SVM模型

```R
# Support Vector Machine
library(e1071)
mymodel = svm(Species~., data = iris, kernel = 'radial') # kernal有几个选项，不同的选项准确度不一样
summary(mymodel)
plot(mymodel, data = iris,
     Petal.Width~Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
```

![SVM模型](https://github.com/lixiang117423/R/raw/master/MachineLearning/SupportVectorMachine/figures/2.png)

# 查看模型准确度

```R
# Confusion Matrix and Misclassification Error
pred = predict(mymodel, iris)
tab = table(Predicted = pred, Actual = iris$Species)
1-sum(diag(tab))/sum(tab)
```

可以看到：

```R
> tab
            Actual
Predicted    setosa versicolor virginica
  setosa         50          0         0
  versicolor      0         48         2
  virginica       0          2        48
> 1-sum(diag(tab))/sum(tab)
[1] 0.02666667
```

可以看到模型预测原数据的错误率为3%左右。

---

# 模型调试

调试模型得到最佳参数。

```R
# Tuning 
set.seed(123)
tmodel = tune(svm, Species~.,data = iris,
              ranges = list(epsilon = seq(0,1,0.1),
                            cost = 2^(2:7)))
plot(tmodel)
```

**cost**参数可以修改，直到模型最佳。

![参数选择](https://github.com/lixiang117423/R/raw/master/MachineLearning/SupportVectorMachine/figures/4.png)

颜色越深的参数越佳。可以直接选择模型给出的最佳参数即可。

---

# 最佳模型选择

```R
# best model
mymodel = tmodel$best.model
summary(mymodel)
plot(mymodel, data = iris,
     Petal.Width~Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))
```

![最佳模型](https://github.com/lixiang117423/R/raw/master/MachineLearning/SupportVectorMachine/figures/3.png)

---

# 最佳模型准确度

```R
# Confusion Matrix and Misclassification Error
pred = predict(mymodel, iris)
tab = table(Predicted = pred, Actual = iris$Species)
1-sum(diag(tab))/sum(tab)
```

运行结果如下：

```R
> tab
            Actual
Predicted    setosa versicolor virginica
  setosa         50          0         0
  versicolor      0         48         0
  virginica       0          2        50
> 1-sum(diag(tab))/sum(tab)
[1] 0.01333333
```

校正后的模型的错误率为1.3%左右。明显比最初的模型准确度更改了。

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！