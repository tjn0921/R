**决策树**（Decision tree）理解成通过一系列的条件判断下一步应该是什么，得到什么样的结果。就像下面这个图：

![决策树例子](https://github.com/lixiang117423/R/raw/master/Machine%20learning/DecisionTree/figures/2.png)

---

本文示例数据来自UCL机器学习中心。下载地址：https://github.com/lixiang117423/R/blob/master/Machine%20learning/DecisionTree/Cardiotocographic.csv

---

# 数据处理

导入数据并随机将数据80%分成训练集，20%分成测试集。

```R
#导入数据
mydata = read.csv("Cardiotocographic.csv")
mydata$NSPF = as.factor(mydata$NSP)
str(mydata)

#划分训练集和测试集
set.seed(1234)
pd = sample(2, nrow(mydata), replace = T, prob = c(.8,.2))

train = mydata[pd == 1,]
test = mydata[pd == 2,]
```

---

# 构建决策树

```R
tree = ctree(NSPF~LB+AC+FM, # 选择其中几个变量进行建树
             data = train,
             controls = ctree_control(mincriterion = 0.99,
                                      minsplit = 500))
plot(tree)
```

![决策树](https://github.com/lixiang117423/R/raw/master/Machine%20learning/DecisionTree/figures/3.png)

从图上看出，决策树从**AC**这个变量开始，依次判断，到最后就看哪个柱子高就是哪个了。

---

# 测试模型

## 用构建的模型测试训练集

```R
tab = table(predict(tree),train$NSPF)
print(tab)
1-sum(diag(tab)/sum(tab))
```

可以看到如下结果：

```
> print(tab)
   
       1    2    3
  1 1222   70  112
  2  126  156   32
  3    0    0    0
  
> 1-sum(diag(tab)/sum(tab))
[1] 0.1979045
```

可以看出用模型预测训练集时，模型的错误率大概是20%。

----

## 用构建的模型预测测试集

```R
testpred = predict(tree,newdata = test)
tab1 = table(testpred,test$NSPF)
tab1
1-sum(diag(tab1)/sum(tab1))
```

可以看到如下结果：

```
> tab1
        
testpred   1   2   3
       1 274  21  28
       2  33  48   4
       3   0   0   0
       
> 1-sum(diag(tab1)/sum(tab1))
[1] 0.2107843
```

 可以看出用模型预测测试集时，模型的错误率大概是21%。

---

# 其他展示树的方法

```R
library(rpart)
tree1 = rpart(NSPF~LB+AC+FM, train)

library(rpart.plot)
rpart.plot(tree1,type = 0)
```

![其他可视化方法](https://github.com/lixiang117423/R/raw/master/Machine%20learning/DecisionTree/figures/4.png)

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！