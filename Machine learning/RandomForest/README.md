**随机森林**（random forest）是一种有监督的学习方法，可以进行**分类**和**回归**。如果目标变量是**分类变量**，随机森林可以进行**分类**；如果目标变量是连续变量，随机森林可以进行**回归预测**。

**随机森林**有如下优势：

- 分类准确率更高；
- 能够有效处理高维数据；
- 擅长处理大数据；
- 可以用于大量缺失值的数据；
- 在分类的同时度量变量对分类的相对重要性。

---

本文数据来自**UCI**的**Cardiotocography Data Set**（[点击下载](https://github.com/lixiang117423/R/blob/master/Machine learning/RandomForest/CTG.csv))。一共有2126个观测值，22个观测变量，其中最后一个变量**NSP**作为目标变量。

---

1. 准备数据，并将变量**NSP**改成分类变量：

   ```R
   data <- read.csv("CTG.csv", header = TRUE) # 读取数据
   str(data)
   data$NSP <- as.factor(data$NSP) # 将变量NSP设置成因子
   table(data$NSP)
   ```

2. 将数据划分成训练集和测试集：

   ```R
   set.seed(123) # 设置随机数种子，保证结果的课重复性
   ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3)) # 数据随机采样设置
   train <- data[ind==1,] # 70%数据用作训练集
   test <- data[ind==2,] # 30%数据用作测试集
   ```

3. 构建随机森林模型
   其中，参数**ntree**和**mtry**分别根据后面的步骤**6**和**7**确定。

   ```R
   library(randomForest) # 加载包
   set.seed(222) # 设置随机数种子，保证结果的课重复性
   rf <- randomForest(NSP~., # NSP作为目标变量
                      data=train, # 使用训练集构建随机森林
                      ntree = 300, # 决策树的数量，默认的是500
                      mtry = 8, # 每个分组中随机的变量数，一般是变量数开根
                      importance = TRUE, # 是否评估预测变量的重要性
                      proximity = TRUE) # 是否计算行之间的接近度
   ```

4. 查看构建的随机森林模型和有贡献的变量

   ```R
   > print(rf)
   
   Call:
    randomForest(formula = NSP ~ ., data = train, ntree = 300, mtry = 8,      importance = TRUE, proximity = TRUE) 
                  Type of random forest: classification
                        Number of trees: 300
   No. of variables tried at each split: 8
   
           OOB estimate of  error rate: 5.35%
   Confusion matrix:
        1   2   3 class.error
   1 1139  15   6  0.01810345
   2   48 164   1  0.23004695
   3    7   3 112  0.08196721
   
   > attributes(rf)
   $names
    [1] "call"            "type"            "predicted"      
    [4] "err.rate"        "confusion"       "votes"          
    [7] "oob.times"       "classes"         "importance"     
   [10] "importanceSD"    "localImportance" "proximity"      
   [13] "ntree"           "mtry"            "forest"         
   [16] "y"               "test"            "inbag"          
   [19] "terms"          
   
   $class
   [1] "randomForest.formula" "randomForest"  
   ```

   可以看出是进行**分类**（classification）。还可以看出目标变量下每个分类值的**预测准确率**以及“有用”的变量是哪些。

5. 使用模型进行预测

   - 使用训练好的模型预测训练集

     ```R
     > library(caret)
     > p1 <- predict(rf, train)
     > confusionMatrix(p1, train$NSP)
     Confusion Matrix and Statistics
     
               Reference
     Prediction    1    2    3
              1 1159    1    0
              2    1  212    0
              3    0    0  122
     
     Overall Statistics
                                               
                    Accuracy : 0.9987          
                      95% CI : (0.9952, 0.9998)
         No Information Rate : 0.7759          
         P-Value [Acc > NIR] : < 2.2e-16       
                                               
                       Kappa : 0.9964          
                                               
      Mcnemar's Test P-Value : NA              
     
     Statistics by Class:
     
                          Class: 1 Class: 2 Class: 3
     Sensitivity            0.9991   0.9953  1.00000
     Specificity            0.9970   0.9992  1.00000
     Pos Pred Value         0.9991   0.9953  1.00000
     Neg Pred Value         0.9970   0.9992  1.00000
     Prevalence             0.7759   0.1425  0.08161
     Detection Rate         0.7753   0.1418  0.08161
     Detection Prevalence   0.7759   0.1425  0.08161
     Balanced Accuracy      0.9981   0.9973  1.00000
     ```

     

   可以看出准确率达到**99.87%**。

   - 使用训练好的模型预测测试集

     ```R
     > p2 <- predict(rf, test)
     > confusionMatrix(p2, test$NSP)
     Confusion Matrix and Statistics
     
               Reference
     Prediction   1   2   3
              1 479  17   3
              2  14  61   2
              3   2   4  49
     
     Overall Statistics
                                               
                    Accuracy : 0.9334          
                      95% CI : (0.9111, 0.9516)
         No Information Rate : 0.7845          
         P-Value [Acc > NIR] : <2e-16          
                                               
                       Kappa : 0.8132          
                                               
      Mcnemar's Test P-Value : 0.7633          
     
     Statistics by Class:
     
                          Class: 1 Class: 2 Class: 3
     Sensitivity            0.9677  0.74390  0.90741
     Specificity            0.8529  0.97086  0.98960
     Pos Pred Value         0.9599  0.79221  0.89091
     Neg Pred Value         0.8788  0.96209  0.99132
     Prevalence             0.7845  0.12995  0.08558
     Detection Rate         0.7591  0.09667  0.07765
     Detection Prevalence   0.7908  0.12203  0.08716
     Balanced Accuracy      0.9103  0.85738  0.94850
     ```

     准确率为**93.34%**。
     **需要注意的是用模型得到的错误率是5.35%，但是用模型去预测训练集的准确率是99.87%，模型预测测试集的准确率是93.34%**。之所以预测训练集那么高是以为我们的模型就是用训练集得到的。会发现模型预测测试集适的准确率和模型给出的准确率差不多。

6. 绘制随机森林模型错误率

   ```R
   plot(rf)
   ```

   ![随机森林模型错误率](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/1.png)

   从图上可以看出，“trees”在300以后，“Error”基本就是稳定的了，所以选择**ntree = 300**是OK的。

7. 调试随机森林模型以获得最佳**mtry**

   ```R
   t <- tuneRF(train[,-22], 
               train[,22],
          stepFactor = 0.5,
          plot = TRUE,
          ntreeTry = 300,
          trace = TRUE,
          improve = 0.05)
   ```

   ![mtry参数确定](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/6.png)

   选择图中位置最低点**8**作为**mtry**参数。

8. 直方图展示树的节点数

   ![树的节点数](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/2.png)

   这个一般只要前面的**mtree**和**mtry**没问题，就没问题。

9. 观测变量的重要性
   查看模型中观测变量的重要性，相当于每个变量在模型中的影响力。

   ```R
   varImpPlot(rf,
              sort = T,
              n.var = 10,
              main = "Top 10 - Variable Importance")
   ```

   ![变量重要性](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/3.png)

   越往上的变量贡献越大。
   还可以用别的方式查看变量的重要性：

   ```R
   importance(rf)
   varUsed(rf)
   ```

10. 每个变量对每个预测值的局部影响力
    下面的代码查看模型中变量**ASTV**对目标变量**NSP为2**的影响：

    ```R
    partialPlot(rf, train, ASTV, "2")
    ```

    ![变量局部重要性](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/4.png)

11. 从模型中提取树

    ```R
    getTree(rf, 1, labelVar = TRUE)
    ```

12. 接近矩阵的多维标度图

    ```R
    MDSplot(rf, train$NSP)
    ```

    ![多维标度图](https://github.com/lixiang117423/R/raw/master/Machine%20learning/RandomForest/figures/5.png)

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！