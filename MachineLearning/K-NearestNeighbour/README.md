**K-Nearest Neighbour**是一种极为常见的分类算法。原理很简单：设定一个**K**值，这个范围内谁多，那么这个新的东西就属于谁。就像下面这个图，如果我们设定**K=黑色实线**，那么<font color=green>绿色的○</font>就属于<font color=red>红色△</font>这个分组。同理，如果我们设定**K=黑色虚线**，那么<font color=green>绿色的○</font>就属于<font color=blue>蓝色△</font>这个分组。

<img src="https://github.com/lixiang117423/R/raw/master/MachineLearning/K-NearestNeighbour/figures/1.png" alt="KNN算法示意图" style="zoom: 33%;" />

---

本教程使用数据是根据gre、GPA和rank决定是否录取一个学生。下载地址：https://github.com/lixiang117423/R/blob/master/MachineLearning/K-NearestNeighbour/binary.csv

---

# 数据预处理

导入数据，将数据划分成训练集和测试集。

```R
# Libraries 加载包
library(caret)
library(pROC)
library(mlbench)

# Example-1 Student Applications (Classification)
data <- read.csv('binary.csv', header = T)
str(data)
# 修改数据
data$admit[data$admit == 0] <- 'No'
data$admit[data$admit == 1] <- 'Yes'
data$admit <- factor(data$admit)

# Data Partition 划分数据为训练集和测试集
set.seed(1234)
ind <- sample(2, nrow(data), replace = T, prob = c(0.7, 0.3))
training <- data[ind == 1,]
test <- data[ind == 2,]
```

---

# 构建模型

```R
# KNN Model
trControl <- trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 3,
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary)
set.seed(222)
fit <- train(admit ~ .,
             data = training,
             method = 'knn',
             tuneLength = 20,
             trControl = trControl,
             preProc = c("center", "scale"),# 数据标准化
             metric = "ROC",
             tuneGrid = expand.grid(k = 1:60)) # 设定K值的范围
```

---

# 查看最佳**K**值：

```R
> fit
k-Nearest Neighbors 

284 samples
  3 predictor
  2 classes: 'No', 'Yes' 

Pre-processing: centered (3), scaled (3) 
Resampling: Cross-Validated (10 fold, repeated 3 times) 
Summary of sample sizes: 256, 256, 256, 256, 255, 256, ... 
Resampling results across tuning parameters:

  k   ROC        Sens       Spec      
   1  0.5412281  0.7100000  0.36962963
   2  0.5631871  0.7019298  0.35703704
   3  0.5822368  0.7973684  0.34074074
   4  0.5622904  0.7770175  0.26111111
   5  0.5878558  0.8065789  0.28518519
   6  0.5911501  0.8241228  0.27740741
   7  0.5884016  0.8585965  0.28296296
   8  0.5892982  0.8550877  0.26851852
   9  0.5965010  0.8725439  0.29074074
  10  0.5899123  0.8657895  0.27444444
  11  0.5955945  0.8799123  0.28555556
  12  0.5876170  0.8695614  0.27740741
  13  0.5942300  0.8659649  0.24185185
  14  0.5949854  0.8885965  0.25703704
  15  0.5987768  0.8816667  0.22777778
  16  0.6125097  0.8974561  0.22074074
  17  0.6257943  0.8976316  0.23629630
  18  0.6256481  0.9008772  0.21481481
  19  0.6286647  0.8976316  0.22851852
  20  0.6372612  0.9027193  0.22185185
  21  0.6386209  0.9064035  0.21074074
  22  0.6393519  0.9113158  0.22518519
  23  0.6397710  0.9200000  0.21444444
  24  0.6438743  0.9268421  0.21740741
  25  0.6460039  0.9215789  0.19962963
  26  0.6482359  0.9252632  0.19925926
  27  0.6554288  0.9287719  0.20296296
  28  0.6596881  0.9356140  0.21037037
  29  0.6670419  0.9356140  0.19925926
  30  0.6726901  0.9426316  0.19888889
  31  0.6659649  0.9460526  0.17740741
  32  0.6657797  0.9478070  0.18148148
  33  0.6623830  0.9495614  0.18851852
  34  0.6627973  0.9461404  0.19259259
  35  0.6652242  0.9495614  0.18851852
  36  0.6679678  0.9442982  0.18148148
  37  0.6682261  0.9477193  0.18148148
  38  0.6704825  0.9426316  0.18555556
  39  0.6682846  0.9459649  0.18925926
  40  0.6617788  0.9442982  0.18185185
  41  0.6624708  0.9408772  0.17481481
  42  0.6600097  0.9443860  0.15703704
  43  0.6633138  0.9478070  0.16037037
  44  0.6617008  0.9495614  0.16000000
  45  0.6627924  0.9478070  0.14962963
  46  0.6605897  0.9478070  0.15037037
  47  0.6598294  0.9494737  0.14629630
  48  0.6578558  0.9512281  0.14962963
  49  0.6581481  0.9511404  0.14259259
  50  0.6615692  0.9546491  0.14592593
  51  0.6638060  0.9546491  0.13148148
  52  0.6637037  0.9581579  0.14222222
  53  0.6655507  0.9546491  0.12814815
  54  0.6666715  0.9528947  0.11703704
  55  0.6662671  0.9597368  0.11370370
  56  0.6639133  0.9614912  0.11777778
  57  0.6619737  0.9597368  0.10333333
  58  0.6612037  0.9632456  0.10333333
  59  0.6640156  0.9597368  0.09592593
  60  0.6681189  0.9597368  0.09962963

ROC was used to select the optimal model using the largest value.
The final value used for the model was k = 30.
```

可以看出模型给出的最佳**K**值是**30**。

----

# 绘图展示K值

```R
plot(fit)
```

<img src="https://github.com/lixiang117423/R/raw/master/MachineLearning/K-NearestNeighbour/figures/2.png" alt="K线" style="zoom: 100%;" />

可以看到，在**K=30**时，**ROC**值最大。所以**K=30**是最佳的选择。

---

# 查看使用的变量在模型中的重要性

```R
> varImp(fit)
ROC curve variable importance

     Importance
gpa      100.00
rank      25.18
gre        0.00
```

可以看出，**gre**这个变量在模型中几乎没有任何作用。

---

# 使用模型预测测试集

```R
> pred <- predict(fit, newdata = test)
> confusionMatrix(pred, test$admit)
Confusion Matrix and Statistics

          Reference
Prediction No Yes
       No  78  30
       Yes  4   4
                                          
               Accuracy : 0.7069          
                 95% CI : (0.6152, 0.7877)
    No Information Rate : 0.7069          
    P-Value [Acc > NIR] : 0.5461          
                                          
                  Kappa : 0.0887          
                                          
 Mcnemar's Test P-Value : 1.807e-05       
                                          
            Sensitivity : 0.9512          
            Specificity : 0.1176          
         Pos Pred Value : 0.7222          
         Neg Pred Value : 0.5000          
             Prevalence : 0.7069          
         Detection Rate : 0.6724          
   Detection Prevalence : 0.9310          
      Balanced Accuracy : 0.5344          
                                          
       'Positive' Class : No 
```

可以看出模型的准确率是**70.69%**左右。

---

# 其他算法的KNN

```R
# Example-2 Boston Housing (Regression)
data("BostonHousing")
data <- BostonHousing 
str(data)

# Data Partition
set.seed(1234)
ind <- sample(2, nrow(data), replace = T, prob = c(0.7, 0.3))
training <- data[ind == 1,]
test <- data[ind == 2,]

# KNN Model
trControl <- trainControl(method = 'repeatedcv',
                          number = 10,
                          repeats = 3)
set.seed(333)
fit <- train(medv ~.,
             data = training,
             tuneGrid = expand.grid(k=1:70),
             method = 'knn',
             metric = 'Rsquared',
             trControl = trControl,
             preProc = c('center', 'scale'))

# Model Performance
fit
plot(fit)
varImp(fit)
pred <- predict(fit, newdata = test)
RMSE(pred, test$medv)
plot(pred ~ test$medv)
```

---

**致谢：**
感谢YouTube博主*Bharatendra Rai*博士的视频！
