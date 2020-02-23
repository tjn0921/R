#Read data file
mydata = read.csv("Cardiotocographic.csv")
mydata$NSPF = as.factor(mydata$NSP)
str(mydata)

#Partition data into Traaining data and test data
set.seed(1234)
pd = sample(2, nrow(mydata), replace = T, prob = c(.8,.2))

train = mydata[pd == 1,]
test = mydata[pd == 2,]

# Decision tree with Party
tree = ctree(NSPF~LB+AC+FM, data = train,
             controls = ctree_control(mincriterion = 0.99,
                                      minsplit = 500))
plot(tree)

# Predict
predict(tree, test, type = 'prob')

# Prediction
predict(tree1, test)

# Misclassification error for 'train' data
tab = table(predict(tree),train$NSPF)
print(tab)
1-sum(diag(tab)/sum(tab))

# Misclassification error with test data
testpred = predict(tree,newdata = test)
tab1 = table(testpred,test$NSPF)
tab1
1-sum(diag(tab1)/sum(tab1))

# Decision tree with rpart
library(rpart)
tree1 = rpart(NSPF~LB+AC+FM, train)

library(rpart.plot)
rpart.plot(tree1,type = 0)