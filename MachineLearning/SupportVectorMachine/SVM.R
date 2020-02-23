# Data
data("iris")

library(ggplot2)
qplot(Petal.Length,Petal.Width,data = iris,color = Species)

# Support Vector Machine
library(e1071)
mymodel = svm(Species~., data = iris, kernel = 'radial')
summary(mymodel)
plot(mymodel, data = iris,
     Petal.Width~Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

# Confusion Matrix and Misclassification Error
pred = predict(mymodel, iris)
tab = table(Predicted = pred, Actual = iris$Species)
1-sum(diag(tab))/sum(tab)

# Tuning 
set.seed(123)
tmodel = tune(svm, Species~.,data = iris,
              ranges = list(epsilon = seq(0,1,0.1),
                            cost = 2^(2:7)))
plot(tmodel)

# best model
mymodel = tmodel$best.model
summary(mymodel)
plot(mymodel, data = iris,
     Petal.Width~Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

# Confusion Matrix and Misclassification Error
pred = predict(mymodel, iris)
tab = table(Predicted = pred, Actual = iris$Species)
1-sum(diag(tab))/sum(tab)
