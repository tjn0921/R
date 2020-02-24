# data
data(iris)
str(iris)

# partition data
set.seed(111)
ind = sample(2, nrow(iris),replace = T, prob = c(.8,.2))
train = iris[ind == 1,]
test = iris[ind == 2,]

# scatter plot and correlations
library(psych)
pairs.panels(train[,-5],
             gap = 0,
             bg = c('red','yellow','blue')[train$Species],
             pch = 21)

# PCA
pc = prcomp(train[,-5],
            center = T,
            scale. = T)
pc$center
pc$scale
print(pc)

# orthogonality of PCs
pairs.panels(pc$x,
             gap = 0,
             bg = c('red','yellow','blue')[train$Species],
             pch = 21)

# bi-plot
library(ggbiplot)
g = ggbiplot(pc,
             obs.scale = 1,
             var.scale = 1,
             groups = train$Species,
             ellipse = T,
             circle = F,
             ellipse.prob = 0.68)
g = g + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal',
                                                legend.position = 'top')
print(g)

# prediction with Principal Components
trg = predict(pc, train)
trg = data.frame(trg, train$Species)

tst = predict(pc,test)
tst = data.frame(tst, test$Species)

# multinomial logistic regression with first two PCs
library(nnet)
trg$species = relevel(trg$train.Species, ref = 'setosa')
mymodel = multinom(train.Species~PC1+PC2,data = trg)
summary(mymodel)

# Confusion matrix
p = predict(mymodel, tst)
tab = table(p, test$Species)
tab
1-sum(diag(tab))/sum(tab)
