rm(list = ls())

data = read.csv('vehicle.csv')
head(data)
pairs(data[3:5])

# multiple linear regression
res = lm(lc ~ Mileage + lh, data = data)
summary(res)

reduced = lm(lc~lh, data = data)
summary(reduced)

anova(res, reduced)

# prediction(default 95% confidence)
predict(reduced, data.frame(lh = 10), interval = 'confidence')
