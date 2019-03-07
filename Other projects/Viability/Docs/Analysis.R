rm(list = ls(all = TRUE))
library(GeneViability)
load('Fully loaded R data.RData')
a = initialization(data = DMDD_train, response = 'Viability')
b = tranSamp(object = a, percent = .3)
c = pcsne(object = b)
#################
## 1.2 Head and Tail
clr = 5
r5 = c$rdata
for (i in 1:(c$dmsne - 1)) {
  for (j in (i + 1):c$dmsne) {
    km = kmeans(r5[, c(i, j)], clr) # HAMED 10
    plot(
      r5[, c(i, j)],
      col = c$inputdata[, 1],
      xlab = paste('V', i),
      ylab = paste('V', j)
    )
    for (k in 1:clr) {
      points(
        x = km$centers[k, 1],
        y = km$centers[k, 2],
        pch = '+',
        col = 'Blue',
        cex = .1
      )
      text(
        x = km$centers[k, 1],
        y = km$centers[k, 2],
        labels = k ,
        col = 'Blue',
        cex = 2
      )
    }
  }
}

set.seed(123456)
km = kmeans(r5, clr)
tbofClust = sapply(1:clr, function(i) {
  (table(c$inputdata[, 1][km$cluster %in% i]))
})
clOfint = c(which.max(tbofClust[1, ] / tbofClust[2, ]),
            which.max(tbofClust[2, ] / tbofClust[1, ]))
tbofClust
clOfint
as.table(tbofClust[, c(clOfint)])
prop.table(as.table(tbofClust[, c(clOfint)]), margin = 2)

scatter3d(
  x = (r5$V1),
  y = (r5$V2),
  z = (r5$V3),
  surface = FALSE,
  point.col =   (km$cluster %in% clOfint) * 1 + 1,
  #km$cluster,
  fit = "smooth",
  grid = FALSE,
  ellipsoid = FALSE
)


# cluster of the tail
km2data = r5[km$cluster %in% clOfint[1], ]
km2orgData = c$inputdata[km$cluster %in% clOfint[1], ]
heatmap(as.matrix(km2data))
plot(km2data$V1, km2data$V2, col = km2orgData[, 1])
set.seed(123456)
km2 = kmeans(km2data, centers = 2)
plot(km2data$V1, km2data$V2, col = km2orgData[, 1], cex = 1.5)
points(km2data$V1,
       km2data$V2,
       col = km2$cluster,
       pch = km2$cluster)
table(km2orgData$Viability)

library(MASS)
fit = qda(
  Viability ~ .,
  data = data.frame(Viability = km2orgData[, 1], km2data),
  na.action = "na.omit",
  CV = FALSE
)
plot(predict(fit)$class)
table(km2orgData[, 1], predict(fit)$class)
table(km2orgData[, 1])
prop.table(table(km2orgData[, 1], predict(fit)$class), margin = 2)

scatter3d(
  x = (km2data$V1),
  y = (km2data$V2),
  z = (km2data$V3),
  surface = FALSE,
  point.col =   km2$cluster,
  #km$cluster,
  fit = "smooth",
  grid = FALSE,
  ellipsoid = FALSE
)
#write.csv(r4.1[km$cluster %in% clOfint[1], ], file = 'd:/head.csv')
table(r4.1$Viability[km$cluster %in% clOfint[1]])
#write.csv(r4.1[km$cluster %in% clOfint[2], ], file = 'd:/tail.csv')
table(r4.1$Viability[km$cluster %in% clOfint[2]])


## Prediction
aTest = initialization(data = DMDD_test, response = 'Viability')
bTest = tranSamp(object = aTest, percent = .3,speed=FALSE)
cTest = pcsne(object = bTest,speed = FALSE)

library(clue)
set.seed(123456)
predictkm = cl_predict(object = km, newdata = cTest$rdata)
tbofClust
table(predictkm)
clOfint
g1 = which(predictkm %in% clOfint[1])
g2 = which(predictkm %in% clOfint[2])


r4.2te = data.frame(cTest$inputdata, 'MGI id' = rownames(cTest$inputdata), check.names = FALSE)
r4.2   = data.frame(c$inputdata  , 'MGI id' = rownames(c$inputdata), check.names = FALSE)
r5.2te = data.frame(cTest$rdata, 'MGI id' = rownames(cTest$rdata), check.names = FALSE)

r5.2   = data.frame(
  via = c$inputdata[,1]  ,
  'MGI id' = rownames(c$inputdata),
  check.names = FALSE
)

m0 = merge(y = r5.2,
           x = r5.2te[g1,],
           by = 'MGI id',
           all = TRUE)

m1 = merge(y = r4.2,
           x = r4.2te[g1,],
           by = 'MGI id',
           all = TRUE)

m2 = merge(x = r4.2te[g2,],
           y = r4.2,
           by = 'MGI id',
           all = TRUE)

table(m1$viability[!is.na(m1[,2])])
table(m2$viability[!is.na(m2[,2])])

## 2nd cluster
pred2 = cl_predict(object = km2, newdata = cTest$rdata[g1,])
g3 = which(pred2 %in% 1)
table(pred2)
m3 = merge(y = r4.2,
           x = r4.2te[g3,],
           by = 'MGI id',
           all = TRUE)
table(m3$viability[!is.na(m3$`4-somite stage.x`)])


### using qdr
library(MASS)
qdapr = predict(fit, cTest$rdata[g1,])
plot(qdapr$class)
table(qdapr$class)
# Discriminat analysis

col = as.integer(m0$via)
col[is.na(col)] = 8
plot(m0$V1, m0$V2, col = col)

##
#write.csv(x = m1,
          file = 'd:predict_lethal.csv')
#write.csv(x = m2,
          file = 'd:predict_viable.csv')



#save(r5te,r4te,r4.1te,file = 'd:r5_r4_r41.rdata')






###############################################################
######################### initial 5
set.seed(123456)
smp = sample(1:length(r5$via),
             size = length(r5$via) / 3 * 2,
             replace = FALSE)
######################### 5
# SVM
library(e1071)


set.seed(123456)
wts = (1 - table(r5$via) / sum(table(r5$via)))
#wts[1] = .001
wts = wts / sum(wts)
svob = tune(
  svm,
  reformulate(termlabels = colnames(r5)[1:dm2], response = colnames(r5)[dm2 +
                                                                          1]),
  data = r5[smp, ],
  kernel = "radial",
  scale = FALSE,
  # HAMED
  probability = FALSE,
  class.weights = wts,
  ranges = list(
    cost =  seq(10 ^ -2, 3, length.out = 5),
    gamma = seq(10 ^ -2, 3, length.out = 5)
  ),
  decision.values = TRUE,
  cachesize = 500
)


plot(svob)
summary(svob)
svob$best.parameters

b = svm(
  formula = reformulate(termlabels = colnames(r5)[1:dm2], response = colnames(r5)[dm2 +
                                                                                    1]),
  data = r5[smp, ],
  scale = FALSE,
  # HAMED
  kernel = 'radial',
  gamma = svob$best.parameters$gamma,
  cost =  svob$best.parameters$cost,
  probability = FALSE,
  class.weights = wts,
  decision.values = TRUE,
  cachesize = 500
)


table(r5$via[smp])
prop.table(table(fitted(b), r5$via[smp]),margin = 2)
prop.table(table(predict(b, newdata = r5[-smp, c(1:dm2)]), r5$via[-smp]),margin = 2)

if (dm2 > 2) {
  match = function(x) {
    apply(x, 1, function(y) {
      y = as.character(y)
      if (y[1] == y[2]) {
        r = (y[1])
      } else {
        r = 1
      }
      return(r)
    })
  }

  mtch = as.numeric(match(cbind(predict(b, newdata = r5[-smp, c(1:dm2)]), r5$via[-smp])))
  table(mtch)
  s3d = scatter3d(
    x = (r5[-smp, 1]),
    y = (r5[-smp, 2]),
    z = (r5[-smp, 3]),
    surface = FALSE,
    fit = "smooth",
    #grid = TRUE,
    ellipsoid = FALSE,
    #text.col = mtch,
    #labels = mtch,
    groups =  predict(b, newdata = r5[-smp, c(1:dm2)])
  )
} else{
  plot(r5$V1[-smp], r5$V2[-smp], col = fitted(b))
  points(r5$V1[-smp],
         r5$V2[-smp],
         col = as.integer(r5$via[-smp]) + 1,
         pch = as.integer(r5$via[-smp]) + 4)
}



######################### 5-2 NN
require(mlbench)
require(mxnet)
set.seed(123456)
nn = mx.mlp(
  as.matrix(r5[smp, 1:(dm2)]),
  r5[smp, dm2+1],
  hidden_node = c(180),
  out_node = 2,
  out_activation = "softmax",
  num.round = 1000,
  #array.batch.size = 15,
  learning.rate = 1*10^-6,
  momentum = 10^-5,
  eval.metric = mx.metric.accuracy
)
#graph.viz(nn$symbol)
predsTR = predict(nn, as.matrix(r5[smp, 1:dm2]))
pred.labelTR = max.col(t(predsTR)) #-1
#prop.table(table(pred.labelTR, as.integer(r5$via[smp])),margin = 2)
prop.table(table(r5$via[smp],c('Viable','Lethal')[pred.labelTR]),margin = 2)
# TEST
preds = predict(nn, as.matrix(r5[-smp, 1:dm2]))
pred.label = max.col(t(preds)) #-- 1
#prop.table(table(r5$via[-smp],pred.label),margin = 2)
prop.table(table(r5$via[-smp],c('Viable','Lethal')[pred.label]),margin = 2)

######################### 6

# f1 = as.formula(paste('viaI~', paste(
# 	'V', 1:dm2, sep = '', collapse =  '+'
# )))
#
# library(neuralnet)
# nn = neuralnet(
# 	formula = f1,
# 	data = r6[smp,],
# 	hidden = c(15, 5,1),
# 	linear.output = TRUE
# 	#stepmax = 10 ^ 6,
# )
# #plot(nn)
# r7 = compute(nn, r6[-smp, 1:dm2])$net.result
# results <- data.frame(actual = r6$viaI[-smp], prediction = r7)
# roundedresults <- sapply(results, round, digits = 1)
# roundedresultsdf = data.frame(roundedresults)
# table(roundedresultsdf$actual, roundedresultsdf$prediction)


##########
library(mlr)
set.seed(123456)
#r6$via[r6$via != 'Lethal'] = 'Viable'
r6 = droplevels(r5)


task = makeClassifTask(data = r6[, -c(dm2 + 2)], target = "via")
lrn = makeLearner('classif.randomForestSRC')#classif.nnet, classif.LiblineaRL1L2SVC, classif.lssvm||,classif.lvq1|,classif.mda,classif.mlp,classif.multinom,classif.nnet|classif.nnTrain,classif.OneR,classif.PART,classif.qda,classif.quaDA,classif.randomForest|||classif.randomForestSRC|||

n = nrow(r6)
train.set = sample(n, size = 2 / 3 * n)
test.set = setdiff(1:n, train.set)
model = train(lrn, task, subset = train.set)
pred = predict(model, task = task, subset = test.set)
performance(pred, measures = list(mmce, acc))
plotLearnerPrediction(lrn, task = task)
calculateConfusionMatrix(pred, relative = TRUE, sums = FALSE)#$relative.row
