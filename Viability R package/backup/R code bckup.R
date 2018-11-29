cat('\014')
rm(list = ls(all = TRUE))
set.seed(123456)
# Loading data
#download.file(url = 'https://www.ebi.ac.uk/~hamedhm/Fully_loaded_R_data.RData','Fully loaded R data.RData')
#download.file(url = 'https://www.ebi.ac.uk/~hamedhm/All_Human_Cells_data.RData','All Human Cells data.Rdata')
load(file = 'All Human Cells data.Rdata')
load(file = 'Fully loaded R data.RData')
overall_human_cell_lines = overall_human_cell_lines[, -c( #Violeta remove two columns
  grep(
    pattern = c('gecko'),
    x = colnames(overall_human_cell_lines),
    ignore.case = TRUE
  ),
  grep(
    pattern = c('shrna'),
    x = colnames(overall_human_cell_lines),
    ignore.case = TRUE
  )
)]
####
# Blomen_HAP1_Q_ortho = Blomen_HAP1_Q_ortho[!duplicated(Blomen_HAP1_Q_ortho$`MGI id`), ]
# Blomen_Q_KBM7_ortho = Blomen_Q_KBM7_ortho[!duplicated(Blomen_Q_KBM7_ortho$`MGI id`), ]
# Hart_BF_ortho = Hart_BF_ortho[!duplicated(Hart_BF_ortho$`MGI id`), ]
# Wang_CS_ortho = Wang_CS_ortho[!duplicated(Wang_CS_ortho$`MGI id`), ]
#
# m1 = merge(Blomen_HAP1_Q_ortho,
# 					 Blomen_Q_KBM7_ortho,
# 					 by = 'MGI id',
# 					 all = TRUE)
# m2 = merge(m1, Hart_BF_ortho, by = 'MGI id', all = TRUE)
# m3 = merge(m2, Wang_CS_ortho, by = 'MGI id', all = TRUE)
# m4 = m3[, -grep('chr', colnames(m3))]
#

####
library(GeneViability)
aM = initialization(data = DMDD_test,
                    response = 'Viability',
                    speed = TRUE)
bM = tranSamp(object = aM,
              percent = 1,
              speed = TRUE)
#cM = pcsne(object = bM, speed = TRUE)
###
aH = initialization(data = overall_human_cell_lines,
                    response = NULL,
                    speed = TRUE)
bH = tranSamp(
  object = aH,
  percent = 1,
  FUN = scale,
  speed = TRUE
)
#cH = pcsne(object = bH, speed = TRUE)
###




f0 = merge(
  data.frame(
    #viability = cM$inputdata$viability, #HAMED
    'MGI id' = rownames(bM$rdata),
    # HAMED
    bM$rdata  ,
    check.names = FALSE
  ),
  data.frame(
    #HAMED
    bH$rdata,
    'MGI id' = rownames(bH$rdata),
    check.names = FALSE
  ),
  by = 'MGI id',
  all = FALSE
)
f0 = merge(
  data.frame(
    'MGI id' = DMDD_train$`MGI id`,
    viability = DMDD_train$Viability,
    check.names = FALSE
  ),
  f0,
  by = 'MGI id',
  all = TRUE
)
final_data = droplevels(f0[!f0$viability %in% c('Conflicting','Subviable'),])
final_data = final_data[complete.cases(final_data[, -c(1:2)]), ]

a=rownames(final_data) = final_data[, 1]
final_data = final_data[, -1]
final_data$viability = as.character(final_data$viability)
final_data$viability[is.na(final_data$viability)] = 'unknown'
final_data$viability = as.factor(final_data$viability)
f00 = final_data

################





dM = pcsne(object = list(
  rdata = final_data,
  response = 'viability'),
  speed = TRUE,
  dmsne = 3
)

cor(dM$rdata)

# kmeans
cl =50
km1 = kmeans(dM$rdata,centers = cl,iter.max = 10^5)
km2 = kmeans(f00[,-1],centers = cl,iter.max = 10^5)
rgl.open()
scatter3d(
  x = dM$rdata$V1,
  y = dM$rdata$V2,
  z = dM$rdata$V3,
  #groups = dM$inputdata$viability,
  point.col = km1$cluster,
  surface = FALSE,
  pch='+'
)

newdata1 = newdata2 = f00
newdata1$viability =as.character(newdata1$viability)
newdata2$viability =as.character(newdata2$viability)

thresh = 20
for(i in 1:cl) {
  t1 = table(final_data$viability[km1$cluster == i])
  t2 = table(final_data$viability[km2$cluster == i])
  #print(t1)
  vl1 = c(t1[1] / t1[3], t1[3] / t1[1])
  vl2 = c(t2[1] / t2[3], t2[3] / t2[1])
  if (max(vl1) > thresh) {
    lbl1 = names(vl1)[which.max(vl1)]
    newdata1$viability[km1$cluster == i &
                         newdata1$viability %in% 'unknown'] = paste(lbl1,'_p',sep = '')
  }
  if (max(vl2) > thresh) {
    lbl2 = names(vl2)[which.max(vl2)]
    newdata2$viability[km2$cluster == i &
                         newdata2$viability %in% 'unknown'] = paste(lbl2,'_p',sep = '')
  }
}
table(newdata1$viability)
table(newdata2$viability)
write.csv(x = newdata1 , file = 'd:tsne_kmeans.csv')
write.csv(x = newdata1 , file = 'd:org_kmeans.csv')


# convert the
newdata1$viability[newdata1$viability %in% 'Lethal_p'] ='Lethal'
newdata1$viability[newdata1$viability %in% 'Viable_p'] ='Viable'
table(newdata1$viability)

###########################
final_data = data.frame(viability = newdata1$viability, dM$rdata)
f2 = f22 = final_data[!is.na(final_data$viability), ]
scatter3d(
  x = f2$V1,
  y = f2$V2,
  z = f2$V3,
  groups = f2$viability,
  surface = FALSE
)
######################### NN
f2$viability[f2$viability %in% 'unknown'] = NA
f2 = f20 = droplevels(data.frame(f2, check.names = TRUE))



f2$viability = as.integer(f2$viability) - 1
hist(f2$viability)

library(neuralnet)
resultTrain = resultTest = list()
for(j in 3){
  set.seed(123456)
  for (i in 1:1) {
    #print(i)
    ### 1000 times bootstrap
    smp = sample(1:nrow(f2),
                 size = nrow(f2) *2/3,
                 replace = FALSE)
    csmp = (1:nrow(f2))[-smp]


    nn <- tryCatch(
      expr = {
        neuralnet(
          formula = as.formula(paste(
            'viability', paste(colnames(f2)[-1], collapse  = '+'), sep  = '~'
          )),
          f2[smp,],
          hidden = c(j),
          linear.output = FALSE
        )
      },
      warning = function(war) {
        warning('\n => Warning : ', war, '\n')
        return(NULL)
      },
      error = function(err) {
        warning('\n => error : ', err, '\n')
        return(NULL)
      }
    )
    if (is.null(nn))
      next

    ### On training
    a = compute(nn, f2[smp, -1])
    a$net.result <- sapply(a$net.result, round, digits = 0)
    resultTrain[[i]] = table(a$net.result, f2[smp, 1])

    ### On test
    b = compute(nn, f2[-smp, -1])
    b$net.result <-
      sapply(b$net.result, function(x) {
        #x[x > 2] = 2
        #x [x < 1] = 1
        round(x, digits = 0)
      })
    resultTest[[i]] = table(b$net.result, f2[-smp, 1])
    #print(table(b$net.result, f2[-smp, 1]))
  }
  #######
  l = list(test=resultTest,train=resultTrain)
  s1 = s2 = table(c(0,1),c(1,0))*0
  counter1=counter2=0
  for (i in 1:(length(l$test))) {
    if (!all(is.null(l$test[[i]])) && all(dim(l$test[[i]]) == c(2, 2))) {
      counter1=counter1+1
      s1 = s1 + l$test[[i]]
      #print(s1)
      #print(counter1)
    }
    if (!all(is.null(l$train[[i]])) && all(dim(l$train[[i]]) == c(2, 2))) {
      counter2=counter2+1
      s2 = s2 + l$train[[i]]
      #print(s2)
      # print(counter2)
    }
    print(i)
  }
  print(j)
  print(prop.table(s1/counter1,margin = 2))
  print(prop.table(s2/counter2,margin = 2))
}






##########
#classif.nnet, classif.LiblineaRL1L2SVC, classif.lssvm||,classif.lvq1|,classif.mda,classif.mlp,classif.multinom,classif.nnet|classif.nnTrain,classif.OneR,classif.PART,classif.qda,classif.quaDA,classif.randomForest|||classif.randomForestSRC|||
l1 = l2 = list()
f20 = f20[!is.na(f20$viability),]
set.seed(123456)
library(mlr)#classif.quaDA *classif.nnTrain
for (i in 1:50) {
  smp = sample(1:nrow(f20),
               size = nrow(f20) / 3 * 2,
               replace = FALSE)
  csmp = (1:nrow(f20))[-smp]

  task1 = makeClassifTask(data = data.frame(f20, check.names = T),
                          target = "viability")
  task2 = makeClassifTask(data = data.frame(f20, check.names = T), target = "viability")
  lrn = makeLearner('classif.nnTrain')

  model1 = train(lrn, task1, subset = smp)
  model2 = train(lrn, task2, subset = smp)

  pred1 = predict(model1, task = task1, subset = csmp)
  pred2 = predict(model2, task = task2, subset = csmp)

  #performance(pred1, measures = list(mmce, acc))
  l1[[i]] = calculateConfusionMatrix(pred1, relative = TRUE, sums = FALSE)$relative.row
  l2[[i]] = calculateConfusionMatrix(pred2, relative = TRUE, sums = FALSE)$relative.row
  print(i)
  #plotLearnerPrediction(lrn, task = task)
}

###


