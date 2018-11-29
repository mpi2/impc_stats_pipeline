suspectRemove = function(replica  = 2500,
                         data,
                         unknown = 'unknown',
                         train.prop = .5,
                         nree = 150,
                         formula = viability ~ V1 + V2 + V3 + V4 + V5,
                         thresh = .5,
                         trim.u = .15,
                         trim.d = .07,
                         seed = 123456,
                         ...) {
  argg <- c(as.list(environment()), list())
  library(randomForestSRC)
  set.seed(123456)
  nx = data[which(!data$viability %in% unknown),]
  n   =  nrow(nx)
  m   =  matrix(ncol = replica, nrow = nrow(data))
  rownames(m) = rownames(data)
  #####
  for (i in 1:replica) {
    cat('\r', round(i / replica * 100), '%')
    #####
    x = sample(1:n, size = train.prop * n)
    train = nx[x , ]
    test  = nx[-x, ]
    #####
    rf = rfsrc(formula = formula,
               data = droplevels(train),
               ntree = nree,
               ...)
    #####
    pr = predict(rf, newdata = data[, -which(names(data) %in% all.vars(formula)[1])])
    m[, i] = pr$predicted[, 1]
  }
  r = apply(m, 1, function(x) {
    sum(x > thresh)
  })

  dow = (replica / 2 * (1 - trim.d))
  up  = (replica / 2 * (1 + trim.u))
  thresh.e = (r > dow & r < up)
  hist(r)
  abline (
    v = c(dow, up),
    lty = 3,
    lwd = 3,
    col = 2
  )

  plot(prcomp(nx[, -which(names(nx) %in% all.vars(formula)[1])])$x[, 1:2], col = thresh.e + 1, pch = thresh.e + 1)
  legend(
    'top',
    legend = c('Normal', 'Suspicious '),
    pch = 1:2,
    col = 1:2
  )
  err  = data[thresh.e,]
  #err2 = subset(err, subset = !all.vars(formula)[1] %in% unknown)
  return(list(
    argg = argg,
    data = data,
    outMatrix = m,
    all.sus   = err#,
    #known.sus = err2
  ))
}
