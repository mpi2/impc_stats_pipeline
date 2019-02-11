SmoothWin = function(object                                           ,
                     data                                             ,
                     t                                                ,
                     m                                                ,
                     l = function(ignore.me.in.default) {
                       r = SmoothWin:::lseq(
                         from = 1                                          ,
                         to = max(abs(t[m] - min(t, na.rm = TRUE)), abs(t[m] - max(t, na.rm = TRUE)), 1),
                         length.out = min(500, max(1, diff(range(
                           t,na.rm = TRUE
                         ))))
                       )
                       r = unique(round(r))
                       return(r)
                     },
                     k = SmoothWin:::lseq(from = .5                   ,
                                          to = 10                     ,
                                          length.out = 50)            ,
                     min.obs   = function(ignore.me.in.default) {
                       lutm = length(unique(t[m]))
                       r = ifelse(lutm > 1, 35, max(pi * sqrt(length(t)), 35))
                       r = max(r * lutm, length(m), na.rm = TRUE)
                       r = min(r       , length(t), na.rm = TRUE)
                       return(r)
                     }                                                ,
                     direction = c(1, 1)                              ,
                     weightFUN = function(x) {
                       x
                     }                                                ,
                     residFun = function(x) {
                       resid(x)
                     }                                                ,
                     predictFun = function(x) {
                       predict(x)
                     }                                                ,
                     weightORthreshold = 'weight'                     ,
                     cdf = plogis                                     ,
                     check = 2                                        ,
                     sensitivity   = c(1, 1, 1, 0)                    ,
                     pvalThreshold = c(0, 0, 0, 0)                    ,
                     threshold     = sqrt(.Machine$double.eps) * 10   ,
                     zeroCompensation = 0                             ,
                     messages      = TRUE                             ,
                     seed          = NULL                             ,
                     simple.output = FALSE                            ,
                     debug         = FALSE                            ,
                     ...) {
  if (!is.null(seed))
    set.seed(seed)
  min.obs = ceiling(is.function0(min.obs))
  l = is.function0(l, decreasing = FALSE)
  k = is.function0(k, decreasing = TRUE)
  argg    = c(as.list(environment()), list())
  if (length(unique(t[m])) > 15)
    message('More than 15 modes detected. The entire procedure can take a longer than usual!')
  
  if (length(m) > min.obs) {
    stop('`min.obs` is less than the total number of treatments!')
  } else {
    msg(argg)
  }
  ### 1. Determining l
  message('\n 1|3 Searching for l ...\n')
  rl = gridSearchModel(
    object = object                       ,
    data = data                           ,
    weightFUN = weightFUN                 ,
    check = check                         ,
    t = t                                 ,
    m = t[m]                              ,
    l = l                                 ,
    k = max(k)                            ,
    threshold = threshold                 ,
    messages = messages                   ,
    onlyOne  = FALSE                      ,
    cdf      = cdf                        ,
    zeroCompensation = zeroCompensation   ,
    weightOrthreshold = weightORthreshold ,
    direction = direction                 , 
    ...
  )
  finall = tv.test(
    obj = rl                              ,
    args = argg                           ,
    name = 'l'                            ,
    residFun = residFun                   ,
    sensitivity = sensitivity             ,
    pvalThreshold = pvalThreshold         ,
    predictFun = predictFun               ,
    debug = debug
  )
  if (is.null(finall$value)){
    finall$value = max(l)
  }
    
  ### 2. Determining k
  message('\n 2|3 Searching for k ...\n')
  rk = gridSearchModel(
    object = object                       ,
    data = data                           ,
    weightFUN = weightFUN                 ,
    check = check                         ,
    t = t                                 ,
    m = t[m]                              ,
    l = finall$value                      ,
    k = k                                 ,
    threshold = threshold                 ,
    messages = messages                   ,
    onlyOne  = FALSE                      ,
    cdf      = cdf                        ,
    zeroCompensation = zeroCompensation   ,
    weightOrthreshold = weightORthreshold ,
    direction = direction                 ,
    ...
  )
  finalk = tv.test(
    obj = rk                              ,
    args = argg                           ,
    name = 'k'                            ,
    residFun = residFun                   ,
    predictFun = predictFun               ,
    sensitivity = sensitivity             ,
    pvalThreshold = pvalThreshold         ,
    debug = debug
  )
  if (is.null(finalk$value))
    finalk$value = max(k)
  
  ##### final model
  message('\n 3|3 Forming the final model ...\n')
  finalr = gridSearchModel(
    object = object                       ,
    data = data                           ,
    weightFUN = weightFUN                 ,
    check = check                         ,
    t = t                                 ,
    m = t[m]                              ,
    l = finall$value                      ,
    k = finalk$value                      ,
    threshold = threshold                 ,
    messages = messages                   ,
    onlyOne  = TRUE                       ,
    cdf      = cdf                        ,
    zeroCompensation = zeroCompensation   ,
    weightOrthreshold = weightORthreshold ,
    direction = direction                 ,
    ...
  )
  if (simple.output) {
    rk = rl = NULL
  }
  out = list(
    object  = object                    ,
    data    = data                      ,
    final.k = finalk                    ,
    final.l = finall                    ,
    finalModel = finalr                 ,
    model.l = rl                        ,
    model.k = rk                        ,
    min.obs = min.obs                   ,
    input   = argg
  )
  class(out) = 'SmoothWin'
  return(out)
}



# Plot windowing object
plot.SmoothWin = function(x, ylab = 'Response', col = NULL ,   ...) {
  if (!is.null(x$finalModel$models)) {
    t  = x$input$t
    y  = x$data[, all.vars(formula(x$finalModel$models))[1]]
    m  = x$input$m
    ly = length(y)
    if (is.unsorted(t, na.rm = TRUE))
      message('To get the right plot, make sure that the dataset is sorted on time!')
    
    if (is.null(col)) {
      col = rgb(
        abs(.8 - x$finalModel$FullWeight),
        1 - x$finalModel$FullWeight      ,
        1 - x$finalModel$FullWeight
      )
    }
    plot(
      t,
      y,
      xlab = 'Time',
      ylab = ylab,
      sub = paste(
        'l='                                                      ,
        round(x$final.l$value, 2)                                       ,
        ', k='                                                    ,
        round(x$final.k$value, 2)                                       ,
        ', '                                                      ,
        ifelse(x$input$weightORthreshold == 'weight', 'ASS=', '#'),
        round(x$finalModel$output$ObsInInterval, 3)               ,
        ' [~'                                                     ,
        round(x$finalModel$output$ObsInInterval / ly * 100)       ,
        '%]'       ,
        ', MaxBW=' ,
        max(x$input$l, na.rm = TRUE),
        ', MinObs=',
        x$min.obs  ,
        '+'        ,
        length(m)  ,
        sep = ''
      ),
      col = col,
      ...
    )
    
    abline(v = unique(t[m]),
           lty = 2,
           col = 'gray')
    wp = x$finalModel$FullWeight
    lines(
      t,
      min(y) + wp * (max(y) - min(y)),
      col = 'gray'                   ,
      lty = 4                        ,
      lwd = 4
    )
    return(invisible(list(
      weight = wp   ,
      k = x$final.k ,
      l = x$final.k ,
      object = x
    )))
  } else{
    message('Windowing failed. No plot available for the failed models.')
  }
  #####
}