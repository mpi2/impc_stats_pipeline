tranSamp = function(object,
                    percent = 1,
                    FUN = function(x) {
                      log(log(x + 1))
                    },
                    speed = TRUE,
                    seed = 123456,
                    ...) {
  #######
  argg <- c(as.list(environment()), list())
  set.seed(seed)

  r2 = object$rdata
  nr2 = nrow(r2)
  smpT  = sample(seq(nr2),
                 size = round(nr2 * percent),
                 replace = FALSE)
  num.cols = sapply(r2, is.numeric)

  if (sum(!num.cols) > 0) {
    r3    = droplevels(r2[smpT, ,drop=FALSE])
  } else{
    r3    = r2[smpT, , drop = FALSE]
  }
  message(
    '\nTotal observations: ',
    nrow(r2),
    '. Selected samples: ',
    length(smpT),
    ' [Shuffled] \n '
  )

  # transformations
  num.cols = sapply(r3, is.numeric)
  r3[, num.cols] = FUN(r3[, num.cols,drop=FALSE])

  if (!speed) {
    plot(
      1,
      type = 'n',
      ylim = c(min(r3[, num.cols]), max(r3[, num.cols])),
      xlim = c(0, sum(num.cols))
    )

    if (!is.null(object$response)) {
      num.cols2 = num.cols
      num.cols2[1] = TRUE
      r4 = r3
      r4[, 1] = as.integer(r4[, 1])
      apply(r4, 1, function(y) {
        lines(spline(y[-1]), lty = 2, col = y[1])
      })
      legend(
        'topright',
        legend = unique(r3[, 1]),
        #pch = as.integer(unique(r4[, 1])),
        fill = as.integer(unique(r3[, 1])),
        bg = 'white'
      )
      ## Table
      print(table(r3[, 1]))
    } else{
      apply(r3[, num.cols], 1, function(y) {
        lines(spline(y), lty = 2)
      })
    }
    # heatmap?
    heatmap(
      x = as.matrix(r3[, num.cols]),
      scale = 'row',
      # verbose = TRUE,
      cexRow = 0.5,
      cexCol = .7
    )
  }
  gc()
  return(list(
    rdata = r3,
    response = object$response,
    smp = smpT,
    arg = argg
  ))
}