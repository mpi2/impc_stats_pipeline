.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste0(
      '\n ========================================================================',
      '\n ~> If you have any question about this package contact us               ',
      '\n ~> hamedhm@ebi.ac.uk                                                    ',
      '\n ~> This project is supported by European Bioinformatic Institute (EBI)  ',
      '\n ~> https://www.mousephenotype.org/                                      ',
      '\n ========================================================================'
    ),
    domain   = NULL,
    appendLF = TRUE
  )
}
###########
expF = function(x, k, l, m, cdf, direction = c(1, 1)) {
  m  = unique(m)
  r  = cdf(x,
          m - l * direction[1],
          1 / k) *
    (1 - (cdf(x,
              m + l * direction[2],
              1 / k)))
  r  = r / ifelse(!is.null(r) &&
                   max(r, na.rm = TRUE) != 0, max(r, na.rm = TRUE), 1)
  return(r)
}
###########
expWeight = function(t                         ,
                     k                         ,
                     l                         ,
                     m                = 0      ,
                     direction        = c(1, 1),
                     plot             = FALSE  ,
                     zeroCompensation = 0      ,
                     cdf              = plogis ,
                     progress         = FALSE  ,
                     ...) {
  requireNamespace('Rfast')
  m  = unique(m)
  lm = length(m)
  r  = matrix(0, ncol = lm, nrow = length(t))
  if (progress)
    pb = txtProgressBar(
      min = 0  ,
      max = lm,
      style = 3,
      width = 50,
      char = '>'
    )
  ############################################################
  r = sapply(X = 1:lm, function(i) {
    expF(
      x = t                 ,
      k = k                 ,
      l = l                 ,
      m = m[i]              ,
      cdf = cdf             ,
      direction = direction
    )
  })
  ############################################################
  # Begin. Just process the rows that are not zero!
  FinalS = rowsums(r)
  r = r[FinalS > zeroCompensation, , drop = FALSE]
  ############################################################
  s = 0
  for (j in 1:lm) {
    val = sapply(
      X = comb_n(lm,  j, simplify = FALSE), 
      FUN = function(i) {
        return( rowprods(r[, i, drop = FALSE]))
      }
    )
    s = s + (-1) ^ (j + 1) * rowsums(val)
    if (progress)
      setTxtProgressBar(pb, j)
    if ( sum(s) == length(s)) {
      break
    }
  }
  ############################################################
  # End. Just process the rows that are not zero!
  FinalS[FinalS > zeroCompensation] = s
  s = FinalS
  ############################################################
  if (progress) {
    setTxtProgressBar(pb, lm)
    close(pb)
    message('# Zero rows = ',
            sum(FinalS <= zeroCompensation),
            ' [',
            round(sum(FinalS <= zeroCompensation) / length(FinalS) * 100, 1),
            '%]')
  }
  ############################################################
  # Scale the weights to [0-1] on the data INTERVAL
  # Note that the weights are always in [0-1] interval but
  # not necessarily in the domain of the data
  # Do not use the scale() to get only a plain vector of values
  ############################################################
  MaxMin = max(s, na.rm = TRUE) - min(s, na.rm = TRUE)
  if (!is.nan(MaxMin) &&
      !is.na(MaxMin) &&
      MaxMin > sqrt(.Machine$double.eps)) {
    s = (s - min(s, na.rm = TRUE)) / MaxMin
  }
  ############################################################
  # For the problems that the machine epsilon needed
  ############################################################
  s[s <= zeroCompensation] =  zeroCompensation
  ############################################################
  # Plot output
  ############################################################
  if (plot) {
    plot.weight(
      weight = s      ,
      t = t           ,
      m = m           ,
      l = l           ,
      ...
    )
  }
  return(s)
}

# Plot weights (only used in expWeight() function)
plot.weight = function(weight, t, m, l, ...) {
  m  = unique(m)
  lm = length(m)
  plot(
    t,
    weight,
    xlab = 'Time',
    ylab = 'Weight',
    xaxt = 'n',
    xlim = c(min(t, na.rm = TRUE) * .85, max(t, na.rm = TRUE) * 1.5),
    ...
  )
  abline(
    v   = c(m)  ,
    lty = 1     ,
    col = 1:lm  ,
    lwd = 2
  )
  abline(v = c(m + l, m - l)   ,
         lty = 2               ,
         col = 1:lm)
  axis(
    side = 1                        ,
    at =   c(m, m - l, m + l)       ,
    labels = c(m, m - l, m + l)     ,
    lwd.ticks = 2
  )
  legend(
    'topright'                                                               ,
    title = '[m, m+L, m-L]'                                                  ,
    legend = paste(round(m, 3), round(m + l, 3), round(m - l, 3), sep = ', '),
    col = 1:lm                                                               ,
    fill = 1:lm
  )
}
###########
lseq = function (from = 1,
                 to = 5,
                 length.out = 6,
                 adj = 1)
{
  r = exp(seq(log(from) / adj, log(to), length.out = length.out))
  return(r)
}
###########
msg = function(args, ...) {
  message(
    '\n Min observation required: ' ,
    args$min.obs + length(args$m),
    '\n  The number of modes: ',
    length(unique(args$t[args$m])),
    '\n  The time      range: ',
    paste(round(range(args$t), 3), collapse = ' to '),
    '\n  The bandwidth range: ',
    paste(round(range(args$l), 3), collapse = ' to '),
    ' [',
    length(args$l),
    ' split]',
    '\n  The shape     range: ',
    paste(round(range(args$k), 3), collapse = ' to '),
    ' [',
    length(args$k),
    ' split]',
    '\n  Threshold          : ',
    args$threshold,
    '\n',
    if (sum(args$pvalThreshold) <= 0) {
      paste0(
        ' Sensitivity: ',
        '\n\t mean = ',
        args$sensitivity[1],
        '\n\t variance = ',
        args$sensitivity[2],
        '\n\t mean*variance = ',
        args$sensitivity[3],
        '\n\t normality test = ',
        args$sensitivity[4]
      )
    } else{
      paste0(
        ' p-Value threshold: ',
        '\n\t mean = ',
        args$pvalThreshold[1],
        '\n\t variance = ',
        args$pvalThreshold[2],
        '\n\t mean*variance = ',
        args$pvalThreshold[3],
        '\n\t normality test = ',
        args$pvalThreshold[4],
        ###################################
        '\n Sensitivity: ',
        '\n\t mean = ',
        args$sensitivity[1],
        '\n\t variance = ',
        args$sensitivity[2],
        '\n\t mean*variance = ',
        args$sensitivity[3],
        '\n\t normality test = ',
        args$sensitivity[4]
      )
    }
  )
}
###########
checkWeightsN <- function(w             ,
                          t             ,
                          check = 1     ,
                          #0 off, #1 all+NoSingleDate, #2 (1) but SingleDates allowed
                          threshold = 10 ^ -18,
                          normalise = FALSE) {
  n  = length(t)
  if (!is.null(w) && check > 0) {
    r      = 1:n
    zw     = which(abs(w) >= threshold)
    #### need at least 2 observations in a group
    if (check == 1) {
      t      = as.character(t)
      tz     = table(t[zw])
      MorTh1 = names(tz)[which(tz > 1)]
      zw     = zw[t[zw] %in% MorTh1] # no singleDay-singleData
    }
    if (length(zw) > 0) {
      r  = r[zw]
      if (normalise) {
        w  = w[zw] / sum(w[zw])
      } else{
        w  = w[zw]
      }
    } else{
      message(
        '\n ~> Windowing weights are ignored! (1) weights may be all close to zero. (2) there may be only one date in dataset [can cause errors in the mixed model]. (3) Setting check = 1 or check = 2 may solve the problem.\n'
      )
      w = (w * 0 + 1)
      r = 1:n
    }
  } else{
    r = 1:n
  }
  return(list(
    w = w,
    wInd = r,
    NoZeWe = sum(w >= threshold)
  ))
}

shapiro.test0 =function(x){
  x  = na.omit(x)
  lx = length(x) 
  
  if(lx>3 && lx<5000 && is.numeric(x)){
    r = shapiro.test(x)$p.value
  }else{
    r = 1
  }
  return(r)
}
###########
tv.test = function(obj                           ,
                   args                          ,
                   name          = 'parameter'   ,
                   sensitivity   = c(1, 1, 1, 0) ,
                   pvalThreshold = c(0, 0, 0, 0) ,
                   residFun      = function(m) {
                     resid(m)
                   }                             ,
                   predictFun    = function(m) {
                     predict(m)
                   },
                   debug         = FALSE         ,
                   ...) {
  if (!is.null(obj$models) && length(obj$models) > 1) {
    ldata  = nrow(obj$data)
    dfv    = lapply(obj$models, residFun)
    dfp    = lapply(obj$models, predictFun)
    tt     = c()
    for (i in 2:length(obj$models)) {
      oi     = obj$weights[[i]]
      oim    = obj$weights[[i - 1]]
      thresh = 0 #args$threshold*0
      # Variations
      vtl  = var.test(dfv[[i]][oi >= thresh], dfv[[i - 1]][oim >= thresh])$p.value
      vtlp = var.test(dfp[[i]][oi >= thresh], dfp[[i - 1]][oim >= thresh])$p.value
      # Means
      ttl  = t.test  (dfv[[i]][oi >= thresh], dfv[[i - 1]][oim >= thresh])$p.value
      ttlp = t.test  (dfp[[i]][oi >= thresh], dfp[[i - 1]][oim >= thresh])$p.value
      # Normality
      nntp  = shapiro.test0(dfp[[i]][oi >= thresh])
      nntl  = shapiro.test0(dfv[[i]][oi >= thresh])
      
      if (sum(pvalThreshold) <= 0) {
        tt[i - 1] =
          sensitivity[4] * (nntl + nntp) +
          sensitivity[3] * ((vtl * ttl) + (vtlp * ttl) + (vtl * ttlp) + (vtlp * ttlp) +
                              (vtl * vtlp * ttl) + (vtl * vtlp * ttlp) +
                              (vtl * ttl * ttlp) + (vtlp * ttl * ttlp) +
                              (vtl * ttl * vtlp * ttlp)
          ) +
          sensitivity[2] * (vtl + vtlp + vtl * vtlp)  +
          sensitivity[1] * (ttl + ttlp + ttl * ttlp)
      } else{
        f01 = function(x, thresh) {
          (x > thresh) + (x <= thresh) * x
        }
        tt[i - 1] =
          sensitivity[4] * (f01(nntp , pvalThreshold[4])  + 
                            f01(nntl , pvalThreshold[4])) +
          sensitivity[3] * (
            f01(vtl * ttl     , pvalThreshold[3]) +
              f01(vtlp * ttl  , pvalThreshold[3]) +
              f01(vtl * ttlp  , pvalThreshold[3]) +
              f01(vtlp * ttlp , pvalThreshold[3])
          ) +
          sensitivity[2] * (
            f01(vtl , pvalThreshold[2]) + f01(vtlp , pvalThreshold[2]) +
              f01(vtl , pvalThreshold[2]) * f01(vtlp , pvalThreshold[2])
          ) +
          sensitivity[1] * (
            f01(ttl , pvalThreshold[1]) + f01(ttlp , pvalThreshold[1]) +
              f01(ttl , pvalThreshold[1]) * f01(ttlp , pvalThreshold[1])
          )
      }
    }
    ####
    al0 = data.frame(
      tp.pval = c(tt, Inf),
      ol      = obj$output$ObsInInterval,
      l       = obj$output[, name]
    )
    mimi      = min(args$min.obs + length(args$m), ldata, na.rm = TRUE)[1]
    pal       =     al0$ol >= ifelse(name == 'l', mimi, 0)
    pmin      = min(al0$tp.pval[pal], Inf, na.rm = TRUE)
    final     =    (al0$l[pal & (al0$tp.pval %in% pmin)])[1]
    
    if (is.null(final)      ||
        is.null(pmin)       ||
        length(final) < 1   ||
        length(pmin)  < 1   ||
        is.infinite(final)  ||
        is.infinite(pmin)   ||
        is.na(final)        ||
        is.na(pmin)) {
      message(paste('\n An optimal', name, 'is not found. Max value will be used.'))
      final = NULL
    } else{
      message ('\n score: ',
               round(pmin, digits = 4),
               '. ',
               name,
               ': ',
               round(final, 4))
      # ##### Debug
      if (debug) {
        plot(
          al0$l[is.finite(al0$tp.pval)]       ,
          al0$tp.pval[is.finite(al0$tp.pval)] ,
          type = 'b'                          ,
          ylab = 'Residual/Predict test pvals',
          xlab = name,
          main = name
        )
        abline(
          v = final,
          col = 2  ,
          lty = 3  ,
          lwd = 3
        )
        legend(
          'bottomright' ,
          fill = 2      ,
          legend = paste0('Final ', name, ' = ', round(final, 4))
        )
        plot(
          al0$l[is.finite(al0$tp.pval)]      ,
          al0$ol[is.finite(al0$tp.pval)]     ,
          type = 'b'                         ,
          ylab = paste0(
            'Sum(weight), MinSumRequiredIncludingTreatment = ' ,
            round(mimi, 4)
          ),
          xlab = name                        ,
          main = name
        )
        abline(
          v = final                          ,
          h = mimi                           ,
          col = c(2, 3)                      ,
          lty = c(3, 4)                      ,
          lwd = c(1, 3)
        )
        legend('bottomright'                      ,
               fill = 2:3,
               legend = c(
                 paste0('MinSumRequiredInclTr = ', round(mimi, 4)),
                 paste0('final ', name, ' = ', round(final, 4))
               ))
      }
      # #####
    }
  } else{
    final = NULL
  }
  return(list(value = final, score = pmin))
}
###########
is.function0 = function(FUN          ,
                        na.rm = TRUE ,
                        sort  = TRUE ,
                        ...) {
  if (is.function(FUN)) {
    r = FUN()
  } else{
    r = FUN
  }
  if (na.rm)
    r = r[!is.na(r)]
  if (sort)
    r = sort(r, ...)
  return(r)
}


