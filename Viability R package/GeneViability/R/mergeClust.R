mergeClust = function(rdata,
                      response = 'viability',
                      cl = 75,
                      thresh = 20,
                      unknown.lev = 'unknown',
                      maxiter = 10 ^ 5,
                      nstart  = 50,
                      suffix = '_p',
                      seed = 123456,
                      minin = 20,
                      method = 'hkmeans',
                      ...) {
  argg <- c(as.list(environment()), list())
  ### kmeans

  if (!response %in% colnames(rdata))
    stop('\nResponse does not found in the data!\n')

  if ((length(levels(rdata[, response])) != 3)  |
      !(unknown.lev %in%  levels(rdata[, response])))
    stop('\nlevels of response different from 3 or ',
         unknown.lev,
         ' not found!')
  message('\n',
          ifelse(method == 'hkmeans', '[h]', ''),
          'kmeans in progress ....\n')
  set.seed(seed)
  if (method == 'hkmeans') {
    km = factoextra::hkmeans(rdata[, which(!colnames(rdata) %in% response), drop =
                                     FALSE],
                             k = cl,
                             iter.max = maxiter)
    message('*info : nstart is not defined for hkmeans method. ')
  } else{
    km = kmeans(rdata[, which(!colnames(rdata) %in% response), drop = FALSE],
                centers = cl,
                iter.max = maxiter,
                nstart = nstart)
  }

  newdata  = rdata
  newdata[, response] = as.character(newdata[, response])
  mergecl = c()

  for (i in 1:cl) {
    t1 = table(rdata[, response][km$cluster == i])
    clo = (1:3)[-which(names(t1) %in% unknown.lev)]
    vl1 = c(t1[clo[1]] /( t1[clo[2]]), t1[clo[2]] / (t1[clo[1]]))
    #print((t1[clo[1]] + t1[clo[2]]))
    #print(max(vl1, na.rm = TRUE))
    if ((t1[clo[1]] + t1[clo[2]]) >minin  && (max(vl1, na.rm = TRUE) > thresh))
         {
      mergecl = c(mergecl, i)
      lbl = names(vl1)[which.max(vl1)]
      newdata[, response][km$cluster == i &
                            newdata[, response] %in% unknown.lev] = paste(lbl, suffix, sep = '')
      message(
        'N:',t1[clo[1]] + t1[clo[2]],
        ', Ratio: ',
        round(max(vl1, na.rm = TRUE), 1),
        '; Cluster ',
        i,
        ' ',
        unknown.lev,
        ' assigned to -> ',
        lbl,
        '; Details: ',
        paste(names(t1), t1, collapse = '; ', sep = ':')
      )
    }

  }
  print(table(newdata[, response]))
  coli = rep(1, cl)
  coli [mergecl] = 2
  plot(
    table(km$cluster),
    ylab = 'The number of points in cluster [including unknown]',
    xlab = 'Cluster',
    las = 3,
    col = coli,
    cex = .5,
    ...
  )
  legend(
    x = 'top',
    legend = c('Not merged', 'Merged'),
    horiz = TRUE,
    inset = -.05,
    xpd = TRUE,
    fill = c(1, 2),
    bg = 'white'
  )

  message(
    '\n Total clusters: ',
    cl,
    '. Min points in clusters [including unknown]: ',
    min(table(km$cluster)),
    '. Min information required in each cluster : ',
    minin,
    ' Points.'
  )

  out =   data.frame(
    'clusters' = km$cluster,
    'merged.clu' = (km$cluster %in% mergecl),
    'assigned.viability' = newdata$viability,
    rdata,
    check.names = FALSE
  )
  rownames(out) = rownames(rdata)
  return(list(out=out,args=argg))
}