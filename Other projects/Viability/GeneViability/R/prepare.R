initialization = function(data,
                          seed = 123456,
                          unique = 'MGI id',
                          response = 'viability' ,
                          shuffle = FALSE,
                          removeDuplicates = TRUE,
                          missingRemove = TRUE,
                          impute = TRUE,
                          max.miss.thresh = .6,
                          filter = c('Lethal', 'Viable'),
                          k = 5,
                          speed = TRUE) {
  argg <- c(as.list(environment()), list())
  if (!is.null(seed))
    set.seed(seed = seed)
  data = as.data.frame(data) # to round tibble
  if (is.null(response) || !response %in% colnames(data)) {
    message('\n Response is not found and then set to NULL\n')
    response = NULL
    #discriptive = FALSE
  }

  message('\n\nInput dataset : ', paste(
    c('rows = ', 'columns = '),
    dim(data),
    collapse = ' | ',
    sep = ''
  ))
  message('\n\nVariable names (',
          length(colnames(data)),
          '): \n',
          paste(colnames(data), collapse = '|'))


  if (!is.null(response))
    plot(data[, response], col = unique(data[, response]), main = 'Original data')

  rdata.back = as.data.frame(data)
  if (removeDuplicates) {
    dpl = duplicated(data[, unique])
    rdata = data[!dpl,]
    if (any(dim(rdata) != dim(data))){
      message(
        '\n\n\n Duplicated data found and removed! [',
        dim(data)[1] - dim(rdata)[1],
        '=',
        round((dim(data)[1] - dim(rdata)[1]) / dim(rdata)[1] * 100, digits = 2),
        '%] check: ',
        paste(data[, unique][dpl], collapse = '|'),
        '\n'
      )}else{
        message('No duplicate found!')
      }
  } else{
    rdata = data
    dpl = NULL
  }

  #rearrangne dataset
  MGIid = which(colnames(rdata) %in% c(unique))
  viability = which(colnames(rdata) %in% c(response))
  Others_col = which(!colnames(rdata) %in% c(response, unique))
  base::rownames(rdata) = rdata[, MGIid]


  if (!is.null(response)) {
    rdata = rdata[, c(viability, Others_col),,drop=FALSE]
    if (!all(is.null(filter))) {
      rdata = subset(rdata, rdata[, 1] %in% filter)
    }
  } else{
    rdata = rdata[, c(Others_col),drop=FALSE]
    message('\n\n\n Response is not included!\n')
  }

  if (shuffle) {
    message('\n\n\n Shuffling in process ... \n')
    rdata = rdata[sample(nrow(rdata)),,drop=FALSE ]
  }

  ####### missings
  # Missing handling [remove the rows with more than xXx missings]
  num.cols = sapply(rdata, is.numeric)
  miss = apply(rdata[, num.cols, drop = FALSE], 1, function(x) {
    (sum(is.na(x)) / length(x)) > max.miss.thresh
  })
  if (!missingRemove) {
    miss = miss & FALSE
  }
  if (sum(miss) > 0 && missingRemove) {
    message(
      '\n\n\n Missing data more than ',
      max.miss.thresh * 100,
      ' percent in rows ',
      '(',
      sum(miss),
      '=',
      round(sum(miss) / length(miss) * 100, digits = 0),
      '%) [removed]: ',
      paste(rownames(rdata)[which(miss == TRUE)], collapse = '|'),
      ' \n\r\n\n'
    )
    rdata = droplevels(rdata[!miss, ])
  }

  # Imputing other missings using kmean
  if (impute && ncol(rdata[, num.cols, drop = FALSE]) > 1) {
    library(impute)
    if (sum(miss) > 0)
      message('\n\n\n [if possible] Missing imputation in progress ... \n','Threshold: ',max.miss.thresh)
    nomissRdata = impute.knn(
      as.matrix(rdata[, num.cols, drop = FALSE]),
      rng.seed = seed,
      rowmax = max.miss.thresh,
      k = k
    )
  } else {
    message('\n\n\n [if possible] missing data removed!\n')
    nomissRdata = rdata[, num.cols, drop = FALSE]
    nomissRdata$data = rdata[, num.cols, drop = FALSE]
  }


  if (!is.null(response) && removeDuplicates) {
    message('Duplicates:\n')
    print(table(data[dpl, response]))
    rdata  = data.frame(viability = rdata[, 1],
                        nomissRdata$data,
                        check.names = FALSE)
  } else{
    rdata = droplevels(rdata[!miss,,drop=FALSE])
    rdata  = data.frame(nomissRdata$data,
                        check.names = FALSE)
  }

  rdata = droplevels(rdata)
  if (!speed & !is.null(response)) {
    slices = table(rdata[, 1])
    print(slices)
    pct = round(slices / sum(slices) * 100)
    lbls = paste(names(slices), ': ', pct , '% [', slices, ']', sep = '') # add percents to labels
    #lbls = paste(lbls , sep = "") # ad % to labels
    pie(x = slices,
        labels = lbls)
  }

  message('\n\n\n Final dataset : ', paste(
    c('rows = ', 'columns = '),
    dim(rdata),
    collapse = ' | ',
    sep = ''
  ))
  message(
    '\n\n\n Variable names : \n length = ',
    length(colnames(rdata)),
    '\n\n\n list: ',
    paste(colnames(rdata), collapse = '|')
  )
  if (!is.null(response))
    plot(rdata[, response], col = unique(rdata[, response]), main = 'Output data')
  gc()
  return(
    list(
      rdata = rdata,
      response = response,
      inputdata = rdata.back,
      duplicates = dpl,
      missings.g.thresh = miss,
      args = argg
    )
  )
}