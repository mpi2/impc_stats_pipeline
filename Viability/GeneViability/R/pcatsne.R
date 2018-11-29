pcsne = function(object,
                 umap = TRUE,
                 supervised = FALSE,
                 perplexity = 30,
                 dmpca = NULL,
                 dmsne = 3,
                 dmumap = 3,
                 seed = 123456,
                 speed = TRUE,
                 theta = .5,
                 unknown = 'unknown',
                 response = 'viability',
                 ...) {
  library(uwot)
  argg <- c(as.list(environment()), list())
  ipar = par()
  rdata = object$rdata
  num.cols = sapply(rdata, is.numeric)
  rdata.back = rdata

  if (!umap) {
    for (i in perplexity) {
      rdata = rdata[!duplicated(rdata[, num.cols, drop = FALSE]), ]
      if (!identical(rdata, rdata.back))
        message('\n duplicated data removed \n')
      if (is.null(dmpca))
        dmpca = min(ceiling(2 * sqrt(sum(num.cols))), ceiling(sum(num.cols) /	2))

      if (!is.null(dmpca) && dmpca > 0) {
        set.seed(seed)
        r5 = r55 = prcomp(rdata[, num.cols], scale. = FALSE, center = FALSE)
        r5 = as.data.frame(r55$x[, 1:dmpca])
      } else{
        r5 = r55 = rdata
      }
      set.seed(seed)
      tsne = Rtsne::Rtsne(
        X = r5,
        dims = dmsne,
        perplexity =  i,
        check_duplicates = TRUE,
        pca = FALSE,
        theta = theta
      )
      r5 = as.data.frame(tsne$Y)
      rownames(r5) = rownames(rdata)
      message('\nPCA dimension : ', dmpca, ', SNE dimension : ', dmsne, '\n')
      if (dmsne > 2) {
        library(car)
        library(rgl)
        if (dmsne > 3)
          message('\n Only the first 3 columns are ploted! \n')
        rgl.open()
        scatter3d(
          x = (r5$V1),
          y = (r5$V2),
          z = (r5$V3),
          surface = FALSE,
          group = if (!is.null(object$response)) {
            rdata[, 1]
          } else{
            NULL
          },
          point.col = ifelse(
            !is.null(object$response),
            as.integer(rdata[, 1]),
            'black'
          ),
          fit = "smooth",
          grid = FALSE,
          ellipsoid = FALSE,
          ...
        )
      } else{
        plot(
          r5$V1,
          r5$V2,
          col = ifelse(!is.null(object$response), as.integer(rdata[, 1]), 'red'),
          main = i,
          pch = ifelse(!is.null(object$response), as.integer(rdata[, 1]), 1),
          ...
        )
        abline(lm(V2 ~ V1, data = r5))
        abline(lm(V1 ~ V2, data = r5))
        ########
      }
      cat('\n perplexity = ', i, '\n')
    }

    if (!speed) {
      if (dmpca > 0) {
        plot(r55, main = 'PCA')
        print(summary(r55))
        print(cor(as.matrix(r5)))
      }
      if (!is.null(object$response)) {
        # 1.3 other plots
        if (dmpca > 0) {
          par(mfrow = c(2, 2))
          for (i in 1:(dmpca - 1)) {
            for (j in (i + 1):dmpca) {
              plot(
                r55$x[, i],
                r55$x[, j],
                col = as.integer(rdata[, 1]),
                pch = as.integer(rdata[, 1]),
                xlab = colnames(r55$x)[i],
                ylab = colnames(r55$x)[j],
                main = 'PCA'
              )
              legend(
                'topleft',
                legend = paste(unique(rdata[, 1])),
                pch = as.integer(unique(rdata[, 1])),
                col = as.integer(unique(rdata[, 1])),
                bg = 'white'
              )
            }
          }
        }
        if (dmsne > 0) {
          par(mfrow = c(2, 2))
          for (i in 1:(dmsne - 1)) {
            for (j in (i + 1):dmsne) {
              plot(
                r5[, i],
                r5[, j],
                col = as.integer(rdata[, 1]),
                pch = as.integer(rdata[, 1]),
                xlab = colnames(r5)[i],
                ylab = colnames(r5)[j],
                #cex=.1,
                main = 'SNE'
              )
              legend(
                'topleft',
                legend = paste(unique(rdata[, 1])),
                pch = as.integer(unique(rdata[, 1])),
                col = as.integer(unique(rdata[, 1])),
                bg = 'white'
              )
            }
          }
        }
      }
    }
    umTest = NULL
  } else{
    set.seed(seed)
    library(parallel)
    message('Umap in progress ....')
    r55 = NULL
    viaCol = which(names(rdata) %in% response)
    if (supervised) {
      # Test
      train = droplevels(rdata[!rdata[, viaCol] %in% unknown,])
      message('Stp 1. training ....')
      umTest = umap(
        X = train[, num.cols, drop = FALSE],
        n_components = dmumap,
        verbose = TRUE,
        y =  train[, viaCol],
        n_threads = detectCores(),
        ret_model = TRUE,
        ...
      )
      message('Stp 2. prediction ....')
      r5 = umap_transform(X = rdata[, -viaCol],
                          model = umTest,
                          #n_components = dmumap,
                          verbose = TRUE
                          #n_threads = detectCores()
                          )
    } else{
      r5 = umap(
        X = rdata[, num.cols, drop = FALSE],
        n_components = dmumap,
        verbose = TRUE,
        n_threads = detectCores(),
        ...
      )
      umTest = NULL
    }
    r5 = data.frame(r5)
    names(r5) = paste0('V', 1:ncol(r5))
    rownames(r5) = rownames(rdata)
  }
  gc()
  par(mfrow = ipar$mfrow)
  return(
    list(
      pca = r55,
      rdata = r5,
      inputdata = rdata,
      response = object$response,
      arg = argg,
      dmpca = dmpca,
      dmsne = dmsne,
      dmumap = dmumap,
      umTest = umTest
    )
  )

}