VectorOutPutNames = function(clean = TRUE,
                             replace = '_',
                             lower = FALSE,
                             ...) {
  r = c(
    'Method',
    'Dependent variable',
    'Batch included',
    "Batch p-val",
    'Residual variances homogeneity',
    "Residual variances homogeneity p-val",
    'Genotype contribution',
    'Genotype estimate',
    'Genotype standard error',
    'Genotype p-Val',
    'Genotype percentage change',
    'Sex estimate',
    'Sex standard error',
    'Sex p-val',
    'Weight estimate',
    'Weight standard error',
    'Weight p-val',
    'Gp1 genotype',
    'Gp1 Residuals normality test',
    'Gp2 genotype',
    'Gp2 Residuals normality test',
    'Blups test',
    'Rotated residuals normality test',
    'Intercept estimate',
    'Intercept standard error',
    'Interaction included',
    'Interaction p-val',
    'Sex FvKO estimate',
    'Sex FvKO standard error',
    'Sex FvKO p-val',
    'Sex MvKO estimate',
    'Sex MvKO standard error',
    'Sex MvKO p-val',
    'Classification tag',
    'Transformation',
    'Additional information'#,
    #'extraColumns'
  )
  if (clean) {
    r = RemoveSpecialChars(r, replaceBy = replace)
    #r = paste('"', r, '"', sep = '')
  } else{
    r = paste('"', r, '"', sep = '')
  }
  if (lower) {
    r = tolower(r)
  }
  return(r)
}
