df = read.delim(file = '~/../Downloads/all/AllDR13SimV1Phenstat.tsvv/AllDR13SimV1Phenstat.tsv', header = FALSE, sep = '\t')

names = c(
  "Analysis date",
  "status",
  "procedure_group",
  "procedure_stable_id",
  "procedure_name",
  "parameter_stable_id",
  "parameter_name",
  "phenotyping_center",
  "allele_symbol",
  "allele_name",
  "allele_accession_id",
  "gene_symbol",
  "pipeline_name",
  "pipeline_stable_id",
  "strain_accession_id",
  "metadata_group",
  "zygosity",
  "colony_id"
)

# numericCols = grepl('p-value',x = names(df),fixed = TRUE)
# df[,numericCols] = apply(df[,numericCols],2,as.numeric)
# names(df)[numericCols]
#
#
#


# numericCols = grepl('q-value',x = names(df),fixed = TRUE)
# df[,numericCols] = apply(df[,numericCols],2,as.numeric)
# names(df)[numericCols]

###############=================
alpha = 10 ^ -4
df0 = df
names(df0)[1:length(names)]=names
f1 = 'NGenotype contribution p-value'
f11 = 'NGenotype p-value'
f2 = 'WGenotype Contribution p-value'
f22 = 'WGenotype p-value'
names(df0)[20] = f1
names(df0)[21] = f11
names(df0)[22] = f2
names(df0)[23] = f22
names(df0)[4] = 'procedure_stable_id'
names(df0)[8] = 'phenotyping_center'
df0$Centre = df0$phenotyping_center
names(df0)[6] = 'parameter_stable_id'
names(df0)[7] = 'procedure_name'
names(df0)[3] = 'procedure_group'
df0 = df0[complete.cases(df0[, f1], df0[, f2]), ]

dfdiff = df0[(df0[, f2] < alpha &
                df0[, f1] > alpha) |
               (df0[, f2] > alpha & df0[, f1] < alpha), ]
loss = df0[(df0[, f2] < alpha & df0[, f1] > alpha), ]
# write.csv(head(loss,1000),file='hamed.csv',row.names = FALSE)
write.csv(dfdiff, file = 'PhenStatDiffBetweenWinAndNormal.csv', row.names = FALSE)

# nna = apply(df0, 1, function(x) {
#   sum(!is.na(x))
# })
df0 = droplevels(subset(
  df0,
  grepl(
    pattern = 'IMPC_',
    x = df0$procedure_group,
    fixed = TRUE
  )
))

# Stp 1. Summary statistics
sort(table(df0$procedure_group), decreasing = TRUE)
length(unique(df0$procedure_group))
length(unique(df0$phenotyping_center))
dim(df0)

sum(na.omit(df0[, f1] < alpha))
sum(na.omit(df0[, f2] < alpha))


sum(df0[, f1] < alpha) - sum(df0[, f2] < alpha)
round((sum(df0[, f1] < alpha) - sum(df0[, f2] < alpha)) /
        sum(df0[, f1] < alpha) * 100)



library(reshape)

# df1  = df0[, c(f1, f2)]
# df11 = melt(data = df1,value.name = 'P.value')
# df11$variable = gsub(pattern = '_',replacement = ' ',x = df11$variable)
#



df0$Pvalue_normal = df0[, f1] #df0[,f1]
df0$Pvalue_window = df0[, f2] #df0[,f2]

df0 = subset(df0, !df0$procedure_group %in% 'BCMLA_XRY_001')

SimSummary  = function(name, df0, alpha) {
  library(plyr)
  sdtable = as.data.frame(
    ddply(
      df0,
      c(name),
      summarize,
      # 'procedure_name' = unique(procedure_name),
      Procedure        = paste(unique(procedure_group), collapse = ', '),
      n                = length(Pvalue_normal),
      NormalFDR        = sum(Pvalue_normal < alpha, na.rm = TRUE) ,
      WindowedFDR      = sum(Pvalue_window < alpha, na.rm = TRUE) ,
      NormalAlpha      = round(NormalFDR    / n, 4),
      WindowedAlpha    = round(WindowedFDR / n, 4),
      WindowedFDRPercent  = round(WindowedFDR / (WindowedFDR +
                                                   NormalFDR) * 100, digits = 2),
      diff             = (WindowedFDR - NormalFDR),
      Gain             = (WindowedFDR < NormalFDR),
      odd              = round(WindowedFDR / NormalFDR, 2),
      .drop = TRUE
    )
  )


  sdtable = sdtable[order(abs(sdtable$diff), decreasing = TRUE),]
  sdtable = sdtable[sdtable$diff != 0,]
  return(sdtable)
}

hist(df0[, f1][df0[, f1] < alpha], col = 2)
hist(df0[, f2][df0[, f2] < alpha],
     add = TRUE,
     col = rgb(0, 1, 0, alpha = .5))



h0 = SimSummary(name = 'parameter_stable_id', df0 = df0)
h = head(SimSummary(name = 'procedure_group', df0 = df0), 100)
h
round(sum((h$diff[h$diff < 0])) / sum(abs(h$diff)) * 100, 2)
round(sum(abs(h$diff[h$diff > 0])) / sum(abs(h$diff)) * 100, 2)

#
#
# SimSummary(name = c('procedure_group', 'Centre'), df0 = df0)
# SimSummary(name = c('Procedure_group'), df0 = df0)
# SimSummary(name = c('Parameter_stable_id'), df0 = df0)
#
# TopLoss = SimSummary(name = c('Parameter_stable_id'), df0 = df0)
#
#
# par(mar = c(9, 3, 3, 3))
# test = SimSummary(name = 'Procedure_group', df0 = df0)
# barplot(
#   t(as.matrix(test[, c('NormalFDR', 'WindowedFDR')])),
#   horiz = 0 ,
#   beside = 1,
#   col = 1:2,
#   names.arg = unique(test$Procedure_group),
#   las = 3,
#   ylab = 'Total FP'
# )
# legend(
#   'top',
#   legend = c('Normal analysis', 'Soft windowed analysis'),
#   fill = 1:2,
#   horiz = TRUE
# )
#
#
# SimBarplot = function(name = '_CBC_', df0) {
#   par(mar = c(13, 6, 5, 3))
#   cbc = SimSummary(name = c('Parameter_stable_id'), df0 = df0)
#   cbc = cbc[grepl(pattern = name, x = cbc$Parameter_stable_id), ]
#   barplot(
#     t(as.matrix(cbc[, c('NormalFDR', 'WindowedFDR')])),
#     horiz = 0 ,
#     beside = 1,
#     col = 1:2,
#     names.arg = unique(cbc$Parameter_stable_id),
#     las = 3,
#     ylab = 'Total FP',
#     cex.axis = 1.5,
#     cex.names = 1.2,
#     cex.lab = 2
#   )
#   legend(
#     'top',
#     legend = c('Normal analysis', 'Soft windowed analysis'),
#     fill = 1:2,
#     horiz = FALSE,
#     cex = 1.5
#   )
# }
#
# SimBarplot(df0 = df0)
# par(mfrow = c(1, 2))
# SimBarplot(df0 = df0, name = '_DXA_')
# SimBarplot(df0 = df0, name = '_CBC_')
#
# SimBarplot(df0 = df0, name = '_EYE_')
# SimBarplot(df0 = df0, name = '_OFD_')
# SimBarplot(df0 = df0, name = '_HEM_')
# SimBarplot(df0 = df0, name = '_XRY_')
