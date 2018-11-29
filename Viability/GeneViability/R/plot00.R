plot00 = function(DR7_viability,
                  data ,
                  unique = 'MGI id',
                  xlab = 'Viability',
                  ...) {
  pM = DR7_viability$Viability[DR7_viability[, unique] %in% data$inputdata[data$duplicates,][, unique]]
  plot(
    pM,
    xlab = xlab,
    ylab = 'Frequency of [removed] duplicates with known label',
    #main = 'DMDD',
    col = unique(pM),
    sub = paste0('Total = ', length(pM)),
    ...
  )
}