# windowing parameters
##########
# windowing parameters
WindowingDetails = function(args) {
  # 1 LIST
  wn = list(
    'l' =  args$r$final.l                   ,
    'k' = args$r$final.k                    ,
    'DOE' = unique(args$tt[args$mm])        ,
    'Min obs required in the window'    = unique(args$r$min.obs)      ,
    'Total obs or weight in the window' = unique(args$r$finalModel$output$ObsInInterval),
    'Threshold'                         = unique(args$threshold)                            ,
    'The number of DOE in the window'   = length(unique(args$tt[args$mm])),
    'DOE note' = ifelse(
      length(unique(args$tt[args$mm])) > args$maxPeaks,
      paste0(
        'More than ',
        args$maxPeaks,
        ' date of experiments but a sample of (',
        args$maxPeaks,
        ') DOE is considered'
      ),
      'no note available'
    ),
    'Window weights' = if (!is.null(args$we)) {
      as.vector(unlist(args$we))
    } else{
      rep(1,length(args$external_sample_ids))
    },
    'external_sample_id' = as.vector(unlist(args$external_sample_ids))
  )
  return(wn)
}