# windowing parameters
##########
# windowing parameters
WindowingDetails = function(args) {
  # 1 LIST
  wn = list(
    l =  args$r$final.l                   ,
    k = args$r$final.k                    ,
    doe = unique(args$tt[args$mm])        ,
    min_obs = unique(args$r$min.obs)      ,
    obs_in_window = unique(args$r$finalModel$output$ObsInInterval),
    threshold = unique(args$threshold)                            ,
    the_number_of_doe_in_window = length(unique(args$tt[args$mm])),
    doe_note = ifelse(
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
    weights = as.vector(unlist(args$we)),
    external_sample_id = as.vector(unlist(args$external_sample_ids))
  )
  return(wn)
}