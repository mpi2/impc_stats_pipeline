# windowing parameters
WindowingDetails_del_me = function(args) {
  wn =
    paste(
      '"window_parameters":{',
      SingleOrMore(
        name = 'l',
        value = args$r$final.l,
        comma = TRUE
      ),
      SingleOrMore(
        name = 'k',
        value = args$r$final.k,
        comma = TRUE
      ),
      SingleOrMore(
        name = 'doe',
        value = unique(args$tt[args$mm]),
        comma = TRUE
      ),
      SingleOrMore(
        name = 'min_obs',
        value = unique(args$r$min.obs),
        comma = TRUE
      ),
      SingleOrMore(
        name = 'obs_in_window',
        value = unique(args$r$finalModel$output$ObsInInterval),
        comma = TRUE
      ),
      SingleOrMore(
        name = 'threshold',
        value = unique(args$threshold),
        comma = TRUE
      ),
      ifelse(
        length(unique(args$tt[args$mm])) > 12,
        '"the_number_of_doe_in_window": "more than 12 date of experiments then a sample of (10) DOE is considered" ,',
        paste0('"the_number_of_doe_in_window": ', length(unique(args$tt[args$mm])), ' ,')
      ),
      SingleOrMore(
        name = 'weights',
        value = args$we,
        comma = FALSE
      ),
      '}',
      sep = ''
    )

  return(wn)
}



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
    weights = as.vector(unlist(args$we))
  )
  return(wn)
}