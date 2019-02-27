vectorOutputAgeing = function(object,
															JSON = FALSE,
															debug = FALSE,
															...) {
	if (is.null(object)) {
		message0('Null object')
		return(NULL)
	}
	out  = tryCatch(
		expr = {
			out = NULL
			suppressMessagesANDErrors(if (is(object, 'PhenStatAgeingMM')) {
				out = vectorOutputCont(object = object,
															 debug = debug)
			} else if (is(object, 'PhenStatAgeingFE')) {
				out = vectorOutputCat(object = object)
			} else if (is(object, 'PhenStatAgeingRR')) {
				out = vectorOutputRR(object = object)
			} else{
				message0('The input object is not of a proper class of PhenStatAgeing')
				out = vectorOutputNULL(object = NULL)
			}, debug = debug)
			
			# JSON engine
			if (JSON && !is.null(out)) {
				for (i in 1:10) {
					out = toJSON(
						out,
						auto_unbox = TRUE,
						null = 'null',
						digits = NA,
						...
					)
					if (i != 10)
						out = fromJSON(txt = out)
				}
			}
			return(out)
		},
		warning = function(war) {
			message0('The functions failed with a warning (see below): ')
			warning(war)
			return(NULL)
		},
		error = function(err) {
			message0('The functions failed with an error (see below): ')
			warning(err)
			return(NULL)
		}
	)
	return(invisible(out))
}
