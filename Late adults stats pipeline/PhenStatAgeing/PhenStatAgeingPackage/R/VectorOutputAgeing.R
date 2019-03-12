vectorOutputAgeing = function(object,
															othercolumns = NULL,
															JSON = FALSE,
															debug = FALSE,
															...) {
	if (!is.null(object$messages)) {
		message0('Null object. Please see the error below:')
		print(object$messages)
		return(NULL)
	}
	##########
	out  = tryCatch(
		expr = {
			out = NULL
			suppressMessagesANDErrors(if (is(object,
																			 'PhenStatAgeingMM')) {
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
			
			#########
			NewNames = variablesInData(df = object$input$PhenListAgeing@datasetPL,
																 names = othercolumns,
																 debug = debug)
			if (!is.null(out)    &&
					!is.null(NewNames)) {
				out$othercolumns = as.list(object$input$PhenListAgeing@datasetPL[, NewNames, drop = FALSE])
			} else{
				out$othercolumns = NULL
			}
			#########
			# JSON engine
			n   = 5
			if (JSON && !is.null(out)) {
				for (i in 1:n) {
					out = toJSON(
						out,
						auto_unbox = TRUE,
						null = 'null',
						digits = NA,
						...
					)
					if (i != n)
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
