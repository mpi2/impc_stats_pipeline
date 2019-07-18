vectorOutputAgeing = function(object,
															othercolumns       = NULL,
															JSON               = FALSE,
															RemoveNullKeys     = FALSE,
															ReportNullSchema   = FALSE,
															...) {
	requireNamespace("rlist")
	debug = FALSE
	if ((is.null(object) ||
			 !is.null(object$messages)) && !ReportNullSchema) {
		message0('Null object')
		if (!is.null(object$messages)) {
			message0('Please see the error below:')
			print(object$messages)
		}
		return(NULL)
	}
	##########
	out  = tryCatch(
		expr = {
			out = NULL
			suppressMessagesANDWarnings(
				if (is(object,
							 'PhenStatAgeingMM')) {
					out = vectorOutputCont(object = object,
																 debug = debug)
				} else if (is(object, 'PhenStatAgeingFE')) {
					out = vectorOutputCat(object = object)
				} else if (is(object, 'PhenStatAgeingRR')) {
					out = vectorOutputRR(object = object)
				} else{
					if (RemoveNullKeys) {
						message0('`RemoveNullKeys` does not apply to the NULL schema.')
						RemoveNullKeys = FALSE
					}
					out = vectorOutputNULL(object = NULL)
					out$`Additional information`$messages = UnlistCall(object$messages)
				},
				sup.messages = !debug,
				sup.warnings = FALSE
			)
			#########
			if (!is.null(object$input$PhenListAgeing)) {
				NewNames = variablesInData(
					df = object$input$PhenListAgeing@datasetPL,
					names = othercolumns,
					debug = debug
				)
				if (!is.null(out)    &&
						!is.null(NewNames)) {
					out$othercolumns = as.list(object$input$PhenListAgeing@datasetPL[, NewNames, drop = FALSE])
				} else{
					out$othercolumns = NULL
				}
			}
			if (RemoveNullKeys && !is.null(out)) {
				for (i in 1:5) {
					out = list.clean(
						out,
						fun = function(x) {
							length(x) == 0L || is.null(x)
						},
						recursive = TRUE
					)
				}
			}
			#########
			# JSON engine
			requireNamespace("jsonlite")
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
			message(err)
			return(NULL)
		}
	)
	return(invisible(out))
}
