# #### Summary core for continues data
# summary.PhenListAgeingModel = function(object, ...) {
# 	if (is.null(object))
# 		stop('\n ~> NULL object\n')
# 	tmpobject  = object
#
# 	ai = formulaTerms(object$initial.model)
# 	af = formulaTerms(object$final.model)
# 	if (is.null(object$final.model$call$random)) {
# 		ar = NULL
# 	} else{
# 		ar = object$final.model$call$random
# 	}
# 	#notSigVars = ai[!(ai %in% af)]	#SigVars    = ai[(ai %in% af)]
# 	tmpobject$final.model$call$data = tmpobject$final.model$call$family = NULL
# 	F.Sum = summary(tmpobject$final.model) #main
# 	message0(
# 		'~> Final model: ',
# 		'\n\t Fixed: ' ,
# 		ifelse(is.null(af), '-', paste(af, collapse = ', ')),
# 		'\n\t Random: ' ,
# 		ifelse(is.null(ar), '-', format(ar)),
# 		'\n\t ----- ' ,
# 		'\n\t Removed term(s): ',
# 		ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
# 		'\n\t ----- ' ,
# 		'\n\t Genotype effect size: ',
# 		ifelse(
# 			object$LifeStage,
# 			paste(
# 				c(
# 					'\n\t  Genotype for Early ~>',
# 					'\n\t  Genotype for Late ~>' ,
# 					'\n\t  LifeStage ~>'
# 				),
# 				object$effect,
# 				collapse = ', '
# 			),
# 			object$effect
# 		),
# 		'\n\t ----- ' ,
# 		'\n\t Variance homogeneity: ',
# 		ifelse(!is.null(object$VarHomo), object$VarHomo, 'Not specified'),
# 		'\n\t ----- ' ,
# 		'\n\n~> Model summary:\n'
# 	)
# 	print(F.Sum)
# 	outp = list(
# 		Initial.Terms = ai   ,
# 		Final.Terms = af     ,
# 		RandoEffect = ar     ,
# 		removed.terms = ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
# 		object  = object     ,
# 		Summary = summary(object$final.model)
# 	)
# 	outp$JSON = toJSONI(outp$Summary)
# 	return(invisible(outp))
# }
#
