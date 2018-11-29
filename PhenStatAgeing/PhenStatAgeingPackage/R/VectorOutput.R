vectorOutput <- function(object, phenotypeThreshold = 0.01)
{
	if (is.null(object))
		return (NULL)
	
	o1          = summary(object)
	depVariable = o1$object$depVariable
	# MM
	if (o1$object$method %in% c("MM", 'GLM')) {
		equation       = paste('Formula: ',format(formula(o1$object$final.model)),sep = '')
		framework      = switch(object$method,
														MM  = "Linear Mixed Model framework",
														GLM = "Generalized Linear Model framework")
		fittingMethod    = "Maximum Likelihood using the Generalised least squares, "
		x                = nlme::getData(o1$object$final.model)
		columnOfInterest = x[, c(depVariable)]
		variability      = paste('variability:',
														 length(unique(columnOfInterest)) / length(columnOfInterest),
														 sep = "")
		DSsize           = paste(ConvDf2Flat(xtabs(
			formula        = reformulate(
				termlabels   = c('Genotype', 'Sex', 'LifeStage'),
				response     = NULL
			),
			data           = x
		)), collapse     = '; ')
		addInfo           = paste("{", DSsize, variability,equation, "}", sep = 	";")
		percentageChanges = NA
		vectorOutput      = list(
			'Method'                               = 	paste(framework, ", ", fittingMethod, format(equation), sep = ""),		
			'Dependent variable'                   =	depVariable,		
			'Batch included'                       =	any(grepl('Batch', x = c(o1$RandoEffect, o1$Final.Terms))),		
			'Residual variances homogeneity'       =	o1$object$VarHomo,		
			'Genotype contribution'                =	ifelse ('Genotype' %in% o1$Final.Terms, TRUE, FALSE),		
			'Genotype estimate'                    =	ifelse ('Genotype' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Genotype') + 1, 1], NA),		
			'Genotype standard error'              =	ifelse ('Genotype' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Genotype') + 1, 2], NA),		
			'Genotype p-Val'                       =	ifelse ('Genotype' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Genotype') + 1, 5], NA),		
			'Genotype percentage change'           =	NA,#as.character(percentageChanges),		
			'Sex estimate'                         =	ifelse ('Sex'    %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms ==	'Sex')+1,1],NA)	,
			'Sex standard error'                   =	ifelse ('Sex'    %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms ==	'Sex')+1,2],NA)	,
			'Sex p-val'                            =	ifelse ('Sex'    %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms == 'Sex') + 1, 5], NA),		
			'Weight estimate'                      =	ifelse ('Weight' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Weight') + 1, 1], NA)	,	
			'Weight standard error'                =	ifelse ('Weight' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Weight') + 1, 2], NA)	,	
			'Weight p-val'                         =	ifelse ('Weight' %in% o1$Final.Terms, o1$Summary$tTable[which(o1$Final.Terms =='Weight') + 1, 5], NA)	,	
			'Gp1 genotype'                         =	PhenStat:::testGenotype(o1$object$data)		,
			'Gp1 Residuals normality test'         =	NA, #RnnnTest	,	
			'Gp2 genotype'                         =	PhenStat:::refGenotype(o1$object$data)		,
			'Gp2 Residuals normality test'         =	NA, #RnotmTest		,
			'Intercept estimate'                   =	as.character(o1$Summary$tTable[1, 1]),		
			'Intercept standard error'             =	as.character(o1$Summary$tTable[1, 2]),		
			'Intercept p-val'                      =	as.character(o1$Summary$tTable[1, 5]),		
			'Interaction(s) included'              =	paste(ifelse(any(grepl(pattern = ':', o1$Final.Terms)), o1$Final.Terms[grep(pattern = ':', o1$Final.Terms)], FALSE), collapse = ','),		
			'Interaction(s) p-val'                 =	paste(as.character(ifelse(any(grepl(pattern = ':', rownames(o1$Summary$tTable))), o1$Summary$tTable[grepl(pattern = ':', rownames(o1$Summary$tTable)), 5], FALSE)), collapse = ','),		
			'Sex FvKO estimate'                    =	NA,# as.character(analysisResults$model.output.summary["sex_FvKO_estimate"]),		
			'Sex FvKO standard error'              =	NA,# as.character(analysisResults$model.output.summary["sex_FvKO_SE"]),		
			'Sex FvKO p-val'                       =	NA,# as.character(analysisResults$model.output.summary["sex_FvKO_p_value"]),		
			'Sex MvKO estimate'                    =	NA,# as.character(analysisResults$model.output.summary["sex_MvKO_estimate"]),		
			'Sex MvKO standard error'              =	NA,# as.character(analysisResults$model.output.summary["sex_MvKO_SE"]),		
			'Sex MvKO p-val'                       =	NA,# as.character(analysisResults$model.output.summary["sex_MvKO_p_value"]),		
			'Classification tag'                   =	NA,# as.character(classificationValue),		
			'Transformation'                       =	FALSE,#transformation(object),		
			'Additional information'               =	addInfo				
		)
	}
	else if (object$method %in% c("FE")) {
		vectorOutput =list(
			'Method'                                     =	'Fisher Exact Test framework',
			'Dependent variable'                         =	depVariable,
			'Batch included '                            =	NA,
			'Residual variances homogeneity'             =	NA,
			'Genotype contribution'                      =	NA,
			'Genotype estimate'                          =	o1$summary$GenEsize,
			'Genotype standard error'                    =	NA,
			'Genotype p-Val'                             =	o1$summary$GenPval,
			'Genotype percentage change'                 =	NA,
			'Sex estimate'                               =	NA,
			'Sex standard error'                         =	NA,
			'Sex p-val'                                  =	NA,
			'Weight estimate'                            =	NA,
			'Weight standard error'                      =	NA,
			'Weight p-val'                               =	NA,
			'Gp1 genotype'                               =	PhenStat:::testGenotype(o1$object$data),
			'Gp1 Residuals normality test'               =	NA,
			'Gp2 genotype'                               =	PhenStat:::refGenotype(o1$object$data),
			'Gp2 Residuals normality test'               =	NA,
			'Blups test'                                 =	NA,
			'Rotated residuals normality test'           =	NA,
			'Intercept estimate'                         =	NA,
			'Intercept standard error'                   =	NA,
			'Interaction included'                       =	NA,
			'Interaction p-val'                          =	NA,
			'Sex FvKO estimate'                          =	o1$summary$SexEsize,
			'Sex FvKO standard error'                    =	NA,
			'Sex FvKO p-val'                             =	NA,
			'Sex MvKO estimate'                          =	NA,
			'Sex MvKO standard error'                    =	NA,
			'Sex MvKO p-val'                             =	NA,
			'Classification tag'                         =	NA,
			'Transformation'                             =	NA,
			'Additional information'                     =	NA
		)
	}
	
	vectorOutput = as.matrix(vectorOutput)
	return(vectorOutput)
}

#-------------------------------------------------------------------------------
vectorOutputMatrices <-
	function(object, vars = c('Sex','Genotype','LifeStage')) {
		
		if(is.null(object))
			return (NULL)
		
		TestSet      = paste('Gp1 Genotype (g1)',PhenStat:::testGenotype(object$data),sep = ': ')
		ReferenceSet = paste('Gp2 Genotype (g2)',PhenStat:::refGenotype (object$data),sep = ': ')
		DVariable    = paste('Dependent variable',object$depVariable,sep = ': ')
		
		opt = sapply(object$data@datasetPL[,vars],levels,simplify = FALSE)
		
		out = NULL
		for(i in seq(length(opt))){
			out=c(out,paste(i,'.',names(opt)[[i]],'|',opt[[i]],sep = ''))
		}
		
		vectorOutput = c(DVariable,TestSet,ReferenceSet,out)
		return (vectorOutput)
		
		
	}
##------------------------------------------------------------------------------
