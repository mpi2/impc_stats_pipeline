vectorOutput <- function(object, phenotypeThreshold = 0.01)
{
	if (is.null(object))
		return (NULL)
	# MM
	if (object$output$Final.Model.Tag %in% c('glm', 'gls','lme')) {
		depVariable = all.vars(formula(object$output$Final.Model))[1]
		equation       = paste('final model: ',format(formula(object$output$Final.Model)),sep = '')
		framework      = switch(object$output$Final.Model.Tag,
														lme  = "Linear Mixed Model framework",
														gls  = "Linear Model Using Generalized Least Squares framework",
														glm  = "Generalized Linear Model framework")
		fittingMethod    = toupper(object$output$Final.Model.Tag)
		x                = object$input$data
		columnOfInterest = x[, c(depVariable)]
		variability      = paste0('variability: ',
														 length(unique(columnOfInterest)) / length(columnOfInterest))
		DSsize            = SummaryStats(x = x,formula = object$input$Fullfixed$initial,label = 'summary_statistics',lower = TRUE,drop = TRUE)
		addInfo           = NULL
		percentageChanges = NA
		o1                = summary(object$output$Final.Model)
		vectorOutput      = list(
			'Method'                               = 	paste(framework, ", ", fittingMethod, ', ', format(equation), sep = ""),		
			'Dependent variable'                   =	depVariable,		
			'Batch included'                       =	object$output$BatchIn,		
			'Residual variances homogeneity'       =	object$output$VarHomoIn,		#modelSummaryPvalueExtract
			'Genotype contribution'                =	ifelse ('Genotype' %in% all.vars0(object$output$Final.Model), TRUE, FALSE),		
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
