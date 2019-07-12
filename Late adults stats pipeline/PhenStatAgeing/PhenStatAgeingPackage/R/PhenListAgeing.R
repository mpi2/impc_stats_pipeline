setClass(
	'PhenListAgeing',
	representation(
		datasetPL                   = 'data.frame'    ,
		refGenotype                 = 'character'     ,
		testGenotype                = 'character'     ,
		hemiGenotype                = 'character'     ,
		clean.dataset               = 'logical'       ,
		dataset.colname.genotype    = 'character'     ,
		dataset.colname.sex         = 'character'     ,
		dataset.colname.batch       = 'character'     ,
		dataset.colname.lifestage   = 'character'     ,
		dataset.colname.weight      = 'character'     ,
		dataset.values.missingValue = 'character'  ,
		dataset.values.male         = 'character'     ,
		dataset.values.female       = 'character'     ,
		dataset.values.early        = 'character'     ,
		dataset.values.late         = 'character'     ,
		datasetUNF                  = 'data.frame'
	)
)

PhenListAgeing =
	function(dataset,
					 testGenotype                = 'experimental',
					 refGenotype                 = 'control'     ,
					 hemiGenotype                = NULL          ,
					 clean.dataset               = TRUE          ,
					 dataset.colname.genotype    = 'biological_sample_group',
					 dataset.colname.sex         = 'sex'                    ,
					 dataset.colname.batch       = 'date_of_experiment'     ,
					 dataset.colname.lifestage   = 'LifeStage'              ,
					 dataset.colname.weight      = 'weight'                 ,
					 dataset.values.missingValue = " ",
					 dataset.values.male         = NULL,
					 dataset.values.female       = NULL,
					 dataset.values.early        = NULL,
					 dataset.values.late         = NULL)
	{
		testGenotype = as.character(testGenotype)
		refGenotype  = as.character(refGenotype)
		if (is.null(dataset) || class(dataset) != "data.frame") {
			message0('error ~> Null dataset or not a data.frame.')
			return(NULL)
		}
		dataset_unfiltered = dataset
		if (clean.dataset) {
			sta.time    = Sys.time()
			message0('Checking the input data in progress ...')
			## Replace missing values specified in the user format with NA
			if (!is.null(dataset.values.missingValue)) {
				message0(
					'Checking for the specified missing values (`',
					dataset.values.missingValue,
					'`) ...'
				)
				dataset[dataset %in%  dataset.values.missingValue] = NA
			}
			dataset  [dataset == ""]     = NA
			dataset                      = droplevels0(dataset)
			chkcols = checkPhenlistColumns(
				dataset = dataset,
				vars    = c(
					dataset.colname.genotype   ,
					dataset.colname.sex        ,
					dataset.colname.batch      ,
					dataset.colname.lifestage  ,
					dataset.colname.weight
				)
			)
			dataset = dataset[, order(names(dataset)), drop = FALSE]
			###############################
			dataset = fIsBecome(dataset,
													is       =  dataset.colname.genotype,
													rename   = 'Genotype',
													ifexist  = 'Genotype')
			if (!'Genotype' %in% names(dataset) ||
					all(!chkcols)              ||
					length(c(testGenotype, refGenotype)) != 2) {
				message0(
					'error ~> Please make sure `dataset.colname.xxx` and/or Genotype levels are properly specified'
				)
				return(NULL)
			}
			###############################
			dataset = fIsBecome(dataset,
													is = dataset.colname.sex,
													rename =  'Sex',
													ifexist = 'Sex')
			dataset = fIsBecome(
				dataset = dataset,
				is      = dataset.colname.lifestage,
				rename  = 'LifeStage',
				ifexist = 'LifeStage'
			)
			dataset = fIsBecome(
				dataset = dataset,
				is      = c(dataset.colname.batch),
				rename  = 'Batch',
				ifexist = 'Batch'
			)
			dataset = fIsBecome(dataset,
													is = dataset.colname.weight,
													rename =  'Weight',
													ifexist = 'Weight')
			dataset = factorise  (dataset, c('Genotype', 'Sex', 'Batch', 'LifeStage'))
			dataset = droplevels0(dataset)
			## Renew levels
			dataset = PhenListAgeingRelabling(
				dataset = dataset,
				col = 'Sex',
				l1 = dataset.values.female,
				rel1 = 'Female',
				l2 = dataset.values.male,
				rel2 = 'Male'
			)
			dataset = PhenListAgeingRelabling(
				dataset = dataset,
				col = 'LifeStage',
				l1 = dataset.values.early,
				rel1 = 'Early',
				l2 = dataset.values.late,
				rel2 = 'Late'
			)
			## Hemi to test genotype replacement
			if (!is.null(hemiGenotype)) {
				if (any(rownames(dataset[dataset$Genotype == hemiGenotype,]))) {
					levels(dataset$Genotype)[levels(dataset$Genotype) == hemiGenotype] =
						testGenotype
					message0(
						"Hemizygotes `",
						hemiGenotype,
						"` have been relabelled to test genotype `",
						testGenotype,
						'`'
					)
				}
			}
			## Clean genotypes
			if (length(setdiff(rownames(dataset),
												 rownames(dataset[dataset$Genotype %in% c(testGenotype, refGenotype),]))) >
					0) {
				if (any(!c(testGenotype, refGenotype) %in% levels(dataset$Genotype))) {
					message0(
						'error ~> Mismatch between `Genotype` levels and input levels.',
						'\n\t Genotype Levels: ',
						pasteComma(sort0(levels(
							dataset$Genotype
						))),
						'\n\t Input levels   : ',
						pasteComma(sort0(c(
							testGenotype, refGenotype
						)), replaceNull = TRUE)
					)
					return(NULL)
				}
				dataset = subset(dataset,
												 dataset$Genotype %in% c(testGenotype, refGenotype))
				message0(
					"Dataset has been cleaned to only keep the `Genotype` values: ",
					pasteComma(testGenotype,
										 refGenotype, replaceNull = TRUE)
				)
			}
			## Clean the empty records!
			dataset = CleanEmptyRecords(dataset , c('Genotype', 'Sex', 'Batch', 'LifeStage'))
			## CHECKS
			dataset = checkDataset(droplevels(dataset),
														 testGenotype,
														 refGenotype)
			if (is.null(dataset))
				return(NULL)
			if ('Weight' %in% colnames(dataset)) {
				if (!is.numeric(dataset$Weight)) {
					message0("`Weight` values are not numeric then renamed to `Weight_labels`")
					colnames(dataset)[colnames(dataset) == 'Weight'] =
						'Weight_labels'
				}
				wsglvls = tapply(X = dataset$Weight, INDEX = interaction(dataset$Genotype, dataset$Sex), function(x) {
					length(na.omit(x))
				}, default = 0)
				message0(
					'Total `Weight` data points for Genotype/Sex interaction.\n\t Level(frequency): ',
					pasteComma(paste0(names(wsglvls), '(', wsglvls, ')'), truncate = FALSE)
				)
				if (min(wsglvls, na.rm = TRUE) <= 2) {
					message0(
						"`Weight` column has (<2) data points for at least one Genotype/Sex interaction, then renamed to `Weight_labels`"
					)
					colnames(dataset)[colnames(dataset) == 'Weight'] =
						'Weight_labels'
					
				}
			}
			message0('Successfully performed checks in ',
							 round(difftime(Sys.time() , sta.time, units = 'sec'), 2),
							 ' seconds.')
		} else{
			message0('No check performed on the input data')
		}
		r = new(
			"PhenListAgeing",
			datasetPL = as.data.frame(dataset),
			refGenotype = as.character(refGenotype),
			testGenotype = as.character(testGenotype),
			hemiGenotype = ifelse(is.null(hemiGenotype), character(0), hemiGenotype),
			dataset.colname.batch = ifelse(
				is.null(dataset.colname.batch),
				character(0),
				dataset.colname.batch
			),
			dataset.colname.lifestage = ifelse(
				is.null(dataset.colname.lifestage),
				character(0),
				dataset.colname.lifestage
			),
			dataset.colname.genotype = ifelse(
				is.null(dataset.colname.genotype),
				character(0),
				dataset.colname.genotype
			),
			dataset.colname.sex = ifelse(
				is.null(dataset.colname.sex),
				character(0),
				dataset.colname.sex
			),
			dataset.colname.weight = ifelse(
				is.null(dataset.colname.weight),
				character(0),
				dataset.colname.weight
			),
			dataset.values.missingValue = ifelse(
				is.null(dataset.values.missingValue),
				character(0),
				dataset.values.missingValue
			),
			dataset.values.male = ifelse(
				is.null(dataset.values.male),
				character(0),
				as.character(dataset.values.male)
			),
			dataset.values.female = ifelse(
				is.null(dataset.values.female),
				character(0),
				as.character(dataset.values.female)
			),
			dataset.values.early = ifelse(
				is.null(dataset.values.early),
				character(0),
				as.character(dataset.values.early)
			),
			dataset.values.late = ifelse(
				is.null(dataset.values.late),
				character(0),
				as.character(dataset.values.late)
			),
			clean.dataset = as.logical(clean.dataset),
			datasetUNF = as.data.frame(dataset_unfiltered)
		)
		return(invisible(r))
	}
#-------------------------------------------------------------------------------
## Check dataset for the minimum required info and additional cleaning steps
checkDataset = function(dataset,
												testGenotype,
												refGenotype = "+/+")
{
	dataset = droplevels(dataset)
	if (all(c('Genotype', 'Sex') %in% colnames(dataset))) {
		InGS = interaction(dataset$Genotype, dataset$Sex)
		tbGS = table(InGS)
		message0(
			'Total samples in Genotype:Sex interaction.\n\t Level(frequency): ',
			pasteComma(paste0(names(tbGS), '(', tbGS, ')'), truncate = FALSE)
		)
		if (min(tbGS) < 1)
			message0(
				'No observations detected in the Genotype:Sex interaction for:\n\t ',
				pasteComma(names(tbGS[tbGS < 1]), truncate = FALSE)
			)
		dataset = droplevels(dataset[InGS %in% names(tbGS[tbGS >= 1]), , drop = FALSE])
		## Check of genotype and sex levels after cleaning
		if (nlevels(dataset$Genotype) != 2) {
			message0(
				"error ~> `Genotype` column must have two levels. Current levels: ",
				pasteComma(levels(dataset$Genotype))
			)
			return(NULL)
		}
		if (nlevels(dataset$Sex) > 2) {
			message0(
				"error ~> `Sex` column must have one or two levels. Current levels: ",
				pasteComma(levels(dataset$Sex))
			)
			return(NULL)
		}
		## Check for sex levels - we want to have 'Female' and/or 'Male' only
		wrong_sex_levels = setdiff(levels(dataset$Sex), c("Female", "Male"))
		if (!length(wrong_sex_levels) == 0) {
			message0(
				'error ~> Sex has undefined levels. See: ',
				pasteComma(wrong_sex_levels, truncate = FALSE)
			)
			return(NULL)
		}
		## Check for reference genotype records
		if (refGenotype %in% levels(dataset$Genotype))
			dataset$Genotype = relevel(dataset$Genotype, ref = refGenotype)
	}
	return(dataset)
}

PhenListAgeingBuilder = function(PhenListobject,
																 DOE = NULL,
																 DOB = NULL,
																 d.threshold = 16 * 7,
																 debug       = TRUE) {
	#### Negative age will be removed
	PhenListobject@datasetPL           = droplevels(PhenListobject@datasetPL)
	
	Ageing = !is.null(DOE) &&
		!is.null(DOB) &&
		all(c(DOE, DOB) %in% names(PhenListobject@datasetPL))
	
	if (Ageing) {
		age.in.day = as.Date(PhenListobject@datasetPL[, DOE]) - as.Date(PhenListobject@datasetPL[, DOB])
		PhenListobject@datasetPL$Age       = as.numeric(age.in.day)
		PhenListobject@datasetPL$LifeStage = ifelse    (age.in.day > d.threshold, 'Late', 'Early')
		PhenListobject@datasetPL$LifeStage = as.factor (PhenListobject@datasetPL$LifeStage)
		PhenListobject@datasetPL           = PhenListobject@datasetPL[PhenListobject@datasetPL$Age > 0,]
		message0 ('Age range: ', paste0(range(age.in.day), collapse = '-'))
	} else{
		message0('DOE and DOB are not specified. Then PhenList is returned.')
	}
	
	LL = levels(PhenListobject@datasetPL$LifeStage)
	LS = levels(PhenListobject@datasetPL$Sex)
	LG = levels(PhenListobject@datasetPL$Genotype)
	
	if (length(LL) != 2 && debug && Ageing) {
		message0(
			'Ageing pipeline requires two levels in the LifeStage. Levels: ',
			pasteComma(LL, replaceNull = FALSE, truncate = FALSE),
			'\nNormal pipeline will apply to this data'
		)
		#return(NULL)
	}
	if (length(LS) != 2 && debug) {
		message0(
			'There should be two levels in Sex. Levels: ',
			pasteComma(LS, replaceNull = FALSE, truncate = FALSE)
		)
		#return(NULL)
	}
	if ('Weight' %in% names(PhenListobject@datasetPL) &&
			sum(is.na(PhenListobject@datasetPL[, 'Weight'])) > 0 && debug)
		message0('There are ', sum(is.na(PhenListobject@datasetPL[, 'Weight'])), ' NAs in body weights.')
	
	if (length(LG) < 2 && debug) {
		message0(
			'Genotype must have two levels. Levels: ',
			pasteComma(LG, replaceNull = FALSE, truncate = FALSE)
		)
		return(NULL)
	}
	PhenListobject        = unclass(PhenListobject)
	class(PhenListobject) = 'PhenListAgeing'
	return(PhenListobject)
}


summary.PhenListAgeing = function(object,
																	vars = NULL,
																	...) {
	requireNamespace("summarytools")
	r  = NULL
	df = SelectVariablesOrDefault(data = object@datasetPL, vars)
	if (!is.null(df))
		r  = summarytools::dfSummary (df, ...)
	return(r)
}


plot.PhenListAgeing = function(x,
															 vars = NULL,
															 ...) {
	requireNamespace("Hmisc")
	#options(grType = 'plotly')
	df = SelectVariablesOrDefault(data = x@datasetPL, vars)
	if (!is.null(df)) {
		r    = Hmisc::describe (df)
		plot = suppressWarnings(plot(r, ...))
		return(plot)
	} else{
		message0('No variable found in the data. Please make sure that `vars` exist in the input data')
		return(NULL)
	}
}


SelectVariablesOrDefault = function(data, vars = NULL) {
	cnames = c('Genotype'     ,
						 'Sex'          ,
						 'LifeStage'    ,
						 'Batch'        ,
						 'age_in_weeks' ,
						 'phenotyping_center',
						 'metadata_group'
						 )
	###################
	cnamesF = varExistsInDF(data = data, if (is.null(vars)) {
		cnames
	} else{
		vars
	})
	if (is.null(cnamesF) || length(cnamesF) < 1) {
		message0(
			'Non of the specified variables exist in the data. See the list below:\n\t',
			pasteComma(cnames, truncate = FALSE)
		)
		return(NULL)
	}
	return(data[, cnamesF, drop = FALSE])
}

varExistsInDF = function(data = NULL, vars = NULL) {
	if (is.null(vars) || is.null(data))
		return(NULL)
	r = vars[vars %in% names(data)]
	return(r)
}
