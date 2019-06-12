NewPhenList =
	function(dataset,
					 testGenotype,
					 refGenotype = '+/+',
					 hemiGenotype = NULL,
					 clean.dataset = TRUE,
					 dataset.colname.batch = NULL,
					 dataset.colname.genotype = NULL,
					 dataset.colname.sex = NULL,
					 dataset.colname.weight = NULL,
					 dataset.values.missingValue = " ",
					 dataset.values.male = NULL,
					 dataset.values.female = NULL)
	{
		if (is.null(dataset) || class(dataset) != "data.frame") {
			message0('error ~> Null dataset or not a data.frame.')
			return(NULL)
		}
		dataset_unfiltered = dataset
		if (clean.dataset) {	
			chkcols = checkPhenlistColumns(
				dataset = dataset,
				vars    = c(
					dataset.colname.genotype,
					dataset.colname.sex     ,
					dataset.colname.batch   ,
					dataset.colname.weight
				)
			)
			if (any(head(!chkcols, 2))) {
				message0('error ~> Please make sure that the `dataset.colname.xxx` is set properly')
				return(NULL)
			}
			dataset = factorise(dataset, c('Genotype', 'Sex', 'Batch'))
			dataset = dataset[, order(names(dataset)), drop = FALSE]
			dataset = fIsBecome(
				dataset = dataset,
				is      = c(dataset.colname.batch, 'Assay.Date', 'AssayDate'),
				rename  = 'Batch',
				ifexist = 'Batch'
			)
			dataset = fIsBecome(dataset,
													is = dataset.colname.genotype,
													rename =  'Genotype',
													ifexist = 'Genotype')
			dataset = fIsBecome(
				dataset,
				is = c(dataset.colname.sex, 'Gender'),
				rename =  'Sex',
				ifexist = 'Sex'
			)
			dataset = fIsBecome(dataset,
													is = dataset.colname.weight,
													rename =  'Weight',
													ifexist = 'Weight')
			
			## Replace missing values specified in the user format with NA
			if (!is.null(dataset.values.missingValue)) {
				message0(
					'Checking for the specified missing values (',
					dataset.values.missingValue,
					') ...'
				)
				dataset[dataset %in%  dataset.values.missingValue] = NA
			}
			dataset  [dataset == ""]                             = NA
			# make Weight column numeric if possible (if there are no strings)
			dataset = droplevels0(dataset)
			## Renew levels
			if ('Sex' %in% colnames(dataset)) {
				## Replace values for sexes with 'Male','Female'
				levels(dataset$Sex) = capitalise(levels(dataset$Sex))
				if (!is.null(dataset.values.female))
					levels(dataset$Sex)[levels(dataset$Sex) %in% dataset.values.female] =
						"Female"
				if (!is.null(dataset.values.male))
					levels(dataset$Sex)[levels(dataset$Sex) %in% dataset.values.male] =
						"Male"
			}
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
						'\n Genotype Levels: ',
						pasteComma(sort(levels(
							dataset$Genotype
						))),
						'\n Input levels   : ',
						pasteComma(sort(c(
							testGenotype, refGenotype
						)))
					)
					return(NULL)
				}
				dataset = subset(dataset,
												 dataset$Genotype %in% c(testGenotype, refGenotype))
				message0(
					"Dataset has been cleaned to only keep the genotype values, ",
					testGenotype,
					" and ",
					refGenotype
				)
			}
			## Clean the empty records  NB - after renaming/cleaning !
			dataset = CleanEmptyRecords(dataset , c('Genotype', 'Sex', 'Batch'))
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
					'Total `Weight` data points for Genotype/Sex interaction.\n\t  Level(frequency): ',
					pasteComma(paste0(names(wsglvls), '(', wsglvls, ')'), truncate = FALSE)
				)
				if (max(wsglvls, na.rm = TRUE) <= 2) {
					message0(
						"`Weight` column has (<2) data points for Genotype/Sex interaction, then renamed to `Weight_labels`",
					)
					colnames(dataset)[colnames(dataset) == 'Weight'] =
						'Weight_labels'
					
				}
			}
		}else{
			message0('No check performed on the input data')
		}
		new(
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
			clean.dataset = as.logical(clean.dataset),
			datasetUNF = as.data.frame(dataset_unfiltered)
		)
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
			'Total samples in Genotype:Sex interaction.\n\t level(frequency): ',
			pasteComma(paste0(names(tbGS), '(', tbGS, ')'), truncate = FALSE)
		)
		if (min(tbGS) <= 2)
			message0(
				'Less than 2 observations detected in the Genotype:Sex interaction for:\n\t ',
				pasteComma(names(tbGS[tbGS <= 2]), truncate = FALSE)
			)
		dataset = droplevels(dataset[InGS %in% names(tbGS[tbGS > 2]), , drop = FALSE])
		## Check of genotype and sex levels after cleaning
		if (nlevels(dataset$Genotype) != 2) {
			message0("error ~> `Genotype` column must have two levels. Current levels: ",
							 pasteComma(levels(dataset$Genotype)))
			return(NULL)
		}
		if (nlevels(dataset$Sex) > 2) {
			message0("error ~> `Sex` column must have one or two levels. Current levels: ",
							 pasteComma(levels(dataset$Sex)))
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
