PhenList0 =
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
		if (is.null(dataset) || class(dataset) != "data.frame")
			stop('Null dataset or not a data.frame.')
		dataset_unfiltered = dataset
		dataset = factorise(dataset, c('Genotype', 'Sex', 'Batch'))
		dataset = dataset[, order(names(dataset)), drop = FALSE]
		if (clean.dataset) {
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
			if (!is.null(dataset.values.missingValue))
				dataset[dataset %in%  dataset.values.missingValue] = NA
			dataset  [dataset == ""]                             = NA
			# make Weight column numeric if possible (if there are no strings)
			if ('Weight' %in% colnames(dataset)) {
				columnName = "Weight"
				checkNA_transformed =
					sum(is.na(suppressWarnings(as.numeric(
						as.character(dataset[, c(columnName)])
					))))
				checkNA_initial = sum(is.na(dataset[, c(columnName)]))
				if (checkNA_transformed == checkNA_initial) {
					dataset[, c(columnName)] =
						as.numeric(as.character(dataset[, c(columnName)]))
				}
			}
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
				if (any(rownames(dataset[dataset$Genotype == hemiGenotype, ]))) {
					levels(dataset$Genotype)[levels(dataset$Genotype) == hemiGenotype] =
						testGenotype
					message0(
						"Information: Hemizygotes '",
						hemiGenotype,
						"' have been relabelled to test genotype '",
						testGenotype,
						"'.\nIf you don't want this behaviour then leave ",
						"'hemiGenotype' to NULL,"
					)
				}
			}
			## Clean genotypes
			if (length(setdiff(rownames(dataset),
												 rownames(dataset[dataset$Genotype %in% c(testGenotype, refGenotype), ]))) >
					0) {
				dataset = subset(dataset,
												 dataset$Genotype %in% c(testGenotype, refGenotype))
				message0(
					"Information: Dataset has been cleaned by ",
					"filtering out records with genotype value other than test ",
					"genotype '",
					testGenotype,
					"' or reference genotype '",
					refGenotype
				)
			}
			## Clean the empty records  NB - after renaming/cleaning !
			dataset = CleanEmptyRecords(dataset , c('Genotype', 'Sex', 'Batch'))
			## CHECKS
			dataset = checkDataset(droplevels(dataset),
														 testGenotype,
														 refGenotype)
			checkWeight = columnChecks0(dataset, "Weight", 2)
			
			if (!checkWeight[1]) {
				message0("Information: Weight column is not present in the database")
			}
			else {
				if (!checkWeight[2]) {
					message0(
						"Information: Weight column values are not numeric.\n",
						"In order to avoid erroneous execution of statistical ",
						"functions column is renamed to 'Weight_labels'"
					)
					colnames(dataset)[colnames(dataset) == 'Weight'] =
						'Weight_labels'
				}
				if (!checkWeight[3]) {
					message0(
						"Information: Weight column does not have enough data points ",
						"for genotype/sex combinations.\n",
						"In order to avoid erroneous execution of statistical ",
						"functions column is renamed to 'Weight_labels'"
					)
					colnames(dataset)[colnames(dataset) == 'Weight'] =
						'Weight_labels'
				}
			}
		}
		new(
			"PhenListAgeing",
			datasetPL = as.data.frame(dataset),
			refGenotype = refGenotype,
			testGenotype = testGenotype,
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
			clean.dataset = clean.dataset,
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
	## Column names should be given
	if (ncol(dataset) < 1 || is.null(colnames(dataset))) {
		stop("Check failed: Dataset with no column names.")
	}
	## Check for mandatory columns: Genotype and Sex
	if (!all(c('Genotype', 'Sex') %in% names(dataset))) {
		stop("Check failed: Dataset's 'Genotype' or 'Sex' columns are missed.")
	}
	## Check for other columns: Weight and Batch
	if (!('Weight' %in% colnames(dataset))) {
		message0(
			"Information: Dataset's 'Weight' column is missed.\n",
			"You can define 'dataset.colname.weight' argument to specify column ",
			"for the weight effect modeling"
		)
	}
	if (!('Batch' %in% colnames(dataset))) {
		message(
			"Information: Dataset's 'Batch' column is missed.\n",
			"You can define 'dataset.colname.batch' argument to specify column ",
			"for the batch effect modeling"
		)
	}
	if (all(c('Genotype', 'Sex') %in% colnames(dataset))) {
		InGS = interaction(dataset$Genotype, dataset$Sex)
		tbGS = table(InGS)
		if (min(tbGS) <= 2)
			message0(
				'Data levels with less than 2 observations in the interaction of Genotype:Sex:\n\t',
				pasteComma(names(tbGS[tbGS <= 2]), truncate = FALSE),
				'\nAll available levels (frequency):\n',
				pasteComma(paste0(names(tbGS), '(', tbGS, ')'), truncate = FALSE)
			)
		dataset = droplevels(dataset[InGS %in% names(tbGS[tbGS > 2]), , drop = FALSE])
		## Check of genotype and sex levels after cleaning
		if (nlevels(dataset$Genotype) != 2) {
			stop(
				"Check failed: Dataset's 'Genotype' ",
				"column has to have two values.\nYou can define 'testGenotype' and ",
				"'refGenotype' arguments to automatically filter out records with ",
				"genotype values other than specified.\nAlternatively you can define ",
				"'hemiGenotype' and 'testGenotype' arguments to relabel hemizygotes ",
				"to homozygotes."
			)
		}
		if (nlevels(dataset$Sex) > 2) {
			stop(
				"Check failed: Dataset's 'Sex' ",
				"column has to have one or two values and currently the data has ",
				"more than two."
			)
		}
		## Check for sex levels - we want to have 'Female' and/or 'Male' only
		wrong_sex_levels = setdiff(levels(dataset$Sex), c("Female", "Male"))
		if (!length(wrong_sex_levels) == 0) {
			stop('Sex has undefined levels. See:',
					 pasteComma(wrong_sex_levels, truncate = FALSE))
		}
		## Check for reference genotype records
		#if (sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1)
		if (refGenotype %in% levels(dataset$Genotype))
			dataset$Genotype = relevel(dataset$Genotype, ref = refGenotype)
		else {
			stop(
				"Check failed: Dataset with not ",
				"enough records for statistical analysis with reference genotype '",
				refGenotype
			)
		}
		## Check for test genotype records
		if (!(testGenotype %in% levels(dataset$Genotype))) {
			stop(
				"Check failed: Dataset ",
				"with not enough records for statistical analysis with test ",
				"genotype '",
				testGenotype
			)
		}
	}
	return(dataset)
}
