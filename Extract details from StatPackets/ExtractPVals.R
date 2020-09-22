args = commandArgs(trailingOnly = TRUE)
library(data.table)
library(jsonlite)
####################
unlist0 = function(x, active = TRUE) {
	x2 = unlist(x)
	if ((is.null(x) || is.null(x2)) && active) {
		return(NA)
	} else{
		return(x2)
	}
}

is.nullorNA = function(x) {
	if (is.null(x) || all(x == '')) {
		return(NA)
	} else{
		return(x)
	}
}

SelectOthers = function(object) {
	paste0(is.nullorNA(c(
		unlist(object$`Concurrent control selection`  )
	)))
}

SelectMISC = function(object) {
  ReferenceGene = c('Nxn', 'Rnf10', 'Ap4e1', 'Prkab1', 'Dnase1l2', 'Dbn1')
  rReference = object$V12   %in% ReferenceGene
  rIgnorome  = object$V11  %in% DRrequired:::ignoromeGenes()
  rBehaviour = object$V6   %in% DRrequired:::BehviourParamters()
  return(c(rReference, rIgnorome, rBehaviour))
}

MakeURL = function(object, objectJSON) {
	r = URLencode(
		paste0(
			"https://wwwdev.ebi.ac.uk/mi/impc/dev/phenotype-archive/media/images/windowing/?alleleSymbol=",
			object$V9,
			"&colonyID=",
			object$V19,
			"&procedure= | ",
			object$V3,
			"&parameter= | ",
			object$V6,
			"&center=",
			object$V8,
			"&zygosity=",
			object$V18,
			"&metadata=",
			object$V17
		)
	)
	r = c(r, URLdecode(unlist(
		objectJSON$Result$Details$`Exported raw data file`
	)))
	return(r)
}

SelectWindowingParameters = function(object) {
  r = c(
    ###### Windowing parameters
    unlist0(object$`Window parameters`$l[1])             ,
    unlist0(object$`Window parameters`$k[1])             ,
    unlist0(object$`Window parameters`$`Min obs required in the window`),
    unlist0(object$`Window parameters`$`The number of DOE in the window`) ,
    unlist0(object$`Window parameters`$Threshold)                   ,
    unlist0(length(object$`Window parameters`$`Window weights`))             ,
    unlist0(object$`Window parameters`$`Total obs or weight in the window`)
  )
  return(r)
}

SelectAnalysis = function(object) {
	r = c(
		######
		unlist0(object$`Applied method`)                ,
		unlist0(object$`Classification tag`$`Classification tag`) ,
		#unlist0(object$formula)                                  ,
		unlist0(object$`Residual variances homogeneity`)  ,
		unlist0(object$`Batch included`)                  ,
		# Sexual dymorphism
		unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Criteria)   ,
		unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Note)   ,
		# effect size
		unlist0(object$`Genotype estimate`$Value)               ,
		unlist0(object$`Sex FvKO estimate`$Value)               ,
		unlist0(object$`Sex MvKO estimate`$Value)               ,
		# Genotype
		unlist0(object$`Genotype p-value`)                  ,
		unlist0(object$`Genotype standard error`)           ,
		# sexDim
		unlist0(object$`Sex FvKO p-value`)                  ,
		unlist0(object$`Sex MvKO p-value`)                  ,
		unlist0(object$`Sex FvKO standard error`)           ,
		unlist0(object$`Sex MvKO standard error`)           ,
		# sex
		unlist0(object$`Sex p-value`)                       ,
		unlist0(object$`Sex estimate`$Value)                ,
		unlist0(object$`Sex standard error`)                ,
		# Bodyweight
		unlist0(object$`Weight p-value`)                    ,
		unlist0(object$`Weight estimate`$Value)             ,
		unlist0(object$`Weight standard error`)
	)
	return(r)
}

tableCount = function(Gen,
                      Sex,
                      levels = c('experimental.male',
                                 'experimental.female',
                                 'control.male',
                                 'control.female')) {
  t = table(interaction(Gen, Sex, sep = '.'))
  t = as.data.frame(t)
  ll = length(levels)
  r = c(rep(NA, ll))
  for (i in 1:ll)
    r[i] =   ifelse(levels[i] %in% t[, 1], t[levels[i] == t[, 1], 2], NA)
  return(r)
}



########## Main function
f = function(start, end, file = 'Index_DR101_V1.txt') {
  if (is.na(end))
    end = start
	ofname = paste0('R', '_', start, '-', end, '_pval.tsv')
	if (file.exists(paste0('./resultF/', ofname)))
		unlink(paste0('./resultF/', ofname))
	library(RJSONIO)
	library(DRrequired)
	library(data.table)
	####
	df = readLines(file)
	for (i in start:end) {
		cat(i, '|')
		file = df[i]
		if (!(
		  grepl(pattern = 'Successful',
		        x = file,
		        fixed = TRUE) ||
		  grepl(pattern = 'NotProcessed',
		        x = file,
		        fixed = TRUE) ||
		  grepl(pattern = 'Failed',
		        x = file,
		        fixed = TRUE)
		))
		  next
		if (!file.exists(file)) {
		  if (file.exists(dirname(file),    'output_Successful.tsv'))
		    file = file.path(dirname(file), 'output_Successful.tsv')
		  else if (file.exists(dirname(file), 'output_NotProcessed.tsv'))
		    file = file.path(dirname(file),   'output_NotProcessed.tsv')
		  else
		    write(file,file =  paste0('DoesNotExistFiles_',Sys.Date(),'.txt'),append = TRUE)
		    next
		}

		r0 = fread(
			file = file,
			header = FALSE,
			sep = '\t',
			quote = "",
			stringsAsFactors = FALSE
		)
		rN = DRrequiredAgeing:::annotationChooser(
		  statpacket = r0,
		  level = .0001
		)
		rW = DRrequiredAgeing:::annotationChooser(
		  statpacket = r0,
		  level = .0001,
		  resultKey = 'Windowed result',
		  TermKey = 'WMPTERM'
		)



		if (ncol(r0) != 20){
			write(file,file = paste0('Error_Overal_',Sys.Date(),'.txt'),append = TRUE)
			next
		}
		r1 =
		  tryCatch(
		    expr =
		      fromJSON(r0$V20, flatten = TRUE),
		    warning = function(w) {
		      write(file,file =paste0('Error_inJSON_',Sys.Date(),'.txt'),append = TRUE)
		      return(NULL)
		    },
		    error = function(e) {
		      write(file,file =paste0('Error_inJSON_',Sys.Date(),'.txt'),append = TRUE)
		      return(NULL)
		    })
		#method = r1$result$detail$applied_method
		if(is.null(r1))
		  next
		#method = r1$result$detail$applied_method


		###### only MM's
		#if (!is.null(method) && method  %in% 'MM') {
			message(paste(i, '|',  end , ':', file))
			x = c(
				unlist(r0[1, 1:19]),
				unlist0(r1$Result$Details$`Response type`) ,
				unlist0(r1$Result$Details$`Applied method`),
				################# Window parameters
				# Have you changed that for new structure of l and k output???
				SelectWindowingParameters(object = r1$Result$Details)                 ,
				################# VectorOutput Results
				SelectAnalysis(r1$Result$`Vector output`$`Normal result`)             ,
				SelectAnalysis(r1$Result$`Vector output`$`Windowed result`)           ,
				SelectAnalysis(r1$Result$`Vector output`$`Full model result`)         ,
				SelectAnalysis(r1$Result$`Vector output`$`Full model windowed result`),
				################# Other results
				SelectOthers(r1$Result$Details)                                      ,
				################ Ignorome/Reference/Behaviour
				SelectMISC(r0)                                                       ,
				##### Variation in response
				unlist0(r1$Result$Details$variation_in_respone_after_preprocess[1])   ,
				unlist0(r1$result$details$variation_in_respone_before_preprocess[1])  ,
				##### Pvals
				unlist0(r1$Result$`Vector output`$`Normal result`$`Genotype p-value`)      ,
				unlist0(r1$Result$`Vector output`$`Windowed result`$`Genotype p-value`)    ,
				unlist0(r1$Result$`Vector output`$`Full model result`$`Genotype p-value`)  ,
				unlist0(r1$Result$`Vector output`$`Full model windowed result`$`Genotype p-value`),
				##### MP TERM
				DRrequiredAgeing:::StratifiedMPTerms(rN),
				DRrequiredAgeing:::StratifiedMPTerms(rW),
				##### URL
				MakeURL(r0, r1),
				tableCount(
				  Gen = r1$Result$Details$Original_biological_sample_group,
				  Sex = r1$Result$Details$Original_sex
				)
			)
			write(
			  x = paste(x, collapse = '\t'),
			  file = DRrequiredAgeing:::file.path0(
			    paste0('./resultF/', r1$Result$Details$`Response type`, '/', ofname),
			    create = TRUE,
			    check  = FALSE,
			    IncludedFileName = TRUE
			  ),
			  append = TRUE,
			  ncolumns = 10 ^ 4
			)
		#}
	}
}

ignore.my.name = f(start =  as.numeric(args[1]), end = as.numeric(args[2]),file = args[3])






# After getting the table
# library(plyr)
# TotalN = length(OnlyChanges_P_Value_results_2018_08_08_10_59_03$Procedure_group)
# cdata <-
# 	ddply(
# 		OnlyChanges_P_Value_results_2018_08_08_10_59_03,
# 		c("Centre","Procedure_group"),
# 		summarise,
# 		N        = length(`Gain (windowing significant but not normal)`),
# 		CentreContribution = paste0(round(N / TotalN * 100,1), '%'),
# 		nGain    = sum(`Gain (windowing significant but not normal)`),
# 		nLoss    = sum(`Loss (Normal significant but not windowing)`),
# 		#GainLossRatio = round(nGain / nLoss, 2),
# 		Gain2TotalInCentre    = paste0(round(nGain / N * 100), '%'),
# 		Loss2TotalInCentre    = paste0(round(nLoss / N * 100), '%')
#
# 	)
# cdata[order(cdata$N),]
