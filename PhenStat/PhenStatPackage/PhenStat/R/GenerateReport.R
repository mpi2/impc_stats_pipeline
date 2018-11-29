# Main function for ReportGenerating function
PhenStatReport = function(PhenlistObject                    ,
													depVariable   = NULL              ,
													other.response = NULL             ,
													update        = TRUE              ,
													Gene.Symbol   = NULL              ,
													Response.name = NULL              ,
													destination   = NULL              ,
													reportTitle   = 'Extended Statistical Report' ,
													DataRelease   = NULL              ,
													Showsource    = FALSE             ,
													open          = FALSE             ,
													clean         = TRUE              ,
													verbos        = FALSE             ,
													...) {
	requireNamespace('knitr')
	requireNamespace('tools')
	requireNamespace('pingr')
	### Output parameters
	knitr::opts_chunk$set(echo = Showsource)
	knitr::opts_chunk$set(results = 'markup')
	knitr::opts_chunk$set(size    = 'small')
	knitr::opts_chunk$set(background = '#fffae6')
	knitr::opts_chunk$set(tidy = TRUE)
	knitr::opts_chunk$set(cache = FALSE)
	#knitr::opts_chunk$set(fig.pos = 'htp')
	knitr::opts_chunk$set(fig.height = 6)
	knitr::opts_knit$set(kable.force.latex = TRUE)
	knitr::opts_chunk$set(fig.align = 'center')
	knitr::opts_chunk$set(autodep = TRUE)
	knitr::opts_chunk$set(dpi = 60)
	knitr::opts_knit$set(width = 95)
	options(width = 95)

	# load(file = file.path(
	#  system.file(package = "PhenStat", mustWork = TRUE),
	#  'report',
	#  'GeCoSy.Rdata'
	#))

	curWD = getwd()
	if (is.null(depVariable)) {
		stop(
			paste(
				"\nPlease specify dependent variable, ",
				"for example: depVariable='Lean.Mass'.\n",
				sep = ""
			)
		)
	}
	### Create a temp directory
	tmp = MakeUniqueString(1)
	if (is.null(destination)) {
		tdes = b2f(file.path(tempdir(), tmp))
	} else {
		tdes = b2f(file.path(dirname(destination), tmp))
	}
	cat(paste("\n Creating destination: ", tdes, "\n"))
	dir.create(tdes, showWarnings = verbos, recursive = TRUE)
	tmp.dir = tdes
	###
	setwd(tmp.dir)
	dir.create(file.path(tmp.dir, "myFigs"),
						 showWarnings = verbos,
						 recursive = TRUE)
	knitr::opts_knit$set(unnamed.chunk.label = paste(s2s(depVariable), "_myRNDid_",
																									 round(runif(1) * 10 ^ 6), sep = "_"))
	knitr::opts_chunk$set(fig.path = file.path(tmp.dir, "//myFigs//"))



	### check the version of the report
	Baseurl = 'https://www.mousephenotype.org/sites/beta.mousephenotype.org/files/mousephenotype_files/'
	if (update) {
		if (pingr::is_online()) {
			current.version = as.numeric(gsub("\\.", "", packageVersion("PhenStat")))
			Allversion = read.table(paste(Baseurl, '/version.txt', sep = ''))
			versions = apply(Allversion, 2, function(x) {
				as.numeric(gsub("\\.", "", x))
			})
			version = tail(sort(versions), 1)
		} else {
			message(
				'Please check your connection! : Error in downloading the new version of the report.
				\n The default version of the report is used.
				'
			)
			update          = FALSE
			version = current.version = as.numeric(gsub("\\.", "", packageVersion("PhenStat")))
		}
	} else{
		version = current.version = as.numeric(gsub("\\.", "", packageVersion("PhenStat")))
	}

	####
	if (update && pingr::is_online()) {
		if (version >= current.version) {
			cat(
				"\n New version of the report found!
				\n Downloading the new version of the report ...\n"
			)
			if (!download.file(
				url = paste(Baseurl, "/report_",
										version,
										".txt",
										sep = ""),
				destfile = file.path(tmp.dir, paste("/report_",
																						version, ".Rnw", sep = "")),
				quiet    = !verbos,
				cacheOK  = FALSE
			) &&
			!download.file(
				url = paste(Baseurl, "/References.txt", sep = ''),
				destfile = file.path(tmp.dir, "maina.bib"),
				quiet    = !verbos,
				cacheOK  = FALSE
			)) {
				filename = paste("report_", version, ".Rnw", sep = "")
				report_file = file.path(tmp.dir, filename, sep = "")
			} else{
				message('Error in downloading the report!')
				message('The default version of the report is used!')
				report_file = file.path(system.file(package = "PhenStat", mustWork = TRUE),
																"report",
																"report.Rnw")
				file.copy(from = report_file,
									to   = file.path(tmp.dir),
									overwrite = TRUE)
				file.copy(
					from = file.path(
						system.file(package = "PhenStat", mustWork = TRUE),
						"/report/maina.bib"
					),
					to = file.path(tmp.dir),
					overwrite = TRUE
				)
				report_file = file.path(tmp.dir, 'report.Rnw')
			}
		} else{
			message('The default version of the report is used!')
			report_file = file.path(system.file(package = "PhenStat", mustWork = TRUE),
															"report",
															"report.Rnw")
			file.copy(from = report_file,
								to   = file.path(tmp.dir),
								overwrite = TRUE)
			file.copy(
				from = file.path(
					system.file(package = "PhenStat", mustWork = TRUE),
					"/report/maina.bib"
				),
				to = file.path(tmp.dir),
				overwrite = TRUE
			)
			report_file = file.path(tmp.dir, 'report.Rnw')
		}
	} else{
		dst = file.path(system.file(package = "PhenStat", mustWork = TRUE),
										"report")
		filename = "report.Rnw"
		file.copy(
			from = file.path(dst, filename),
			to   = file.path(tmp.dir, filename),
			overwrite = TRUE
		)
		file.copy(
			from = file.path(
				system.file(package = "PhenStat", mustWork = TRUE),
				"/report/maina.bib"
			),
			to = file.path(tmp.dir, "maina.bib"),
			overwrite = TRUE
		)
		report_file = file.path(tmp.dir, filename)
	}
	######
	save(PhenlistObject, file = paste(report_file, '.Rdata', sep = ''))
	PhenListInput             = PhenlistObject
	texfile = paste(report_file, "_", s2s(depVariable), ".tex", sep = "")
	pdffile = paste(report_file, "_", s2s(depVariable), ".pdf", sep = "")
	###### 1.
	cat("\n 1/3. Processing the report ...\n")
	tex = knitr::knit(
		input = report_file,
		output = texfile,
		tangle = FALSE,
		quiet = !verbos,
		encoding = "UTF-8"
	)
	##### 2.
	cat("\n 2/3. Processing the pdf file ...\n")
	tools::texi2pdf(
		file = tex,
		clean = TRUE,
		quiet = !verbos,
		index = TRUE,
		texinputs = tmp.dir
	)
	#### 3.
	if (!is.null(destination)) {
		file.copy(from = pdffile,
							to = destination,
							overwrite = TRUE)
		pdffile = destination
	}

	if (clean && dir.exists(tmp.dir)) {
		if (is.null(destination)) {
			file.copy(
				from = pdffile,
				to = curWD ,
				overwrite = TRUE,
				recursive = TRUE
			)
			pdffile = file.path(curWD, basename(pdffile))
		}
		tex = NULL
		unlink(tmp.dir, recursive = TRUE, force = TRUE)
	}

	cat(paste("\n 3/3. pdf file is stored in ", pdffile, sep = ""), "\n")
	if (open)
		system(paste("open \"", pdffile, "\"", sep = ""))

	cat("\n Done! \n")
	setwd(curWD)

	return(
		list(
			PhenlistObject = PhenlistObject,
			update = update,
			depVariable = depVariable,
			other.response = other.response,
			destination = destination,
			verbos = verbos,
			open = open,
			Showsource = Showsource,
			version = version,
			texfile = tex,
			pdffile = pdffile,
			output.dir = tmp.dir
		)
	)
}
