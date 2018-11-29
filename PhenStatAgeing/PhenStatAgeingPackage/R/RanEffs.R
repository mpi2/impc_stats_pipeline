rndProce = function(procedure) {
	if (length(procedure) < 1 ||	is.null(procedure) ||
			!(procedure %in% lop()) || length(procedure) > 1)
		stop ('\n=> Error in inputing the procedure symbol \n  => Current input : ',
					procedure,
					'\n')
	##############################################################
	if (procedure == 'TYPICAL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ABR') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'ACS') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'ALZ') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'BLK') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'BWT') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'CAL') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'CBC') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'CHL') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'CSD') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'DXA') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'ECG') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'ECH') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'ELZ') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVL') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVM') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVO') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVP') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'EYE') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'FER') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEL') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEM') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEO') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEP') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPL') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPM') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPO') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'GRS') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'HEM') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'HIS') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'HWT') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'IMM') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'INS') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'IPG') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'OFD') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'PAT') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'VIA') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	} else if (procedure == 'XRY') {
		random = reformulate  (' 1 |  Batch / Age', response = NULL, intercept = TRUE)
	}
	
	return(random)
}

