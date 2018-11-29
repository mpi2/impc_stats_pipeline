PhenSmoothWin = function(testDatasetObj ,
												 l = 4:5,
												 k = 1:2,
												 min.obs   = 50,
												 criteria  = 'AICc',
												 method    = 'enet',
												 threshold = 10 ^ -18,
												 plot      = FALSE) {
	message(
		'This function is deactivated in this version of PhenStat.\n Please refer to the SmoothWin repository on github for more information.\n https://github.com/cran/SmoothWin'
	)
	# #Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.21/bin/gswin64.exe")
	# requireNamespace("SmoothWin")
	# min.obs = ceiling(min.obs)
	# if (testDatasetObj@method != 'MM')
	# 	stop(
	# 		'This method is only available for the data that are analysed by Linear Mixed Model in PhenStat'
	# 	)
	#
	# k = sort(k[k > 0 & !is.na(k)], decreasing = TRUE)
	# l = sort(l[l > 0 & !is.na(l)])
	# #if (length(k) <= 1 | length(l) <= 1)
	# #print('l and k must have at least 2 values each!')
	#
	# t = as.numeric(as.Date(testDatasetObj@analysedDataset[,
	# 																											which(names(testDatasetObj@analysedDataset) == 'Batch')]))
	#
	#
	# y = testDatasetObj@analysedDataset[, testDatasetObj@depVariable]
	# x = testDatasetObj@analysedDataset[, which(names(testDatasetObj@analysedDataset) != testDatasetObj@depVariable)]
	# m = which(testDatasetObj@analysedDataset$Genotype == testDatasetObj@testGenotype)
	#
	# if (length(unique(t[m])) > 10)
	# 	message('More than 10 modes detected. The entire procedure may take a long time!')
	# if (length(m) > min.obs) {
	# 	stop('min.obs is less than the number of mutants!')
	# } else{
	# 	cat(
	# 		'\n Minimum observations in the model : (#mut =',
	# 		length(m),
	# 		', #con =',
	# 		min.obs,
	# 		') =',
	# 		min.obs + length(m)
	# 	)
	# }
	# ### 1. Determining l
	# cat('\n 1|3 Searching for l ...\n')
	# rl = PhengridSearch(
	# 	PhenStatTestDataObject = testDatasetObj,
	# 	t = t,
	# 	x = x,
	# 	y = y,
	# 	m = t[m],
	# 	l = l,
	# 	k = max(k),
	# 	mutInd    = m,
	# 	plot      = plot,
	# 	threshold = threshold
	# )
	# ####
	# dcpt        = suppressMessages(SmoothWin::scale21(rl$output$esd) * 100)
	# if (length(unique(dcpt)) > 1) {
	# 	cpt         = SmoothWin::penCPD(
	# 		dcpt,
	# 		threshold = threshold,
	# 		method = method,
	# 		criteria = criteria,
	# 		plot  = plot,
	# 		main = 'l - CPD solution path'
	# 	)
	# } else{
	# 	cpt      = length(dcpt)
	# }
	# ####
	# finall      = rl$output$l[cpt][which(rl$output$Obs.in.Interval[cpt] >=
	# 																		 	(min.obs + length(m)))][1]
	# if (is.na(finall)) {
	# 	cat('\n An optimal l is not found. Max l will be used.')
	# 	finall = max(l)
	# }
	#
	# ### 2. Determining k
	# cat('\n 2|3 Searching for k ... \n')
	# rk    = PhengridSearch(
	# 	PhenStatTestDataObject = testDatasetObj,
	# 	t = t,
	# 	x = x,
	# 	y = y,
	# 	m = t[m],
	# 	l = finall,
	# 	k = k,
	# 	plot = 0,
	# 	mutInd = m,
	# 	threshold = threshold
	# )
	#
	# dcptk = suppressMessages(SmoothWin::scale21(rk$output$esd) * 100)
	# if (length(unique(dcptk)) > 1) {
	# 	cptk  = SmoothWin::penCPD(
	# 		x = dcptk,
	# 		threshold = threshold,
	# 		method = method,
	# 		criteria = criteria,
	# 		plot  = plot,
	# 		main = 'k - CPD solution path'
	# 	)
	# } else{
	# 	cptk = 1
	# }
	# finalk = rk$output$k[cptk][which(rk$output$Obs.in.Interval[cptk] >=
	# 																 	(min.obs + length(m)))][1]
	# if (is.na(finalk)) {
	# 	cat('\n An optimal k is not found. Max k will be used.')
	# 	finalk = max(k)
	# }
	#
	#
	#
	# ##### Plots
	# if (plot && !is.na(finall) && length(l) > 1) {
	# 	plot(
	# 		rl$output$l,
	# 		SmoothWin::scale21(rl$output$esd),
	# 		ylim = c(0, 3),
	# 		xaxt = 'n',
	# 		ylab = 'Standardized sd',
	# 		xlab = 'l',
	# 		main = paste('Final l = ', round(finall, 5), sep = '')
	# 	)
	#
	# 	abline(v = rl$output$l[cpt],
	# 				 col = 'grey',
	# 				 lty = 3)
	# 	axis(
	# 		side = 1,
	# 		at = rl$output$l[cpt],
	# 		labels = round(rl$output$l[cpt], 3),
	# 		las = 3
	# 	)
	# 	abline(
	# 		v = finall,
	# 		col = 2,
	# 		lwd = 2,
	# 		lty = 3
	# 	)
	# 	text(
	# 		x      = rl$output$l[cpt],
	# 		y      = 2 + (cos(rl$output$Ind[cpt] / 2 / pi)),
	# 		col    = 2 + ((-1) ^ rl$output$Ind[cpt]) ,
	# 		labels = rl$output$Obs.in.Interval[cpt],
	# 		cex    = .5
	# 	)
	#
	# }
	# if (plot && !is.na(finalk) && length(k) > 1) {
	# 	plot(k,
	# 			 rk$output$BIC,
	# 			 xlab = 'k',
	# 			 ylab = 'BIC',
	# 			 main = 'k - BIC')
	# 	abline(v = finalk)
	#
	# 	plot(
	# 		rk$output$k,
	# 		SmoothWin::scale21(rk$output$esd),
	# 		ylim = c(0, 3),
	# 		xaxt = 'n',
	# 		ylab = 'Standardized sd',
	# 		xlab = 'k',
	# 		main = paste('Final k = ', round(finalk, 5), sep = '')
	# 	)
	#
	# 	abline(v = rk$output$k[cptk],
	# 				 col = 'grey',
	# 				 lty = 3)
	# 	axis(
	# 		side = 1,
	# 		at = rk$output$k[cptk],
	# 		labels = round(rk$output$k[cptk], 3),
	# 		las = 3
	# 	)
	# 	abline(
	# 		v = finalk,
	# 		col = 4,
	# 		lwd = 2,
	# 		lty = 3
	# 	)
	# 	text(
	# 		x      = rk$output$k[cptk],
	# 		y      = 2 + (cos(rk$output$Ind[cptk] / 2 / pi)),
	# 		col    = 2 + ((-1) ^ rk$output$Ind[cptk]) ,
	# 		labels = rk$output$Obs.in.Interval[cptk],
	# 		cex    = .5
	# 	)
	# }
	#
	#
	# ##### final model
	# cat('\n 3|3 Forming the final model \n')
	# finalr = PhengridSearch(
	# 	PhenStatTestDataObject = testDatasetObj,
	# 	t = t,
	# 	x = x,
	# 	y = y,
	# 	m = t[m],
	# 	l = finall,
	# 	k = finalk ,
	# 	mutInd = m,
	# 	plot = plot,
	# 	threshold = threshold
	# )
	# finalr$weights = finalr$weights[-1] #Remove index from the final weights
	# return(
	# 	list(
	# 		final.k = finalk,
	# 		final.l = finall,
	# 		finalModel = finalr,
	# 		model.l = rl,
	# 		model.k = rk,
	# 		x = x,
	# 		y = y,
	# 		t = t,
	# 		l = l,
	# 		k = k,
	# 		criteria  = criteria,
	# 		method    = method,
	# 		min.obs   = min.obs,
	# 		threshold = threshold,
	# 		mutInd    = m,
	# 		plot      = plot
	# 	)
	# )
}
