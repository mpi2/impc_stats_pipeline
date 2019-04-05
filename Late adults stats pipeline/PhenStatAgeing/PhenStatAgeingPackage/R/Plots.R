plotFERR = function(x, l1, l2, main, ...) {
	if (!is.null(x$output$SplitModels[[l1]][[l2]]$table))
		mosaicplot(x$output$SplitModels[[l1]][[l2]]$table, main = main, ...)
	else
		message0('Did not find [',
						 l1,
						 ' X ',
						 l2,
						 '] table')
}

###############################################
# Plot RR
###############################################
plot.PhenStatAgeingRR = function(x,
																 main = 'Mosaic plot',
																 ask = FALSE         ,
																 mfrow = c(2, 2)     ,
																 ...) {
	if (!is.null(x$messages) || is.null(x)) {
		message0('Due to error(s), no plot available')
		message0(x$messages)
		stop()
	}
	p = par()
	par(ask = ask, mfrow = mfrow)
	
	Labels  = PhenListAgeingLevels(x)
	LowRes  = pasteUnderscore('Low' , Labels$response, 'discretised')
	HighRes = pasteUnderscore('High', Labels$response, 'discretised')
	
	plotFERR (
		x = x,
		l1 = LowRes,
		l2 = Labels$Genotype$Genotype,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = HighRes,
		l2 = Labels$Genotype$Genotype,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = pasteUnderscore('Low', Labels$Genotype$Genotype),
		l2 = Labels$Sex$Sex,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = pasteUnderscore('Low', Labels$Genotype$Genotype),
		l2 = Labels$LifeStage$LifeStage,
		main = main,
		...
	)
	
	par(ask = p$ask, mfrow = p$mfrow)
}

###############################################
# Plot FE
###############################################
plot.PhenStatAgeingFE = function(x,
																 main = 'Mosaic plot',
																 ask = FALSE         ,
																 mfrow = c(2, 2)     ,
																 ...) {
	if (!is.null(x$messages) || is.null(x)) {
		message0('Due to error(s), no plot available')
		message0(x$messages)
		stop()
	}
	p = par()
	par(ask = ask, mfrow = mfrow)
	
	Labels = PhenListAgeingLevels(x)
	plotFERR (
		x = x,
		l1 = Labels$response,
		l2 = Labels$Genotype$Genotype,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = Labels$response,
		l2 = Labels$Sex$Sex,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = Labels$Genotype$Genotype,
		l2 = Labels$Sex$Sex,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = Labels$response,
		l2 = Labels$LifeStage$LifeStage,
		main = main,
		...
	)
	par(ask = p$ask, mfrow = p$mfrow)
}

###############################################
# Plot MM
###############################################
plot.PhenStatAgeingMM = function (x                   ,
																	main = 'Final Model',
																	ask = FALSE         ,
																	mfrow = c(2, 2)     ,
																	...) {
	if (!is.null(x$messages) || is.null(x)) {
		message0('Due to error(s), no plot available')
		message0(x$messages)
		stop()
	}
	p = par()
	par(ask = ask, mfrow = mfrow)
	
	predR  = predict(x$output$Final.Model)
	residR = resid(x$output$Final.Model)
	plot(predR ,
			 residR,
			 xlab = 'Fitted values',
			 ylab = 'Residuals'    ,
			 main = main           ,
			 ...)
	abline(h = 0, lwd = 3, lty = 2)
	hist(
		residR,
		xlab = 'Residuals',
		main = paste0(main, ': Histogram of residuals'),
		probability = TRUE,
		...
	)
	qqnorm(residR, main = paste0(main, ': Normal Q-Q plot of residuals'), ...)
	qqline(residR, ...)
	plot(density(getData(x$output$Final.Model)[, x$input$depVariable]), main = 'Density of the response', ...)
	par(ask = p$ask, mfrow = p$mfrow)
}
