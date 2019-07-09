plot.NULL = function(x, ...) {
	message0('No plot available for a NULL object')
}

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
	LowRes  = pastedot('Low' , Labels$response, 'Genotype')
	HighRes = pastedot('High', Labels$response, 'Genotype')
	
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
		l1 = LowRes,
		l2 = Labels$Sex$Sex,
		main = main,
		...
	)
	plotFERR (
		x = x,
		l1 = HighRes,
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
	requireNamespace("car")
	fm        = x$output$Final.Model
	formula   = formula(fm)
	transData = applyFormulaToData(formula = formula, getData(fm))
	
	if (!is.null(x$messages) || is.null(x)) {
		message0('Due to error(s), no plot available')
		message0(x$messages)
		stop()
	}
	p = par()
	par(ask = ask, mfrow = mfrow)
	
	predR            = predict(fm)
	residR           = resid(fm)
	residShapiroTest = shapiro.test0(residR)
	plot(
		predR ,
		residR,
		xlab = 'Fitted values',
		ylab = 'Residuals'    ,
		main = main           ,
		#ellipse = TRUE        ,
		#reset.par = FALSE     ,
		...
	)
	abline(h = 0, lwd = 3, lty = 2)
	densityPlot(
		residR,
		xlab = ifelse(
			is.null(residShapiroTest$'P-value'),
			'Residuals'              ,
			paste0(
				'Residual - '          ,
				residShapiroTest$'Test',
				' p-value = '          ,
				round(residShapiroTest$'P-value', 8)
			)
		),
		main = paste0(main, ': Density of the residuals'),
		...
	)
	qqPlot(
		as.vector(residR),
		ylab = 'Residuals',
		main = paste0(main, ': Normal Q-Q of the residuals'),
		grid = FALSE ,
		col.lines = 1,
		...
	)
	#qqline(residR, ...)
	
	densityPlot(transData$data[, transData$names[1]],
			 main = 'Density of the response',
			 xlab = transData$names[1],
			 ...)
	par(ask = p$ask, mfrow = p$mfrow)
}
