WhatIsNew = function() {
	cat(
		'\n
	What is new in Version', as.character(packageVersion("PhenStat")), ':\n
	1. Due to the complexity of the method, Soft windowing function is REMOVED,
		see https://github.com/cran/SmoothWin for the new implementation of the method\n
	2. VectorOutput now lets user defined variables in the output\n
	3. Several improvements on VectorOutput including summary statistics and so on\n
	4. Bug fixed and minor improvements\n
	'
	)

}
