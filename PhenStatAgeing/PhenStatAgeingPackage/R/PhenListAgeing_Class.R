setClass(
	'PhenListAgeing',
	representation(
		datasetPL    = 'data.frame'    ,
		refGenotype  = 'character'     ,
		testGenotype = 'character'     ,
		hemiGenotype = 'character'     ,
		dataset.colname.batch    = 'character'     ,
		dataset.colname.genotype = 'character'     ,
		dataset.colname.sex      = 'character'     ,
		dataset.colname.weight   = 'character'     ,
		dataset.values.missingValue = 'character'  ,
		dataset.values.male      = 'character'     ,
		dataset.values.female    = 'character'     ,
		dataset.clean = 'logical'             ,
		datasetUNF = 'data.frame'             ,
		reserve1   = 'list'                   ,
		reserve2   = 'data.frame',
		reserve3   = 'character'
	)
)
