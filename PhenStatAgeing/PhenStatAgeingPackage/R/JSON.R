toJSONI = function(object, ...) {
	object$call$data = 	object$data = 	object$OriginalDataSet = object$object$data =
		object$object$OriginalDataSet = NULL
	#class(object) = NULL
	object = lapply(object, replace_sub_dataframes)
	r = jsonlite::toJSON(
		object,
		pretty = TRUE,
		auto_unbox = TRUE,
		factor = 'string',
		null  = 'null',
		simplifyVector = 1,
		simplifyDataFrame = 1,
		simplifyMatrix = 1,
		na   = 'null',
		dataframe = 'columns',
		matrix = 'columnmajor',
		force = TRUE,
		...
	)
	r <- gsub(pattern = "^\\[",
						replacement = "",
						x = r)
	r <- gsub(pattern = "\\]$",
						replacement = "",
						x = r)
	return(r)
}


replace_sub_dataframes <- function(x)
{
	if (any(class(x) == "list"))
	{
		x_copy = x
		attrs = setdiff(names(attributes(x)), "names")
		x = lapply(x, replace_sub_dataframes)
		if (length(attrs) > 0)
		{
			for (i in 1:length(attrs))
			{
				attr(x, attrs[i]) <- replace_sub_dataframes(attr(x_copy, attrs[i]))
			}
		}
		return(x)
	}
	else
	{
		if (any(class(x) == "matrix"))
			return(as.data.frame(x))
		else
			return(x)
	}
}
