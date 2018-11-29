grepl0 = function(pattern,x,...){
  s = rep(FALSE,length(x))
  for (pat in pattern){
    s = s | grepl(pattern = pat,x = x,...)
  }
  return(s)
}

pie2 = function(x,...){
  x= table(x)
  piepercent<- paste(names(x),' [#',x,'] ',round(100*x/sum(x), 1),'%',sep = '')
  pie(x, labels = piepercent, ...)
}