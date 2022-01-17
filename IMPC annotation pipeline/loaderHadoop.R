file = commandArgs(trailingOnly = TRUE)
library('RPostgreSQL')
library('data.table')
########################### Annotation pipeline #################################
##############################
library(data.table)
library(jsonlite)
library(rlist)
library(Tmisc)
library(rwebhdfs)
##############################
shuffle = function (x = as.numeric(Sys.time()) * 100000,
                    replace = TRUE,
                    length = 5) {
  r = as.numeric(unlist(strsplit(as.character(x), "")))
  r = sample(r, length, replace = replace)
  return(paste0(r, collapse = ''))
}


###########################################
load('config.Rdata')
mp_chooser_file = configlist$mp_chooser_file
host =  configlist$host
path = configlist$path
prefix = configlist$prefix
port = configlist$port
user = configlist$user
password = configlist$password
level = configlist$level
rrlevel = configlist$rrlevel
###########################################

flist = readLines(con = file[1])
lflist = length(flist)
id = 1
statpackets_out = NULL
for (i in 1:lflist) {
  cat('\r', i, '/', lflist)
  file = flist[i]
  cat("\n", i, "/", lflist, "~>", file, "")
  if (file.exists(file) &&
      (grepl(pattern = 'NotProcessed', x = file) ||
       grepl(pattern = 'Successful'  , x = file))) {
    df = data.table::fread(
      file = file,
      header = FALSE,
      sep = '\t',
      quote = "",
      stringsAsFactors = FALSE
    )
    if (ncol(df) != 20)
      next
    ###################
    rN = DRrequiredAgeing:::annotationChooser(statpacket = df,
                                              level = level,
                                              rrlevel = rrlevel,
                                              mp_chooser_file = mp_chooser_file)
    rW = DRrequiredAgeing:::annotationChooser(
      statpacket = rN$statpacket,
      level = level,
      rrlevel = rrlevel,
      resultKey = 'Windowed result',
      TermKey = 'WMPTERM',
      mp_chooser_file = mp_chooser_file
    )
    statpackets_out = c(statpackets_out,rW$statpacket)
  }
  if (!dir.exists("log")) {
    dir.create("log")
  }
  write(
    paste(file, sep = "\t-->", collapse = "\t-->"),
    file = paste0("./log/", basename(file[1]), "__Log.log"),
    append = TRUE
  )
  
  hdfs <- webhdfs(namenode_host = host, namenode_port = port,hdfs_username =  user)
  write_file(hdfs, file.path(path, format(Sys.time(), '%a%b%d_%Y'), paste0(basename(file[1], '_.json'))))
  
  gc()
}
