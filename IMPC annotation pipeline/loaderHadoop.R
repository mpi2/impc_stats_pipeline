orgfile = commandArgs(trailingOnly = TRUE)
file = orgfile
library('data.table')
########################### Annotation pipeline #################################
##############################
library(data.table)
library(jsonlite)
library(rlist)
library(Tmisc)
library(rwebhdfs)
###########################################
load('configHadoop.Rdata')
mp_chooser_file = configlist$mp_chooser_file
host =  configlist$host
path = configlist$path
prefix = configlist$prefix
port = configlist$port
user = configlist$user
password = configlist$password
level = configlist$level
rrlevel = configlist$rrlevel
transfer = configlist$transfer
###########################################
today = format(Sys.time(), '%d%m%Y')
flist = readLines(con = file[1])
lflist = length(flist)


# store StatPackets temporary
if (!dir.exists("tmp")) {
  dir.create("tmp")
}

tmplocalfile    =  file.path('tmp', paste0(basename(orgfile[1]), '_.statpackets'))


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
    if (ncol(df) != 20 ||
        nrow(df) > 1) {
      message('file ignored (!=20 columns): ', file)
      next
    }
    ###################
    rN = DRrequiredAgeing:::annotationChooser(
      statpacket = df,
      level = level,
      rrlevel = rrlevel,
      mp_chooser_file = mp_chooser_file
    )
    rW = DRrequiredAgeing:::annotationChooser(
      statpacket = rN$statpacket,
      level = level,
      rrlevel = rrlevel,
      resultKey = 'Windowed result',
      TermKey = 'WMPTERM',
      mp_chooser_file = mp_chooser_file
    )

    write(paste0(as.character(rW$statpacket$V20), collapse = ''),
          file = tmplocalfile,
          append = TRUE)
    #statpackets_out = c(statpackets_out, rW$statpacket)
    #writeLines(statpackets_out, con = tmplocalfile)
  }
}

# statpackets need to be stored as characters
# statpackets_out = as.character(statpackets_out)

if (transfer) {
  # Prepare and transfer files to hadoop
  hadoopPath = file.path(path,
                         prefix,
                         today,
                         paste0(basename(orgfile[1]), '_.statpackets'))

  hdfs <-
    webhdfs(
      namenode_host = host,
      namenode_port = port,
      hdfs_username =  user
    )
  rwebhdfs::mkdir(hdfs, dirname(hadoopPath))

  transfered = rwebhdfs::write_file(
    fs = hdfs,
    targetPath = hadoopPath,
    srcPath = tmplocalfile,
    sizeWarn = 10 ^ 16,
    append = FALSE,
    overwrite = TRUE
  )
  gc()
  if (transfered) {
    message('Transfer successful.')
    message('Compressing the temp statpacket file.')
    system(command = paste0('gzip ', tmplocalfile),
           wait = TRUE)
    message('Done!')
  }  else{
    stop('Transfered not successful!')
  }
}
gc()
