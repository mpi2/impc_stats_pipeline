file = commandArgs(trailingOnly = TRUE)
library('RPostgreSQL')
library('data.table')
########################### Annotation pipeline #################################
##############################
library(data.table)
library(jsonlite)
library(rlist)
library(Tmisc)
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
tablename = configlist$tablename
dbname = configlist$dbname
port = configlist$port
user = configlist$user
password = configlist$password
level = configlist$level
pplevel = configlist$pplevel
###########################################

flist = readLines(con = file[1])
lflist = length(flist)
id = 1
rnd <-
  DRrequiredAgeing:::RandomRegardSeed(
    n = 1,
    decimal = 0,
    stringOutput = 1,
    max = 99999,
    round = 1
  )
for (i in 1:lflist) {
  cat('\r', i, '/', lflist)
  status = FALSE
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
                                              pplevel = pplevel,
                                              mp_chooser_file = mp_chooser_file)
    rW = DRrequiredAgeing:::annotationChooser(
      statpacket = rN$statpacket,
      level = level,
      pplevel = pplevel,
      resultKey = 'Windowed result',
      TermKey = 'WMPTERM',
      mp_chooser_file = mp_chooser_file
    )

    df = rW$statpacket
    ###################
    df = cbind(as.numeric(paste0(
      shuffle(),
      format(Sys.time(), "%H%S")
    )), df)


    names(df)  =  c(
      "id",
      "analysisTime",
      "status",
      "procedure_group",
      "procedure_stable_id",
      "procedure_name",
      "parameter_stable_id",
      "parameter_name",
      "phenotyping_center",
      "allele_symbol",
      "allele_name",
      "allele_accession_id",
      "gene_symbol",
      "gene_accession_id",
      "pipeline_name",
      "pipeline_stable_id",
      "strain_accession_id",
      "metadata_group",
      "zygosity",
      "colony_id",
      "statpacket"
    )
    status  = DRrequiredAgeing:::Write2Postg(
      df = df[1, ],
      host =  host,
      tablename = tablename,
      dbname = dbname,
      port = port,
      user = user,
      password = password
    )
    id = id + 1
    #print(file)
  }
  if (!dir.exists("log")) {
    dir.create("log")
  }
  write(
    paste(status, file, sep = "\t-->", collapse = "\t-->"),
    file = paste0("./log/", rnd, "__Log.log"),
    append = TRUE
  )
  gc()
}
