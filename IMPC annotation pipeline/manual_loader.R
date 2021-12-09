file = commandArgs(trailingOnly = TRUE)
changeRpackageDirectory = function(path = '~/DRs/R/packages') {
  v = paste(
    R.version$major,
    gsub(
      pattern = '.',
      replacement = '_',
      R.version$minor,
      fixed = TRUE
    ),
    sep = '_',
    collapse = '_'
  )
  wdirc = file.path(path, v)
  if (!dir.exists(wdirc))
    dir.create(wdirc,
               showWarnings = FALSE,
               recursive = TRUE)
  .libPaths(new = wdirc)
  message(' => new package path set to: ', wdirc,'. Please add this to the .bash_profile under R_LIBS_USER  environment variable. [file = manual_loader.R]')
}
changeRpackageDirectory()
library('DRrequiredAgeing')
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


# file = 'https://www.ebi.ac.uk/~hamedhm/windowing/DR11/jobs/Results_DR11V2OpenStats/RBRC/IMPC_CSD/IMPC_CSD_050_002/RIKEN_Fbn2_A06_/homozygote/978482d700f332623fc65f3e729dea93/output_Successful.tsv'
# file ='https://www.ebi.ac.uk/~hamedhm/windowing/DR11/jobs/Results_DR11V2OpenStats/ICS/ESLIM_008/ESLIM_008_001_014/EPD0060_2_H09/heterozygote/d41d8cd98f00b204e9800998ecf8427e/output_Successful.tsv'
# file = 'https://www.ebi.ac.uk/~hamedhm/windowing/DR12/jobs/Results_DR12V1OpenStats/ICS/ESLIM_008/ESLIM_008_001_014/EPD0060_2_H09/heterozygote/d41d8cd98f00b204e9800998ecf8427e/output_Successful.tsv'
# download.file(url = file, destfile = 'delme.tsv')
# file = 'delme.tsv'
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
                                              level = .0001)
    rW = DRrequiredAgeing:::annotationChooser(
      statpacket = rN$statpacket,
      level = .0001,
      resultKey = 'Windowed result',
      TermKey = 'WMPTERM'
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
    status  = DRrequiredAgeing:::Write2Postg(df = df[1, ], host =  "hh-yoda-05-01",tablename = 'db_xxxxxx')
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
