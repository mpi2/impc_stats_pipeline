args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Two arguments are required: the file path and the mp_chooser_file.")
}
file <- args[1]
mp_chooser_file <- args[2]

# Load necessary libraries
library(data.table)
library(jsonlite)
library(rlist)
library(Tmisc)
library(rwebhdfs)

# Set levels
level <- .0001
rrlevel <- .0001

# Start annotation pipeline
today <- format(Sys.time(), "%d%m%Y")
flist <- readLines(con = file)
lflist <- length(flist)

# Store StatPackets temporary
if (!dir.exists("tmp")) {
  dir.create("tmp")
}

tmplocalfile <- file.path('tmp', paste0(basename(file), '_.statpackets'))

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
  }
}
