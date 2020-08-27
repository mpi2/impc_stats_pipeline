rm(list = ls(all = TRUE))
# B6N strains can be compared to C57BL/6NJ` BUT NOT`  C57BL/6J
# B6N Bottom line, our wildtype control mice are still the best option available to use for the statistical comparison.
# heard back from Lynette (attached).  The stray het can be excluded for this line (CR1760,1xHet and  16xHoms)
library(knitr)
library(PhenStat)
library(DRrequiredAgeing)
library(nlme)
library(plyr)
library(OpenStats)
library(DescTools)
library(reshape2)
library(digest)
library(sjPlot)
library(lattice)
library(reshape2)
library(abind)
dataFramerows <- function(x) {
  xx <- as.data.frame(t(x))
  return(xx)
}
plotGene <- function(td, etc, centre = "UCD") {
  par(mar = c(10, 5, 5, 5))
  allTerms <- OpenStats:::formulaTerms(td$input$fixed)
  omitingTerm <- allTerms[grepl(pattern = "Genotype|Group", x = allTerms)]
  td2 <- OpenStatsAnalysis(
    OpenStatsListObject = td$input$OpenStatsList,
    method = "MM",
    MM_fixed = update.formula(
      old = td$input$fixed,
      new = paste0(
        ". ~ . -",
        paste(
          omitingTerm,
          sep = "-", collapse = "-"
        )
      )
    ),
    MM_random = ~ 1 | id / Batch,
    MM_BodyWeightIncluded = FALSE,
    MM_weight = NULL,
    MM_optimise = c(0, 0, 0, 0, 0, 0),
    correlation = if (centre %in% "UCD") NULL else corSymm(form = ~ 1 | id / Batch)
  )
  td2 <- td2$output$Final.Model
  #########
  x <- getData(td2)
  y <- resid(td2)
  main <- paste(
    "Centre: ",
    unique(x$centre)[1],
    ", Gene symbol:",
    unique(x$Gene_symbol)[1],
    ifelse(length(etc) > 0, paste(etc, collapse = ", "), ""),
    collapse = " "
  )
  boxplot(
    y ~ x$Biological_sample_group * x$Group,
    col = 2:3,
    main = main,
    xlab = "",
    #' Genotype',
    ylab = paste0("Residuals from [", OpenStats:::printformula(formula(td2)), "]"),
    sub = "",
    las = 3,
    cex.main = .7
  )
  abline(
    v = c(2 * 1:2 + .5),
    lty = 4,
    col = "gray"
  )
  # abline(lm(y~x$Genotype),lty = 2,col='gray',lwd=3)
}
############################
plotGeneSex <- function(td, etc, centre = "UCD") {
  par(mar = c(10, 5, 5, 5))
  allTerms <- OpenStats:::formulaTerms(td$input$fixed)
  omitingTerm <- allTerms[grepl(pattern = "Genotype|Group|Sex", x = allTerms)]
  td2 <- OpenStatsAnalysis(
    OpenStatsListObject = td$input$OpenStatsList,
    method = "MM",
    MM_fixed = update.formula(
      old = td$input$fixed,
      new = paste0(
        ". ~ . -",
        paste(
          omitingTerm,
          sep = "-", collapse = "-"
        )
      )
    ),
    MM_random = ~ 1 | id / Batch,
    MM_BodyWeightIncluded = FALSE,
    MM_weight = NULL,
    MM_optimise = c(0, 0, 0, 0, 0, 0),
    correlation = if (centre %in% "UCD") NULL else corSymm(form = ~ 1 | id / Batch)
  )
  td2 <- td2$output$Final.Model
  x <- getData(td2)
  y <- resid(td2)
  main <- paste(
    "Centre: ",
    unique(x$centre)[1],
    ", Gene symbol:",
    unique(x$Gene_symbol)[1],
    ifelse(length(etc) > 0, paste(etc, collapse = ", "), ""),
    collapse = " "
  )
  boxplot(
    y ~ x$Biological_sample_group * x$Sex * x$Group,
    col = 2:5,
    main = main,
    xlab = "",
    #' Genotype',
    ylab = paste0("Residuals from [", OpenStats:::printformula(formula(td2)), "]"),
    sub = "",
    las = 3,
    cex.main = .7
  )
  abline(
    v = c(2 * 1:5 + .5),
    lty = 4,
    col = "gray"
  )
  # abline(lm(y~x$Genotype),lty = 2,col='gray',lwd=3)
}
# library(summarytools)
############################
files <- c(
  "https://www.ebi.ac.uk/~hamedhm/PainWG/Round%2012%20-28-07-2020/Data/csv/JAX_VFR_data_collection_QCedJune.zip",
  "https://www.ebi.ac.uk/~hamedhm/PainWG/Round%2012%20-28-07-2020/Data/csv/HAR_VFR_data_SUDO_converted.zip",
  "https://www.ebi.ac.uk/~hamedhm/PainWG/Round%2012%20-28-07-2020/Data/csv/TCP_VFR_data_collection_sheet_2020-06-15.zip",
  "https://www.ebi.ac.uk/~hamedhm/PainWG/Round%2012%20-28-07-2020/Data/csv/UCD_VFR_data_collection_sheet%20022820%20lrb.zip"
)

for (i in 1:length(files)) {
  df <- read.csv(
    file = DRrequiredAgeing::UnzipAndfilePath(files[i]),
    check.names = FALSE,
    skip = 1,
    stringsAsFactors = TRUE
  )
  df <- OpenStats:::trimColsInDf(df)
  ###### Discussed in the group, HOMS and HEMIS could be combined
  levels(df$Zygosity)[levels(df$Zygosity) %in% "hemizygous"] <- "homozygous"
  # UCD only
  levels(df$Zygosity)[levels(df$Zygosity) %in% "Hemizygous"] <- "homozygous"
  levels(df$Zygosity)[levels(df$Zygosity) %in% "Homozygous"] <- "homozygous"
  levels(df$Zygosity)[levels(df$Zygosity) %in% "Heterozygous"] <- "heterozygous"
  levels(df$Zygosity)[levels(df$Zygosity) %in% "Wild type"] <- "baseline"
  table(df$Zygosity)
  ######
  df <- DRrequiredAgeing:::removeLeadingSpaceFromDataFrameFactors(df)
  df <- df[!apply(df == "", 1, all), ]
  df <- df[!duplicated(df), ]
  # Change manually!
  # JAX  '%m/%d/%Y'
  # HRWL '%d/%m/%Y'
  # TCP  '%Y-%m-%d'
  # UCD  '%d/%m/%Y''
  #### Never change below names!
  centre <- c("JAX", "MRC Harwell", "TCP", "UCD")[i]
  TryTheseFormatsForDates <- c("%m/%d/%Y", "%d/%m/%Y", "%Y-%m-%d", "%d/%m/%Y")[i]
  df$age_in_weeks <- floor((
    as.Date(df$`Date of experiment`, tryFormats = TryTheseFormatsForDates) -
      as.Date(df$`Date of birth`, tryFormats = TryTheseFormatsForDates)) / 7)
  ###############################
  # any negative age!
  if (all(dim(df[df$age_in_weeks < 0, ]) > 0)) {
    df[df$age_in_weeks < 0, c("Date of experiment", "Date of birth", "Is Baseline?")]
  }
  df.org <- df
  df <- df[df$age_in_weeks > 0, ]
  if (!identical(df, df.org)) {
    message("There are some negative ages in the data!")
  }
  ###############################
  addmargins(table(df$age_in_weeks), margin = c(1))
  ###############################
  pie2 <- function(x, main = "", minThreshold = 3, ...) {
    toCondense <- names(which(x < minThreshold))
    if (length(toCondense) > 1) {
      Others <- sum(x[names(x) %in% toCondense])
      x <- c(x[!names(x) %in% toCondense], Other = Others)
    }
    pie(
      x,
      labels = paste("", names(x), " weeks [", unname(round(x / sum(
        x
      ) * 100, 1)), "%]", sep = ""),
      main = main,
      sub = ifelse(length(toCondense) > 1, paste0(
        "Other = ", paste(toCondense, sep = ",", collapse = ","), " weeks"
      ), ""),
      ...
    )
  }
  #############################
  # par(mfrow=c(2,2))
  # pie2(table(df$age_in_weeks[df$`Is Baseline?`  != TRUE]	),main=paste(centre,'- KO'),col=1:5,minThreshold = 10)
  # pie2(table(df$age_in_weeks[df$`Is Baseline?`  == TRUE]	),main=paste(centre,'- WT'),col=6:11)
  # table(df$`Gene symbol`)
  # addmargins(table(df$`Colony name`))
  # table(df$`Colony name`,df$Zygosity)
  # table(df$`Colony name`,df$Strain)
  # addmargins(table(df$`Is Baseline?`))
  # nlevels(droplevels(df$`Colony name`[df$`Is Baseline?`== FALSE & df$`Colony name` !='']))
  # nrow(df)
  # ##############################
  # addmargins(table(df$Sex,df$`Is Baseline?`),margin = c(1,2))
  # addmargins(table(df$Sex,df$age_in_weeks),margin = c(1,2))
  # addmargins(table(df$Sex,df$Zygosity),margin = c(1,2))
  # #dfSummary(df,style = 'grid')
  # ###############################
  # if (sum(df$`Bodyweight (g)`, na.rm = TRUE) > 0) {
  #   # plot(df$age_in_weeks, df$`Bodyweight (g)`)
  # }
  ###############################
  # fix a QC issue in the data
  df <- df[!((df$`Colony name` %in% "CR1760") & (df$Zygosity %in% "heterozygous")), ]

  # Just to make sure that the controls are labeled in data
  table(df$`Colony name`, df$`Is Baseline?`)
  levels(df$`Colony name`)[levels(df$`Colony name`) %in% c("", "C57BL/6N Tac-USA", "109B6NCrlMVP")] <- "Control"
  nlColony <- levels(df$`Colony name`)
  table(df$`Colony name`, df$`Is Baseline?`)
  df$`Is Baseline?`[df$`Colony name` %in% "Control"] <- TRUE
  table(df$`Colony name`, df$`Is Baseline?`)
  ############################
  # table(df$`Colony name`)
  # table(df$`Gene symbol`)
  ############################
  # Fix typos (Roby)
  replaceLevel <- function(x, level, replaceby) {
    if (any(levels(x) %in% level)) {
      message(
        "Levels found in the data levels! ",
        OpenStats:::pasteComma(level),
        " replace by ",
        OpenStats:::pasteComma(replaceby)
      )
    }
    levels(x)[levels(x) %in% level] <- replaceby
    return(x)
  }
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Gapvd2", replaceby = "Gapvd1")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Ppp2r 5c", replaceby = "Ppp2r5c")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Piezo2 ", replaceby = "Piezo2")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Stk36 ", replaceby = "Stk36")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Trpa1 ", replaceby = "Trpa1")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "Emp", replaceby = "Emp1")
  df$`Gene symbol` <- replaceLevel(df$`Gene symbol`, level = "PINK1", replaceby = "Pink1")
  df <- droplevels(df)
  table(df$`Gene symbol`)
  sort(names(table(df$`Gene symbol`)))
  ################
  df$Sex <- replaceLevel(df$Sex, level = "M", replaceby = "Male")
  df$Sex <- replaceLevel(df$Sex, level = "male", replaceby = "Male")
  df$Sex <- replaceLevel(df$Sex, level = "F", replaceby = "Female")
  df$Sex <- replaceLevel(df$Sex, level = "female", replaceby = "Female")
  df <- droplevels(df)
  table(df$Sex)
  table(df$Zygosity)
  ##############################
  length(df$`Colony name`[df$`Is Baseline?` == "TRUE"])
  length(unique(df$`Colony name`[df$`Is Baseline?` != "TRUE"]))
  ################
  # Only Maged1 can be hemizygous or homozygous. The other genes are of one type
  # Can combine, said Jacqu (Roby)
  df <- df[!(df$Zygosity %in% "homozygous" & df$`Gene symbol` == "Pah"), ]
  df <- droplevels(df)
  t(table(df$Zygosity, df$`Gene symbol`))
  t(table(df$Zygosity, df$`Colony name`))
  ############################
  ############################
  # Detect column of interest
  ############################
  cols <- rep(FALSE, length(names(df)))
  # Included
  for (pInc in c(
    # UCD
    # JAX
    "JAX_VFR_003_001",
    "JAX_VFR_007_001",
    "JAX_VFR_011_001",
    # TCP
    "TCP_VFR_005_001",
    "TCP_VFR_010_001",
    "TCP_VFR_015_001",
    # HRWL
    "HRWL_VFR_005_001",
    "HRWL_VFR_011_001",
    "HRWL_VFR_017_001",
    # UCD
    "UCD_VFR_004_001-INCREMENT:1",
    "UCD_VFR_009_001-INCREMENT:1",
    "UCD_VFR_014_001-INCREMENT:1"
  )) {
    cols <- cols | grepl(
      x = names(df),
      pattern = pInc,
      fixed = TRUE
    )
  }
  # pairs(df[,cols],panel = function(x,y){points(x,y);abline(lm(y~x),lwd=2,col=2,lty=2)})
  for (pExc in c("Eception-not-required-here")) {
    cols <- cols & !grepl(
      x = names(df),
      pattern = pExc,
      fixed = TRUE
    )
  }
  cols <- names(df)[cols]
  message("Response variables: ", paste(cols, collapse = ", "))
  ############################
  # Area Under the Curve (AUC) calculation (for Hamed)
  AUC <- apply(df[, head(cols, 3)], 1, function(x) {
    AUC(x = c(0, 24, 48), y = x[1:3])
  })
  ############################
  # Hash metadata to make it easier to work with
  metadata_colsIndex <- which(grepl(
    x = names(df),
    pattern = "(HRWL_VFR_019_001)|(JAX_VFR_013_001)|(TCP_VFR_016_001)|(UCD_VFR_016_001)",
    fixed = FALSE
  )):ncol(df)
  exclude <- which(grepl(pattern = "(JAX_VFR_039_001)|(age_in_weeks)|(HRWL_VFR_041_001)|(TCP_VFR_043_001)|(TCP_VFR_033_001)|(UCD_VFR_042_001)", x = names(df)))
  metadata_colsIndex <- metadata_colsIndex[!metadata_colsIndex %in% exclude]
  Metadata_group <- apply(df[, metadata_colsIndex], 1, function(x) {
    digest(paste0(trimws(x), collapse = ""))
  })
  # write(unique(Metadata_group),file = 'd:/metadata.txt')
  table(Metadata_group)
  ###########################
  # Creating the dataset for the analysis
  ###########################
  # Because the age is at least 100 days, subtract age by 100 to avoid scaling issues
  ###########################

  df_initial <- data.frame(
    id                      = df$`Animal name`,
    Colony_name             = df$`Colony name`,
    Biological_sample_group = df$`Is Baseline?`,
    Sex                     = df$Sex,
    Zygosity                = df$Zygosity,
    Strain                  = df$Strain,
    Weight                  = df$`Bodyweight (g)`, ### !
    Gene_symbol             = df$`Gene symbol`,
    Batch                   = df$`Date of experiment`,
    Metadata_group          = Metadata_group,
    age                     = as.integer(df$age_in_weeks),
    centre                  = centre,
    AUC                     = AUC,
    df[, cols],
    check.names = FALSE
  )
  ###############################
  # dfSummary(df_initial)
  df_initial$Biological_sample_group <- as.factor(ifelse(as.logical(df_initial$Biological_sample_group), "control", "mutant"))
  df_initial <- droplevels(df_initial)
  ###############################
  ##############################
  # summary(df_initial)
  # par(mfrow=c(3,3))
  # lapply(head(cols, 3), function(x) {
  #   message('Parameter: ',x)
  #   density0= function(x,...){
  #     r = density(x = x,bw=.15,...)
  #     return(r)
  #   }
  #   plot(density0(log(df_initial[, x]), na.rm = TRUE),
  #        main = paste0('Log (', x, ')'),
  #        cex.main = .85)
  #   plot(
  #     density0(log(df_initial[df_initial$Sex %in% 'Male', x]), na.rm = TRUE),
  #     main = paste0('Male - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 2
  #   )
  #   plot(
  #     density0(log(df_initial[df_initial$Sex %in% 'Female', x]), na.rm = TRUE),
  #     main = paste0('Female - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 3
  #   )
  #   ################################
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'control', x]), na.rm = TRUE),
  #     main = paste0('Control - Log (', x, ')'),
  #     cex.main = .85
  #   )
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'control' &
  #                              df_initial$Sex %in% 'Male', x]), na.rm = TRUE),
  #     main = paste0('Control*Male - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 2
  #   )
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'control' &
  #                              df_initial$Sex %in% 'Female', x]), na.rm = TRUE),
  #     main = paste0('Control*Female - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 3
  #   )
  #   ################################
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'mutant', x]), na.rm = TRUE),
  #     main = paste0(centre,': Mutant - Log (', x, ')'),
  #     cex.main = .85
  #   )
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'mutant' &
  #                              df_initial$Sex %in% 'Male', x]), na.rm = TRUE),
  #     main = paste0(centre,': Mutant*Male - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 2
  #   )
  #   plot(
  #     density0(log(df_initial[df_initial$Biological_sample_group %in% 'mutant' &
  #                              df_initial$Sex %in% 'Female', x]), na.rm = TRUE),
  #     main = paste0(centre,': Mutant*Female - Log (', x, ')'),
  #     cex.main = .85,
  #     col = 3
  #   )
  # })
  ##############################
  ##############################
  # Melt dataset to make the levels right for the repeated measures
  df_melt <- reshape2::melt(df_initial,
    measure.vars = head(cols, 3),
    variable.name = "Group"
  )
  levels(df_melt$Group) <- c("Test0", "Test1", "Test2")
  df_melt <- droplevels(df_melt)
  # Discussed with Sonia to remove Test0 as it is inhabituation and not challenged in MRC Harwell
  if (centre %in% "MRC Harwell") {
    # Jacquie and Janine email
    # https://mail.google.com/mail/u/0/#search/janine/FMfcgxwJWrVzHFXDrKbzkMvLgxcnhVrQ
    # df_melt <- droplevels(subset(df_melt, !df_melt$Group %in% "Test0"))
  }
  ###############################
  # Separate controls from mutants
  df2_control <- droplevels(subset(
    df_melt,
    df_melt$Biological_sample_group %in% "control"
  ))
  df2_mutant <- droplevels(subset(
    df_melt,
    df_melt$Biological_sample_group %in% "mutant"
  ))
  table(df2_control$Strain)
  table(df2_mutant$Strain)
  table(df2_mutant$Group)
  ###############################
  # Just for HAMED
  # par(mfrow = c(1, 1), mar = c(15, 5, 3, 3))
  # boxplot(
  #   log(df2_control$value) ~
  #     df2_control$Group*df2_control$Sex  ,
  #   col = 2:4,
  #   las = 3,
  #   main = 'Baseline',
  #   ylab = 'Log response',
  #   horizontal = FALSE,
  #   xlab = ''
  # )
  # boxplot(
  #   log(df2_mutant$value) ~
  #     df2_mutant$Group * df2_mutant$Sex + 1,
  #   col = 4:6,
  #   las = 3,
  #   main = paste0(centre,': Mutant'),
  #   horizontal = FALSE,
  #   ylab = 'Log Response',
  #   xlab = ''
  # )
  # boxplot(log(df_melt$value) ~
  #           df_melt$Biological_sample_group*df_melt$Group*df_melt$Sex + 1,
  #         col = 2:5,
  #         las = 3,
  #         main = 'Dataset',
  #         horizontal = FALSE,
  #         ylab = 'Log Response',
  #         xlab = '',
  #         srt=45
  # )
  # dfSummary(df2_control)
  # dfSummary(df2_mutant)
  summary(df_melt$value ~ df_melt$Sex + df_melt$Group)
  ############################
  ###############################
  ###############################
  message("Checking for the age variation in the data ...")
  Onlycontrols <- data.frame(df2_control, FakeId = rbinom(
    n = nrow(df2_control),
    size = 1,
    prob = .5
  ))
  Onlycontrols$age <- as.factor(Onlycontrols$age)
  plAgeVariation <- OpenStatsList(
    dataset = Onlycontrols,
    dataset.colname.batch = "Batch",
    dataset.colname.genotype = "FakeId",
    dataset.colname.sex = "Sex",
    dataset.colname.weight = "Weight",
    refGenotype = unique(Onlycontrols$FakeId)[1],
    testGenotype = unique(Onlycontrols$FakeId)[2]
  )
  table(
    plAgeVariation@datasetPL$Biological_sample_group,
    plAgeVariation@datasetPL$Sex
  ) / 3 # 3 times an animal
  ###############################
  plAgeVariation@datasetPL$logValue <- log(plAgeVariation@datasetPL$value)
  ###############################
  # Here is the analysis engine
  # Model: Optimised Linear mixed model (MM) using nlme package in R (optimisation using AICc)
  tdAgeVariation <- OpenStatsAnalysis(
    OpenStatsListObject = plAgeVariation,
    method = "MM",
    MM_fixed = logValue ~ Sex * Group,
    MM_random = ~ 1 | id / Group / Batch,
    correlation = corSymm(form = ~ 1 | id / Group / Batch),
    MM_weight = NULL,
    MM_lower = ~1,
    MM_optimise = c(1, 1, 1, 0, 0, 0)
  )
  #  pie2(x = table(plAgeVariation@datasetPL$age))
  #  plot(tdAgeVariation)
  gd <- getData(tdAgeVariation$output$Final.Model)
  gd$resids <- resid(tdAgeVariation$output$Final.Model)
  tsdGroup <- TukeyHSD(aov(
    reformulate(termlabels = "age", response = "resids"),
    data = gd
  ))
  plot(tsdGroup, sub = centre)
  tsdGroupTable <- dataFramerows(tsdGroup$age)[4, , drop = FALSE]
  table(gd$age)
  write.csv(
    x = data.frame(
      tsdGroupTable,
      check.names = FALSE
    ),
    file = paste0("Age_comparision_", centre, ".csv"),
    row.names = FALSE
  )
  ###############################
  # The order of the loops is: Filter on Zygosity, Strain, Metadata and colonies.
  Rlist <- NULL
  counter <- 1
  ###############################
  for (zyg in unique(df2_mutant$Zygosity)) {
    for (stra in unique(df2_mutant$Strain)) {
      for (meta in unique(df2_mutant$Metadata_group)) {
        for (colony in unique(df2_mutant$Colony_name)) {
          for (cstrain in unique(df2_control$Strain)) {
            if(centre %in% 'JAX' && cstrain %in% 'C57BL/6J')
              next
            ###############################
            # Select the mutants###########
            df_mut_filtered <- droplevels(
              subset(
                df2_mutant,
                df2_mutant$Colony_name %in% colony &
                  df2_mutant$Zygosity %in% zyg &
                  df2_mutant$Strain %in% stra
              )
            )
            message("##########################################")
            message(paste(zyg, stra, meta, colony, sep = " ~> "))
            message("##########################################")
            ###############################
            # This is the final dataset that would be fed into the model
            df_both <- rbind(
              df_mut_filtered,
              subset(df2_control, df2_control$Strain %in% cstrain)
            )

            # for (ageGroup in unique(df_both$age)) {
            # adf_both = subset(df_both, df_both$age %in% ageGroup)
            df_both <- droplevels(df_both)
            ageGroup <- paste(unique(df_both$age), sep = "~", collapse = "~")
            ###############################
            # Just to skip errors
            if (nlevels(df_both$Colony_name) < 2) {
              message("Only on level in genotype .... ")
              write(
                paste(zyg, stra, meta, colony, ageGroup, sep = ","),
                file = "SkipFromAnalysisNoControlOrMutant.csv",
                append = TRUE,
                ncolumns = 10^3
              )
              next
            }
            ################################
            tsdGroup <- TukeyHSD(aov(
              reformulate(termlabels = c("Group*Sex*Biological_sample_group"), response = "value"),
              data = df_both
            ))
            ff <- function(x, col = 4) {
              f1 <- unlist(lapply(x, function(x) {
                x[, col]
              }))
              f2 <- data.frame(t(f1))
              names(f2) <- names(f1)
              return(f2)
            }
            dfExperimentDetails <- data.frame(
              colony_id = colony,
              control_strain = cstrain,
              mutant_strain = stra,
              zygosity = zyg,
              metadata = meta
            )
            if (counter <= 1) {
              comparisons <- cbind(dfExperimentDetails, ff(tsdGroup, 4))
              comDifferences <- cbind(dfExperimentDetails, ff(tsdGroup, 1))
            } else {
              comparisons <- rbind.fill(comparisons, cbind(dfExperimentDetails, ff(tsdGroup, 4)))
              comDifferences <- rbind.fill(comDifferences, cbind(dfExperimentDetails, ff(tsdGroup, 1)))
            }
            ################################
            ###############################
            # plot(
            #   y = tapply(
            #     df_both[, c('value')],
            #     INDEX = df_both$id,
            #     FUN = function(x) {
            #       mean(x)
            #     }
            #   ),
            #   x = tapply(
            #     df_both[, c('age')],
            #     INDEX = df_both$id,
            #     FUN = function(x) {
            #       mean(x)
            #     }
            #   ),
            #   xlab='age mean',
            #   ylab = 'Response'
            # )

            # plot(
            #   as.integer(as.Date(df_both$Batch, tryFormats = TryTheseFormatsForDates)),
            #   df_both$value,
            #   xlab =
            #     'days from the Origin 1900',
            #   col = abs(as.integer(df_both$Biological_sample_group) + 2),
            #   pch = as.integer(df_both$Biological_sample_group),
            #   ylab='Response'
            # )
            # abline(h=mean(df_both$value[df_both$Biological_sample_group %in% 'control']),col=1,lty=3,lwd=3)
            # abline(h=mean(df_both$value[df_both$Biological_sample_group %in% 'mutant']),col=2,lty=2,lwd=3)
            # ##############################
            # We process the data using the new IMPC stats pipeline.
            ###############################
            # Creating a PhenList object. Do not know what is PhenList? Check PhenStat package (?PhenList)
            pl <- OpenStatsList(
              dataset = df_both,
              dataset.colname.batch = "Batch",
              dataset.colname.genotype = "Colony_name",
              dataset.colname.sex = "Sex",
              dataset.colname.weight = "Weight",
              refGenotype = unique(df2_control$Colony_name),
              testGenotype = colony
            )
            table(pl@datasetPL$Biological_sample_group, pl@datasetPL$Sex) / 3 # 3 times an animal
            ###############################
            pl@datasetPL$logValue <- log(pl@datasetPL$value)
            ###############################
            # Here is the analysis engine
            # Model: Optimised Linear mixed model (MM) using nlme package in R (optimisation using AICc)
            td <- OpenStatsAnalysis(
              OpenStatsListObject = pl,
              method = "MM",
              MM_fixed = logValue ~ Genotype * Sex * Group ,#+ age,
              MM_random = ~ 1 | id / Group / Batch,
              correlation = corSymm(form = ~ 1 | id / Group / Batch),
              MM_weight = varIdent(~ 1 | Genotype),
              MM_optimise = c(1, 1, 1, 1, 1, 1)
            )
            ###############################
            # Just checking the results ...
            if (is.null(td$messages)) {
              message("plotting in progress ...")
              plotGene(td, etc = c(
                paste0("\nC B.strain = ", cstrain),
                paste0("M B.strain = ", stra),
                paste0("Colony = ", colony),
                paste0("Zyg = ", zyg),
                paste0("\nage = ", ageGroup, "w")
              ))
              plotGeneSex(td, etc = c(
                paste0("\nC B.strain = ", cstrain),
                paste0("M B.strain = ", stra),
                paste0("Colony = ", colony),
                paste0("Zyg = ", zyg),
                paste0("age = ", ageGroup, "w")
              ))
              # plot(td,
              #      main = paste0(unique(df_both$Gene_symbol)[1]),
              #      mfrow = c(2, 2)
              # )
              # Summary from PhenStatAgeing           ####################
              summary(td)
              # Normal Summary/anova from lme package ####################
              summary(td$output$Final.Model)
              # anova(td$output$Final.Model)
              report <- OpenStatsReport(
                td,
                JSON = TRUE,
                RemoveNullKeys = FALSE,
                ReportNullSchema = TRUE,
              )
              getValue <- function(x, col = "Gene_symbol", sep = "~", splitGenotype = FALSE) {
                if ("Biological_sample_group" %in% names(x)) {
                  rM <- paste(sort(unique(x[, col][x$Biological_sample_group %in% "mutant"])), sep = sep, collapse = sep)
                  rC <- paste(sort(unique(x[, col][!x$Biological_sample_group %in% "mutant"])), sep = sep, collapse = sep)
                  if (splitGenotype) {
                    return(paste0("Controls:", rC, " | Mutants: ", rM))
                  } else {
                    return(rM)
                  }
                } else {
                  return(NA)
                }
              }
              tsdGroupR <- TukeyHSD(aov(
                reformulate(termlabels = "Group", response = "logValue"),
                data = pl@datasetPL # dcast(data = pl@datasetPL,formula = id~Group)
              ))
              # Because MRC Harwell does not have Test0 included in the analysis
              if (nrow(tsdGroupR$Group) < 3) {
                # For Harwell temporary suppresed
                # https://mail.google.com/mail/u/0/#search/janine/FMfcgxwJWrVzHFXDrKbzkMvLgxcnhVrQ
                #tsdGroupR$Group <- rbind(NA, NA, tsdGroupR$Group)
                #rownames(tsdGroupR$Group) <- c("Test1-Test0", "Test2-Test0", "Test2-Test1")
              }

              tsdGroupTableR <- dataFramerows(tsdGroupR$Group)[4, , drop = FALSE]

              write(
                names(tsdGroupTableR),
                file = "just_order_of_the_testAnova.log",
                append = TRUE,
                ncolumns = 100
              )
              report <- c(
                centre,
                getValue(df_both, col = "Gene_symbol", splitGenotype = FALSE),
                colony,
                zyg,
                stra,
                cstrain,
                meta,
                getValue(x = df_both, col = "age", splitGenotype = TRUE),
                tsdGroupTableR[1],
                tsdGroupTableR[2],
                tsdGroupTableR[3],
                report
              )
              write(
                paste(report, collapse = "\t"),
                file =
                  DRrequiredAgeing:::RemoveSpecialChars(
                    paste(round(runif(1, 1, max = 10^5)), centre, colony, zyg, stra, meta, ageGroup, ".tsv", sep = "-"),
                    what = "[^0-9A-Za-z.]"
                  ),
                ncolumns = 10^5,
                append = TRUE
              )


              # just for Hamed ##############################
              # xyplot(
              # 	value ~ Group,
              # 	groups = Biological_sample_group,
              # 	df_both,
              # 	panel = function(x, y, subscripts, groups, ...) {
              # 		yhat = fitted(lm(y ~ groups * x))
              # 		#panel.superpose(x, y, subscripts, groups, pch = 1:10,type = 'b')
              # 		panel.superpose(x, yhat, subscripts, groups, type = "l")
              # 	},
              # 	xlab = "Group",
              # 	ylab = 'value'
              # )
              # plot_model(td$output$Final.Model)
              # plot_model(td$output$Final.Model, type = 'eff')
              # plot_model(td$output$Final.Model, type = 'resid')
              # interaction.plot(
              # 	td$input$data$Group,
              # 	x.factor = td$input$data$id,
              # 	response = td$input$data$value,
              # 	col = td$input$data$id,
              # 	lty = td$input$data$Biological_sample_group + 1,
              # 	lwd = abs(td$input$data$Biological_sample_group - 2) * 2
              # )
              ################
              Rlist[[counter]] <- td
              counter <- counter + 1
            }
            # }
          }
        }
      }
    }
  }
  save(Rlist, file = paste0(centre, "-AllResults.Rdata"))
  rna = paste0(sample(x = letters,size = 3),collapse = '')
  write.csv(
    comparisons,
    file =
      DRrequiredAgeing:::RemoveSpecialChars(paste(rna,"_MutualComparisonPvalue", centre, ".csv", sep = "-"),
        what = "[^0-9A-Za-z.]"
      ),
    row.names = FALSE
  )
  write.csv(
    comDifferences,
    file =
      DRrequiredAgeing:::RemoveSpecialChars(paste(rna,"_MutualComparisonDifferences", centre, ".csv", sep = "-"),
        what = "[^0-9A-Za-z.]"
      ),
    row.names = FALSE
  )
}
graphics.off()
###############################
# The object containing all results
# Rlist
###############################
# suspicious animals
# suspData = lapply(Rlist, function(x) {
# 	getData(x$output$Final.Model)[resid(x$output$Final.Model) > .04,]
# })
# susCases = abind(suspData, along = 1)
# write.csv(a[!duplicated(susCases), ], file = 'd:/Suspicious.csv',row.names = FALSE)
###############################
