\name{PhenStatReport}
\alias{PhenStatReport}
\title{
  Generating a pdf report from a PhenList object
}
\description{
  This function takes a Phenlist object and generates a pdf report containing several statistical methods.
}
\usage{
  PhenStatReport(
    PhenlistObject           ,
    depVariable = NULL       ,
    other.response = NULL    ,
    update = TRUE            ,
    Gene.Symbol = NULL       ,
    Response.name = NULL     ,
    destination = NULL       ,
    reportTitle = "Extended Statistical Report",
    DataRelease = NULL       ,
    Showsource = FALSE       ,
    open = FALSE             ,
    clean = TRUE             ,
    verbos = FALSE           ,
    ...
  )
}
\arguments{
  \item{PhenlistObject}{
    A phenlist object that is already created by PhenList() function in PhenStat.
  }
  \item{depVariable}{
    String. Name of the dependent variable.
  }
  \item{other.response}{
    The vector of strings. A vector of names containing the other dependent variables in the data set. Default NULL
  }
  \item{update}{
    Logical flag. Set to TRUE to get the latest version of the report on the fly. Default is TRUE
  }
  \item{Gene.Symbol}{
    Optional string. Gene symbol. Default NULL
  }
  \item{Response.name}{
    Optional string. Name of the dependent variable. Default NULL
  }
  \item{destination}{
    Location of the final report file. The default is the working directory of R.
  }
  \item{reportTitle}{
    The title of the report that is printed on the top of the first page of the report.
  }
  \item{DataRelease}{
    Optional flag. Data release version. Default NULL
  }
  \item{Showsource}{
    Logical flag. Set to TRUE to see the actual R codes that generate each section of the report. Default FALSE
  }
  \item{open}{
    Logical flag. Set to TRUE to open the report after it is generated. Default FALSE
  }
  \item{clean}{
    Logical flag. Set to TRUE to remove the auxiliary files after successfully generating the report. Default TRUE
  }
  \item{verbos}{
    Logical flag. Setting to TRUE shows the details and progress of the report generating function on screen. Default FALSE
  }
  \item{\dots}{
    List of other parameters that can be passed to the functions.
  }
}
\details{
  See \href{https://www.mousephenotype.org/phenstat-report}{PhenStatReport page on the IMPC website}
}
\value{
  \item{PhenlistObject}{The input object}
  \item{depVariable}{Dependent variable}
  \item{texfile}{The location of the .tex file, given Clean=FALSE}
  \item{pdffile}{The location of the final pdf file}
  \item{...}{The same as the input parameters}
}
\author{
  Hamed Haselimashhadi <hamedhm@ebi.ac.uk>
}



  \seealso{
    \code{\link{PhenList}}
  }
\examples{
    file <- system.file("extdata", "test1.csv", package = "PhenStat")
    test = PhenStat:::PhenList(dataset = read.csv(file,na.strings = '-'),
                               testGenotype = "Sparc/Sparc")
 \dontrun{
    PhenStatReport(test,
                   depVariable = 'Bone.Area',
                   open = TRUE)
  }
}


