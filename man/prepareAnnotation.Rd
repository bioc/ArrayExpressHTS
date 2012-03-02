\name{prepareAnnotation}
\alias{prepareAnnotation}
\title{ Prepare annotation data for the RNA-Seq Pipeline }
\description{
    \code{prepareAnnotation} downloads the required annotation file for the selected organism 
    from Ensembl and processes it so that it can be used by the pipeline.
    \code{prepareAnnotation} requires an Internet connection.
}
\usage{
    prepareAnnotation(organism, version = "current", 
        location = getDefaultReferenceDir(), refresh = FALSE, run = TRUE)
}
\arguments{
    \item{organism}{ supported organism names can be viewed in the Ensemble database. Check 'ftp://ftp.ensembl.org/pub'. }
    \item{version}{ \code{"current"} or other appropriate version. Check 'ftp://ftp.ensembl.org/pub'. }
    \item{location}{ indicates where the annotation data should be stored. }
    \item{refresh}{ if TRUE, existing annotation data will be rebuilt. }
    \item{run}{ if FALSE, the commands to obtain and process the annotation will not be executed. }
}
\value{
    The output is the version of the organism annotation that has been downloaded and processed.
    The annotation files are kept in the folder defined in \code{location} parameter.
}
\seealso{
    \code{\link{ArrayExpressHTS}}, 
    \code{\link{ArrayExpressHTSFastQ}}, 
    \code{\link{prepareReference}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    
if (isRCloud()) { # disabled on local configs 
    # so as not to affect package building process
    
    par(ask = FALSE)

    # the following piece of code will take ~1.5 hours to complete
    #
    
    # if executed on a local PC, make sure tools are available
    # to the pipeline.
    # 
    
    # create directory
    #
    # Please note, tempdir() is used for automamtic test 
    # execution. Select directory more appropriate and 
    # suitable for keeping reference data.
    #
    referencefolder = paste(tempdir(), "/reference", sep = "")
    
    dir.create(referencefolder)
    
    # download and prepare annotation
    prepareAnnotation("Homo_sapiens", "current", location = referencefolder)
    prepareAnnotation("Mus_musculus", "NCBIM37.61", location = referencefolder)

}

}

