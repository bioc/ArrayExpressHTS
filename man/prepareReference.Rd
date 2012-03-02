\name{prepareReference}
\alias{prepareReference}
\title{ Prepare reference data for the RNA-Seq Pipeline }
\description{
    \code{prepareReference} downloads reference genome or transcriptome for 
    the selected organism from the Ensembl database and processes it so that 
    it can be used by the pipeline.
    \code{prepareReference} requires an Internet connection.
}
\usage{
    prepareReference(organism, version = "current", 
        type = c("genome", "transcriptome"), 
        location = getDefaultReferenceDir(), 
        aligner = c("bwa", "bowtie", "tophat"), 
        refresh = FALSE, run = TRUE)
}
\arguments{
    \item{organism}{ supported organism names can be viewed in the Ensemble database. Check 'ftp://ftp.ensembl.org/pub'. }
    \item{version}{ \code{"current"} or other appropriate version. Check 'ftp://ftp.ensembl.org/pub'. }
    \item{type}{ two values are supported: \code{"genome"}, \code{"transcriptome"} }
    \item{location}{ indicates where the reference data should be stored. }
    \item{aligner}{ 3 types of aligners are supported: \code{"bwa"}, \code{"bowtie"} and \code{"tophat"}. }
    \item{refresh}{ if TRUE, existing reference data will be rebuilt. }
    \item{run}{ if FALSE, the downloading and processing commands will not be executed. }
}
\value{
    The output is the version of the organism reference that has been downloaded and processed.
    The reference files are kept in the folder defined in \code{location} parameter.
}
\seealso{
    \code{\link{ArrayExpressHTS}}, 
    \code{\link{ArrayExpressHTSFastQ}}, 
    \code{\link{prepareAnnotation}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{

if (isRCloud()) {
        
    par(ask = FALSE)

    # the following piece of code will take ~3 hours to complete
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
    
    # download and prepare reference
    prepareReference("Homo_sapiens", version = "GRCh37.61", 
        type = "genome", aligner = "bowtie", location = referencefolder )
    prepareReference("Homo_sapiens", version = "GRCh37.61", 
        type = "transcriptome", aligner = "bowtie", location = referencefolder )
    prepareReference("Mus_musculus", version = "current", 
        type = "genome", aligner = "bowtie", location = referencefolder )
    prepareReference("Mus_musculus", version = "current", 
        type = "transcriptome", aligner = "bowtie", location = referencefolder )
}

}

