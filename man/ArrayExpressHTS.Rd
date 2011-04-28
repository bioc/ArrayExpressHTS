\name{ArrayExpressHTS}
\alias{ArrayExpressHTS}
\title{ ExpressionSet for RNA-Seq experiment submitted in ArrayExpress and ENA }
\description{
    \code{ArrayExpressHTS} runs the RNA-Seq pipeline on a transcription profiling experiment 
    available on the ArrayExpress database and produces an 
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} R object. 
    \code{ArrayExpressHTS} requires an Internet connection.
}
\usage{
    ArrayExpressHTS(accession, usercloud = TRUE, 
        options = getAEDefaultOptions(), nnodes = 10, 
        pool = "32G", attempts = 4, dir = getwd(), 
        refdir = getDefaultReferenceDir(), filter = TRUE, 
        want.reports = TRUE )
}
\arguments{
    \item{accession}{ an ArrayExpress experiment accession identifier. }
    \item{usercloud}{ if TRUE, the R-Cloud will be used to schedule and process the experiment 
        in parallel, otherwise data files are processed sequentially. }
    \item{options}{ default options .}
    \item{nnodes}{ if set, the selected amount of nodes will be allocated when 
        R-Cloud cluster is created. Not used when \code{usercloud} is set to FALSE. }
    \item{pool}{ server pool, from which cluster nodes are allocated. Allowed values are: 
        'default', '4G', '8G', '16G', '32G'. Not used when \code{usercloud} is set to FALSE. }
    \item{attempts}{ number of attempts the package uses to allocate server node before giving up.
         Not used when \code{usercloud} is set to FALSE. }
    \item{dir}{ folder where experiment data will be stored and processed. Default is current directory. }
    \item{refdir}{ the directory where reference data is located. }
    \item{filter}{ if TRUE, data filtering will be used as part of the the pipeline. }
    \item{want.reports}{ if TRUE, quality reports are produced, however, it usually takes longer and more memory is used. For faster computation, set to FALSE. }
}
\value{
    The output is an object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing expression values in 
    assayData (corresponding to the raw sequencing data files), the information contained in 
    the .sdrf file in phenoData, the information in the adf file in featureData and the idf 
    file content in experimentData.
    
    If executed on a local PC, make sure that tools are available to the pipeline. Check 
    \code{\link{prepareAnnotation}} to see what needs to be done to make tools available.
    
}
    \seealso{\code{\link{ArrayExpressHTSFastQ}}, \code{\link{prepareReference}}, \code{\link{prepareAnnotation}}, \code{\link{prepareAnnotation}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{

if (isRCloud()) { # disabled on local configs so as not to kill package building process

    # if executed on a local PC, make sure tools
    # are available to the pipeline. Check initDefaultEnvironment 
    # help page for instructions.
    expfolder = tempdir();
    
    # run the pipeline
    #
    aehts = ArrayExpressHTS("E-GEOD-16190", dir = expfolder);
    
    # load the expression set object
    loadednames = load(paste(expfolder, "/E-GEOD-16190/eset_notstd_rpkm.RData", sep=""));
    loadednames;

    get('library')(Biobase);

    # print out the expression values
    #
    head(assayData(eset)$exprs);

    # print out the experiment meta data 
    experimentData(eset);
    pData(eset);
}

}


