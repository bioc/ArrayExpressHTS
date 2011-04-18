\name{ArrayExpressHTSFastQ}
\alias{ArrayExpressHTSFastQ}
\title{ ExpressionSet for RNA-Seq raw data files }
\description{
    \code{ArrayExpressHTSFastQ} runs the RNA-Seq pipeline on raw RNA-Seq data 
    files and an .sdrf experiment descriptor and produces an 
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
    R object.
}
\usage{
    ArrayExpressHTSFastQ(accession, organism, 
        usercloud = TRUE, options=getAEDefaultOptions(), 
        nnodes = 10, pool = "32G", attempts = 4, dir = getwd(), 
        refdir = getDefaultReferenceDir(), filter = TRUE, 
        want.reports = TRUE)
}
\arguments{
    \item{accession}{ name of the folder where experiment data is stored. Experiment description files .sdrf and .idf should be stored in this folder. 
        The actual data files should be stored in the 'data' that should be created in this folder. }
    \item{organism}{ this parameter is used to select appropriate alignment reference. }
    \item{usercloud}{ if TRUE, the R-Cloud will be used to schedule and compute the experiment 
    in parallel, otherwise the data files are computed sequentially. }
    \item{options}{ default options .}
    \item{nnodes}{ if set, the selected amount of nodes will be allocated when the
        R-Cloud cluster is created. Not used when \code{usercloud} is set to FALSE. }
    \item{pool}{ server pool, from which cluster nodes are allocated. Allowed values are: 
        'default', '4G', '8G', '16G', '32G'. Not used when \code{usercloud} is set to FALSE. }
    \item{attempts}{ number of attempts the package uses to allocate server node before giving up.
        Not used when \code{usercloud} is set to FALSE. }
    \item{dir}{ folder where experiment data will be stored and processed. Default is current directory. }
    \item{refdir}{  directory where the reference data is located. }
    \item{filter}{ if TRUE, data filtering will be used as part of the the pipeline. }
    \item{want.reports}{ if TRUE, quality reports are produced, however, it usually takes longer and more memory is used. For faster computation, set to FALSE. }
}
\value{
     The output is an object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing expression
    values in assayData (corresponding to the raw sequencing data files), 
    the information contained in the .sdrf file in phenoData, the information in the adf file in 
    featureData and the idf file content in experimentData.
    
}
    \seealso{\code{\link{ArrayExpressHTS}}, \code{\link{prepareReference}}, \code{\link{prepareAnnotation}}}
\author{
    Andrew Tikhonov
    Angela Goncalves
    Maintainer: <andrew@ebi.ac.uk>  
    Maintainer: <angela.goncalves@ebi.ac.uk>  
}
\examples{

if (isRCloud()) { # disabled on local configs so as not to kill package building process
    
    # In ArrayExpressHTS/expdata there is testExperiment, which is
    # a very short version of E-GEOD-16190 experiment, placed there 
    # for testing reasons.
    #
    # Experiment in ArrayExpress:
    # http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-16190
    #
    # the following piece of code will take ~1.5 hours to compute
    # on local PC and ~10 minutes on R-Cloud
    #
    
    # if executed on a local PC, make sure tools are available
    # to the pipeline. Check initDefaultEnvironment help page 
    # for instructions.
    # 
    
    # create a temporary folder where experiment will be copied
    # computing experiment in the package folder may cause issues
    # with file permissions and therefore failures.
    # 
    #
    srcfolder <- system.file("expdata", "testExperiment", package="ArrayExpressHTS");
    
    dstfolder <- tempdir();
    
    file.copy(srcfolder, dstfolder, recursive = TRUE);
    
    # run the pipeline
    #
    # set usercloud = FALSE if executing on local PC, 
    # therefore parallel computation will be disabled
    #
    aehts = ArrayExpressHTSFastQ(accession = "testExperiment", organism = "Homo_sappiens", dir = dstfolder);
        
    # load the expression set object
    loadednames = load(paste(dstfolder, "/testExperiment/eset_notstd_rpkm.RData", sep=""));
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

