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
    ArrayExpressHTS( accession,
            options = list (
                    stranded           = FALSE,
                    insize             = NULL,
                    insizedev          = NULL,
                    reference          = "genome",
                    aligner            = "tophat",
                    aligner_options    = NULL,
                    count_feature      = "transcript",
                    count_options      = "",
                    count_method       = "cufflinks",
                    filter             = TRUE,
                    filtering_options  = NULL ),
            usercloud = TRUE,
            rcloudoptions = list (
                    nnodes     = "automatic",
                    pool       = c("4G", "8G", "16G", "32G", "64G"),
                    nretries   = 4 ),
            steplist = c("align", "count", "eset"),
            dir = getwd(),
            refdir = getDefaultReferenceDir(),
            want.reports = TRUE, 
            stop.on.warnings = FALSE )
    
}
\arguments{
    \item{accession}{ an ArrayExpress experiment accession identifier, e.g. "E-GEOD-16190" }
    \item{options}{ defines pipeline options. See \code{\link{getDefaultProcessingOptions}}. }
    \item{usercloud}{ defines if the R Cloud will be used to parallel experiment computation, 
        if \code{FALSE}, experiment data files will be processed sequentially. }
    \item{rcloudoptions}{ defines R Cloud options. See \code{\link{getDefaultRCloudOptions}}. }
    \item{steplist}{ defines the steps the pipeline will perform. }
    \item{dir}{ folder where experiment data will be stored and processed. Default is current directory. }
    \item{refdir}{ the directory where reference data is located. }
    \item{want.reports}{ defines if quality reports are produced. Reports usually make 
        computation longer and eat up more memory. For faster computation use \code{FALSE}. }
    \item{stop.on.warnings}{ self explanatory. Warnings are normally producesd when there 
        are inconsistencies, which however would allow the result to be produced. }
}
\value{
    The output is an object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing expression values in 
    \code{assayData} (corresponding to the raw sequencing data files), the information contained in 
    the .sdrf file in \code{phenoData}, the information in the adf file in featureData and the idf 
    file content in experimentData.
    
    If executed on a local PC, make sure that tools are available to the pipeline. Check 
    \code{\link{prepareAnnotation}} to see what needs to be done to make tools available.
    
}
\seealso{
    \code{\link{ArrayExpressHTSFastQ}},
    \code{\link{prepareReference}},
    \code{\link{prepareAnnotation}},
    \code{\link{prepareAnnotation}}
    \code{\link{getDefaultProcessingOptions}},
    \code{\link{getPipelineOptions}}, 
}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{

if (isRCloud()) { # disabled on local configs
    # so as not to affect package building process

    # if executed on a local PC, make sure tools
    # are available to the pipeline.
    expfolder = tempdir();
    
    # run the pipeline
    #
    aehts = ArrayExpressHTS("E-GEOD-16190", dir = expfolder);
    
    # load the expression set object
    loadednames = load(paste(expfolder, 
                        "/E-GEOD-16190/eset_notstd_rpkm.RData", sep=""));
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


