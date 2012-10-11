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
    ArrayExpressHTSFastQ( accession, 
            organism = c("automatic", "Homo_sapiens", "Mus_musculus" ), 
            quality = c("automatic", "FastqQuality", "SFastqQuality"), 
            options = list(
                    stranded          = FALSE,
                    insize            = NULL,
                    insizedev         = NULL,
                    reference         = "genome",
                    aligner           = "tophat",
                    aligner_options   = NULL,
                    count_feature     = "transcript",
                    count_options     = "",
                    count_method      = "cufflinks",
                    filter            = TRUE,
                    filtering_options = NULL),
            usercloud = TRUE,
            rcloudoptions = list(
                    nnodes        = "automatic",
                    pool          = c("4G", "8G", "16G", "32G", "64G"),
                    nretries      = 4),
            steplist = c("align", "count", "eset"),
            dir = getwd(),
            refdir = getDefaultReferenceDir(),
            want.reports = TRUE,
            stop.on.warnings = FALSE )
}
\arguments{
    \item{accession}{ name of the folder where experiment data is kept and computaton data is stored.
        The actual data files and descriptor .sdrf and .idf fles should be stored in the 'data' folder, which
        should be created in this folder. }
    \item{organism}{ defines organism filter, \code{"automatic"} means all organism found in the descript files will be used.
        If a specific organism or an array e.g. c("Homo_sapiens", "Mus_musculus") is defined, those organisms will 
        be used to filter out other unwanted organisms and associated data files from the result object. The value 
        canot be \code{"automatic"} if no .sdrf is found and no automatic discovery can be performed.
    }
    \item{quality}{ this parameter is used to explicitly select the quality scale used in the fastq data files. 
        By default \code{"automatic"} is used. In case quality scale cannot be automatically detected, the user will be 
        prompted to increase the detection depth or set the quality scale manually. "FastqQuality" corresponds 
        to \code{Phred+33}, "SFastqQuality" corresponds to \code{Phred+64}. 
    }
    \item{options}{ defines pipeline options. See \code{\link{getDefaultProcessingOptions}} }
    \item{usercloud}{ defines if the R Cloud will be used to parallel experiment computation. 
        If set to FALSE, the experiment files will be processed sequentially. }
    \item{rcloudoptions}{ defines R Cloud options. See \code{\link{getDefaultRCloudOptions}} }
    \item{steplist}{ defines the steps the pipeline will perform }
    \item{dir}{ folder where experiment data will be stored and processed. Default is current directory. }
    \item{refdir}{ the directory where reference data is located. }
    \item{want.reports}{ defines if quality reports are produced. Reports usually make 
        computation longer and eat up more memory. For faster computation use FALSE. }
    \item{stop.on.warnings}{ self explanatory. Warnings are normally producesd when there 
        are inconsistencies, which however would allow the result to be produced. }
    
}
\value{
     The output is an object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing expression
    values in assayData (corresponding to the raw sequencing data files), 
    the information contained in the .sdrf file in phenoData, the information in the adf file in 
    featureData and the idf file content in experimentData.
    
}
\seealso{
    \code{\link{ArrayExpressHTS}},
    \code{\link{prepareReference}},
    \code{\link{prepareAnnotation}},
    \code{\link{prepareAnnotation}}
    \code{\link{getDefaultProcessingOptions}},
    \code{\link{getPipelineOptions}}, 
}
    
\author{
    Andrew Tikhonov
    Angela Goncalves
    Maintainer: <andrew@ebi.ac.uk>  
    Maintainer: <angela.goncalves@ebi.ac.uk>  
}
\examples{

if (isRCloud()) { # disabled on local configs 
    # so as not to affect package building process
    
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
    # to the pipeline.
    # 
    
    # create a temporary folder where experiment will be copied
    # computing experiment in the package folder may cause issues
    # with file permissions and therefore failures.
    # 
    #
    srcfolder <- system.file("expdata", 
                "testExperiment", package="ArrayExpressHTS");
    
    dstfolder <- tempdir();
    
    file.copy(srcfolder, dstfolder, recursive = TRUE);
    
    # run the pipeline
    #
    # set usercloud = FALSE if executing on local PC, 
    # therefore parallel computation will be disabled
    #
    aehts = ArrayExpressHTSFastQ(accession = "testExperiment", 
        organism = "Homo_sapiens", dir = dstfolder);
        
    # load the expression set object
    loadednames = load(paste(dstfolder, 
                        "/testExperiment/eset_notstd_rpkm.RData", sep=""));
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

