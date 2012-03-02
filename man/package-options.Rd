\name{package-options}
\title{Functions to work with package level options}
\alias{getPipelineOptions}
\alias{getPipelineOption}
\alias{setPipelineOptions}
\description{
    Functions to work with ArrayExpressHTS package level options.
}
\usage{
    getPipelineOptions()
    getPipelineOption(name, default = NULL, options = defaultOptions)
    setPipelineOptions(...)
}
\arguments{
    \item{...}{ a list ot options to set. }
    \item{name}{ name of the option to return. In order to see all options, 
        use \code{\link{getPipelineOptions}}; }
    \item{default}{ default values in case the option is not defined }
    \item{options}{ options environment, by default the ArrayExpresHTS 
        environment; }
}
\details{
    These are the basic functions to get and set pakcage options and 
    perform necessary reinitialization following the change in required.
    
    \code{getPipelineOptions} returns an environment where ArrayExpressHTS 
    options are stored.
    
    \code{setPipelineOptions} sets one or a number of ArrayExpressHTS 
    options. 
    
    Please note that the following options are used to store locations of 
    external tools used by the package: 
    \code{"ArrayExpressHTS.bowtie"}, \code{"ArrayExpressHTS.tophat"}, 
    \code{"ArrayExpressHTS.bwa"}, \code{"ArrayExpressHTS.cufflinks"}, 
    \code{"ArrayExpressHTS.mmseq"}, \code{"ArrayExpressHTS.samtools"}, etc. 
    
    If any of the options needs to be changed, use \code{setPipelineOption(...)} 
    to make necessary chnages. Following the change, PATH environment varibale 
    automatically gets updated as necessary.
    
    Please note that the package supports only certain versions of tools. 
    Other versions may not be fully compatible, and if used, this can result 
    in errors on certain steps of the pipeline. The list of supported 
    versions can be obtained using \code{getPieplineOption(...)} with the 
    corresponding tool option, e.g. \code{"ArrayExpressHTS.bowtie"}. Not 
    supported versions can be used at user's own risk.
    
    
    \code{getPipelineOption} returns ArrayExpressHTS option specified by 
    option name. If option is not defined the specified default value will 
    be returned.
}

\value{
    The output is the environment where pipeline options are stored.
}
\seealso{
    \code{\link{getDefaultProcessingOptions}}, 
    \code{\link{ArrayExpressHTS}},
    \code{\link{ArrayExpressHTSFastQ}}
}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    # get options
    pipelineOptions = getPipelineOptions()
    
    # list all names
    ls(pipelineOptions)
    
    # set pipeline options
    setPipelineOptions("trace" = "disabled", "memorymonitor" = "disabled");
    
    # get options
    getPipelineOption("trace");
    getPipelineOption("memorymonitor");
    
}

