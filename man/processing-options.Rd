\name{processing-options}
\title{Functions to work with ArrayExpressHTS processing options}
\alias{getDefaultFilteringOptions}
\alias{getDefaultProcessingOptions}
\alias{getDefaultRCloudOptions}
\description{
    Functions to get default ArrayExpressHTS processing options.
}
\usage{
    getDefaultFilteringOptions()
    getDefaultProcessingOptions()
    getDefaultRCloudOptions()
}
\arguments{
    No arguments.
}
\value{
    These functions return default options used for ArrayExpressHTS 
    processing.
    
    Function \code{getDefaultFilteringOptions} returns a \code{list} of 
    filtering options. These options control filtering of BAM files after 
    alignemnt.
    
    \item{\code{mismatches}}{ \code{2} by default, number of mismatches
    }
    \item{\code{chr_ignore}}{ \code{c("MT")} by default, reference names to ignore
    }
    \item{\code{minqual}}{ \code{10} by default
    }
    \item{\code{minmapq}}{ \code{1} by default
    }
    \item{\code{maxN}}{ \code{2} by default
    }
    \item{\code{maxpol}}{ \code{0.75} by default, ither a proportion (<0) or an integer
    }
    \item{\code{duplicates}}{ \code{"remove"} by default, supported value: \code{"keep"}, \code{"remove"}
    }
    \item{\code{multihits}}{ \code{"remove"} by default, supported value: \code{"keep"}, \code{"remove"}
    }
    \item{\code{gapped}}{ \code{"remove"} by default, supported value: \code{"keep"}, \code{"remove"}
    }
    ---
    
    Function \code{getDefaultProcessingOptions} returns a \code{list} of 
    options used throught the processing. These options control alignment and 
    counting.
    
    These options include:
    \item{\code{stranded}}{ \code{FALSE} by default, set to \code{TRUE} if a strand 
        specific protocol was used 
    }
    \item{\code{insize}}{ \code{NULL} by default, the insert size, an integer, 
        which will be automatically determined if set to \code{NULL}
    }
    \item{\code{insizedev}}{ \code{NULL} by default, insert size deviation, an integer, 
        which will be automatically determined if set to \code{NULL}
    }
    \item{\code{reference}}{ \code{"genome"} by default, controlls whether alignment 
        is performed using \code{"genome"} or \code{"transcriptome"} reference.
    }
    \item{\code{aligner}}{ \code{"tophat"} by default, supported values are: \code{"tophat"}, 
        \code{"bowtie"}, \code{"bwa"} or \code{"custom"}
    }
    \item{\code{aligner_options}}{ \code{NULL} by default, string of options to be 
        passed to the aligner according to manual.
    }
    \item{\code{count_feature}}{ \code{"transcript"} by default, count over \code{"gene"} or \code{"transcript"}
    }
    \item{\code{count_method}}{ \code{"cufflinks"} by default, supported: \code{"cufflinks"}, \code{"mmseq"} or 
        \code{"count"}
    }
    \item{\code{count_options}}{ options for count software
    }
    \item{\code{standardise}}{ \code{FALSE} by default, supported: \code{TRUE}, \code{FALSE} 
    }
    \item{\code{normalisation}}{ \code{"rpkm"} by default, supported: \code{"none"} or \code{"rpkm"} 
    }
    \item{\code{filter}}{ \code{"TRUE"} by defult, supported: \code{"TRUE"}, \code{"FALSE"}
    }
    \item{\code{filtering_options}}{ these are the default filtering options, see above.
    }
    ---
    
    Function \code{getDefaultRCloudOptions} returns a set of R Cloud options
    related to R Cloud resources and execution flow.
    
    These include:
    \item{\code{nnodes}}{ \code{"automatic"} by default, supported: \code{"automatic"} or a numeric value, e.g. \code{10} 
    }
    \item{\code{pool}}{ \code{"32G"} by default, supported: \code{"4G"}, \code{"8G"}, \code{"16G"}, \code{"32G"}, \code{"64G"}
    }
    \item{\code{nretries}}{ \code{4} by default, a numberic value.
    }
}
\seealso{
    \code{\link{getPipelineOptions}}, 
    \code{\link{ArrayExpressHTS}},
    \code{\link{ArrayExpressHTSFastQ}}
}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    
    # get filtering options
    getDefaultFilteringOptions()
    
    # get alignment and counting options
    getDefaultProcessingOptions()
    
    # get R Cloud options
    getDefaultRCloudOptions()

}

