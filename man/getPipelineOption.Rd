\name{getPipelineOption}
\alias{getPipelineOption}
\title{ Get ArrayExpressHTS option }
\description{
    \code{getPipelineOption} returns ArrayExpressHTS option specified by option name.
}
\usage{
    getPipelineOption(name, options = defaultOptions)
}
\arguments{
    \item{name}{ name of the option to return. In order to see all options, use \code{\link{getPipelineOptions}}; }
    \item{options}{ options environment, by default the ArrayExpresHTS option environment; }
}
\value{
    The result is the value of the option.
}
\seealso{
    \code{\link{getPipelineOptions}}, 
    \code{\link{setPipelineOptions}}, 
    \code{\link{initDefaultEnvironment}}, 
    \code{\link{ArrayExpressHTS}}}
\author{
    
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    
    getPipelineOption("trace");
    getPipelineOption("cufflinks");
    getPipelineOption("bowtie");
    getPipelineOption("tophat");
    getPipelineOption("samtools");
}
