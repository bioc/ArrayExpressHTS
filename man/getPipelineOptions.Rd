\name{getPipelineOptions}
\alias{getPipelineOptions}
\title{ Get ArrayExpressHTS options }
\description{
    \code{getPipelineOptions} returns an environment where ArrayExpressHTS options are stored.
}
\usage{
    getPipelineOptions()
}
\value{
    The output is the environment where pipeline options are stored.
}
    \seealso{\code{\link{getPipelineOption}}, \code{\link{setPipelineOptions}}, \code{\link{initDefaultEnvironment}}, \code{\link{ArrayExpressHTS}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    # get options
    pipelineOptions = getPipelineOptions()
    
    # list all names
    ls(pipelineOptions)
}
