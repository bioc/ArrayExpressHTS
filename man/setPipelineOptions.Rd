\name{setPipelineOptions}
\alias{setPipelineOptions}
\title{ Set ArrayExpressHTS options }
\description{
    \code{setPipelineOptions} sets one or a number of ArrayExpressHTS options.
}
\usage{
    setPipelineOptions(...)
}
\arguments{
    \item{...}{ a list ot options to set. }
}
\value{
    No output.
}
\seealso{
    \code{\link{getPipelineOption}}, 
    \code{\link{getPipelineOptions}}, 
    \code{\link{initEnvironmentVariables}}, 
    \code{\link{ArrayExpressHTS}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    setPipelineOptions("trace" = "disabled", "memorymonitor" = "disabled");
}
