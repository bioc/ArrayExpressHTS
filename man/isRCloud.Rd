\name{isRCloud}
\alias{isRCloud}
\title{ Check the code is running on R-Cloud }
\description{
    \code{isRCloud} returns TRUE/FALSE indicating whether the configuration is an R-Cloud.
}
\usage{
    isRCloud()
}
\value{
    If TRUE, the configuration is the R-Cloud.
}
\seealso{
    \code{\link{setPipelineOptions}}, 
    \code{\link{ArrayExpressHTS}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    
    if ( isRCloud() ) {
        
        # we're on the R-Cloud
        print("R-Cloud configuration");
    } else {
    
        # we're somewhere else
        print("Other configuration");
    }
}
