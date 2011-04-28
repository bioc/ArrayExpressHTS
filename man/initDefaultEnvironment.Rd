\name{initDefaultEnvironment}
\alias{initDefaultEnvironment}
\title{ Initialize ArrayExpressHTS default environment }
\description{
    \code{initDefaultEnvironment} initializes 'PATH' environment variable with paths to the tools used by ArrayExpressHTS.
    The function can be called to reinitialize the environment if any of the paths needs to be changed or if the package 
    is executed in a local PC configuration. 
    
    Paths to pipeline tools are stored in the following pipeline options 'bowtie', 'tophat', 'bwa', 'cufflinks', 'mmseq', 
    'samtools', etc. If you need to change any of the paths, use \code{\link{setPipelineOptions}}. To check current values, 
    use \code{\link{getPipelineOption}}. To see all options, use \code{\link{getPipelineOptions}}.
    
    When options, namely the paths, have been changed call 'initDefaultEnvironment' to reinitialize environment. If the
    pipeline is executed on the R-Cloud, there is generally no necessity to setup anything since everything has been 
    setup to a proper working configuration.
    
    Plase note that some 'tophat' and 'cuflinks' versions are not fully compatible and cannot be automatically executed 
    by the pipeline without intervension. Use the versions specified in the options, these versions are used in the R-Cloud.
    Other versions can be used at user''s own risk.
    
}
\usage{
    initDefaultEnvironment()
}
\arguments{
    No arguments.
}
\value{
    No output.
}
\seealso{
    \code{\link{getPipelineOption}}, 
    \code{\link{getPipelineOptions}}, 
    \code{\link{setPipelineOptions}}, 
    \code{\link{ArrayExpressHTS}}}
\author{
    Andrew Tikhonov <andrew@ebi.ac.uk>, Angela Goncalves <angela.goncalves@ebi.ac.uk>
}
\examples{
    initDefaultEnvironment();
}
