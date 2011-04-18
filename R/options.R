#
#   O P T I O N S 
#
#

defaultOptions <- NULL

initPipelineOptions <- function() {

    options <- list(
        aligner       = "tophat",
        trace         = "disabled",
        memorymonitor = "disabled",
        ebilocalmode  = FALSE,
        cufflinks = "/ebi/microarray/home/biocep/local/tools/cufflinks-0.8.2.Linux_x86_64",
        samtools  = "/ebi/microarray/home/biocep/local/tools/samtools-0.1.8",
        bwa       = "/ebi/microarray/home/biocep/local/tools/bwa-0.5.7",
        mmseq     = "/ebi/microarray/home/biocep/local/tools/mmseq_0.9.8",
        tophat    = "/ebi/microarray/home/biocep/local/tools/tophat-1.0.14.Linux_x86_64",
        bowtie    = "/ebi/microarray/home/biocep/local/tools/bowtie-0.12.5",
        fastx     = "/ebi/microarray/home/biocep/local/tools/fastx-toolkit-0.0.13",
        referenceDir  = "/ebi/microarray/home/biocep/service/sequencing",
        pathbackup    = Sys.getenv('PATH'));
    
    assignInNamespace('defaultOptions', 
        addOptions(emptyenv(), options), ns="ArrayExpressHTS");
    
}

addOptions <- function(options, new) {
    if (! is.null(new)) {
        options <- new.env(parent = options)
        names <- names(new)
        for (i in seq(along = new))
        assign(names[i], new[[i]], env = options)
    }
    options
}

getPipelineOptions <- function() {
    return(defaultOptions);
}

getPipelineOption <- function(name, options = defaultOptions) {
    get(name, env = options)
}

setPipelineOptions <- function(...) {
    list <- list(...)
    names <- names(list)
    for (i in seq(along = list))
    assign(names[i], list[[i]], env = defaultOptions)
}


