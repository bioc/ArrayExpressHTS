getDefaultReferenceDir <- function() {
    trace.enter("getDefaultReferenceDir");
    on.exit({ trace.exit() })
    
    return(getPipelineOption('referenceDir'));

}

initDefaultEnvironment <- function() {

    trace.enter("initDefaultEnvironment");
    on.exit({ trace.exit() })
    
    Sys.setenv('PATH' = 
        paste(
            getPipelineOption('pathbackup'),
        
            "/ebi/microarray/home/biocep/local/lib64/R/lib:",
            "/ebi/microarray/home/biocep/local/lib64/R/bin:",
            "/ebi/microarray/home/biocep/local/bin:",
            "/ebi/microarray/sw/bin:",
            "/ebi/research/software/Linux_x86_64/opt/java/jdk1.6/bin:",
            "/usr/local/bin:",
            "/bin:",
            "/usr/bin:",
            "/usr/X11R6/bin:",
            "/usr/lib64:",
            "/lib64:",
        
            getPipelineOption('cufflinks'), ":",
            getPipelineOption('samtools'), ":",
            getPipelineOption('bwa'), ":",
            getPipelineOption('mmseq'), ":",
            getPipelineOption('tophat'), ":",
            getPipelineOption('bowtie'), ":",
            getPipelineOption('fastx'), ":",
        
        sep=""));
    
    Sys.setenv("LD_LIBRARY_PATH" = paste(
        "/ebi/microarray/home/biocep/local/tools/samtools-0.1.8:", Sys.getenv("LD_LIBRARY_PATH"), sep=""))

}
