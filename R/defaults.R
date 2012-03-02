getDefaultReferenceDir <- function() {
    trace.enter("getDefaultReferenceDir");
    on.exit({ trace.exit() })
    
    optval = getPipelineOption('ArrayExpressHTS.reference');
    
    for (path in optval) {
        if (file.exists(path)) {
            return(path);
        }
    }
    
    return("");
}

getOptionPathComponent <- function(optname) {
    optvalue = getPipelineOption(optname)
    
    if( is.null(optvalue) || nchar(optvalue) == 0) {
        return("");
    } else {
        return(paste(optvalue, ":", sep=""));
    }
}

initEnvironmentVariables <- function() {

    trace.enter("initEnvironmentVariables");
    on.exit({ trace.exit() })
    
    addPathToPATH("/ebi/microarray/home/biocep/local/lib64/R/lib");
    addPathToPATH("/ebi/microarray/home/biocep/local/lib64/R/bin");
    addPathToPATH("/ebi/microarray/home/biocep/local/bin");
    addPathToPATH("/ebi/microarray/sw/bin");
    addPathToPATH("/ebi/research/software/Linux_x86_64/opt/java/jdk1.6/bin");
    addPathToPATH("/usr/local/bin");
    addPathToPATH("/bin");
    addPathToPATH("/usr/bin");
    addPathToPATH("/usr/X11R6/bin");
    addPathToPATH("/usr/lib64");
    addPathToPATH("/lib64");
    
    addPathToPATH(getPipelineOption("ArrayExpressHTS.cufflinks"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.samtools"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.bwa"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.mmseq"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.bam2hits"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.bowtie"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.tophat"));
    addPathToPATH(getPipelineOption("ArrayExpressHTS.fasta_formatter"));
    
}

removePathFromPATH = function(pathtoremove) {
    
    # get the array of existing records
    pathArray = unlist(strsplit(Sys.getenv('PATH'), ":"));
    
    # get all but the path
    pathArray = pathArray[ is.na(match(pathArray, pathtoremove)) ];
    
    # assemble back and set the PATH
    Sys.setenv('PATH' = paste(pathArray, collapse=":"));
}


addPathToPATH = function(newpath) {
    
    # get the array of existing records
    pathArray = unlist(strsplit(Sys.getenv('PATH'), ":"));
    
    # get all but the path
    # make sure the path newpath is not duplicated
    pathArray = pathArray[ is.na(match(pathArray, newpath)) ];
    
    # assemble back, add the path at the 
    # beginning and set the PATH variable
    # 
    Sys.setenv('PATH' = paste(newpath, paste(pathArray, collapse=":"), sep=":"));
}


replacePathInPATH = function(oldpath, newpath) {
    
    # get the array of existing records
    pathArray = unlist(strsplit(Sys.getenv('PATH'), ":"));
    
    # get all but the oldpath
    # remove duplicates of oldpath
    pathArray = pathArray[ is.na(match(pathArray, oldpath)) ];
    
    # get all but the newpath
    # remove duplicates of newpath
    pathArray = pathArray[ is.na(match(pathArray, newpath)) ];
    
    # assemble back, add the path at the 
    # beginning and set the PATH variable
    # 
    Sys.setenv('PATH' = paste(newpath, paste(pathArray, collapse=":"), sep=":"));
}

