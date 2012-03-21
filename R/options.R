#
#   O P T I O N S 
#
#

defaultOptions <- NULL

initPipelineOptions <- function() {
    
    # create defult options
    #
    options <- list(
        aligner               = "tophat",
        trace                 = "disabled",
        memorymonitor         = "enabled",
        ebilocalmode          = FALSE,
        fastqreadmax          = 10000,
        ignorequalityerrors   = FALSE,
        defaultquality        = "FastqQuality",
        defaultquality        = "",
        pathbackup            = Sys.getenv('PATH'));
    
    optionsenv = createOptions(emptyenv(), options);
    
    # define them
    #
    assignPipelineOptions(optionsenv);
    
    # define reference and toolbase
    #
    refoptionlist = list(
        "ArrayExpressHTS.reference" = c("/ebi/microarray/home/biocep/service/sequencing", "/ebi/rcloud/sequencing"));
    
    
    log.info("Setting Options Step 1");
    
    # setup paths to toolbase
    checkAndDefineReferencePath(refoptionlist);
    
    baseoptionlist = list(
        "ArrayExpressHTS.toolbase"  = c("/ebi/rcloud/local/tools", "/ebi/microarray/home/biocep/local/tools"));
    
    
    log.info("Setting Options Step 2");
    
    # setup paths to toolbase
    checkAndDefineBasePath(baseoptionlist);
    
    # experimental versions
    tooloptionlist = list(
            "ArrayExpressHTS.fasta_formatter" = "/ebi/microarray/home/biocep/local/tools/fastx-toolkit-0.0.13",
            "ArrayExpressHTS.cufflinks"       = "/ebi/microarray/home/biocep/local/tools/cufflinks-1.1.0.Linux_x86_64",
            "ArrayExpressHTS.samtools"        = "/ebi/microarray/home/biocep/local/tools/samtools-0.1.18",
            "ArrayExpressHTS.bwa"             = "/ebi/microarray/home/biocep/local/tools/bwa-0.5.9",
            "ArrayExpressHTS.mmseq"           = "/ebi/microarray/home/biocep/local/tools/mmseq_0.9.14",
            "ArrayExpressHTS.bam2hits"        = "/ebi/microarray/home/biocep/local/tools/mmseq_0.9.14",
            "ArrayExpressHTS.bowtie"          = "/ebi/microarray/home/biocep/local/tools/bowtie-0.12.7",
            "ArrayExpressHTS.tophat"          = "/ebi/microarray/home/biocep/local/tools/tophat-1.3.2.Linux_x86_64");
    
    
    log.info("Setting Options Step 3");
    
    checkAndDefineToolOptions(tooloptionlist);
    
    return(optionsenv);
}

getToolOptionNames = function() {
    return (c("ArrayExpressHTS.fasta_formatter",
            "ArrayExpressHTS.cufflinks",
            "ArrayExpressHTS.samtools",
            "ArrayExpressHTS.bwa",
            "ArrayExpressHTS.mmseq",
            "ArrayExpressHTS.bam2hits",
            "ArrayExpressHTS.bowtie",
            "ArrayExpressHTS.tophat"));
}

checkAndDefineReferencePath <- function(refpathlist){
    for(optname in names(refpathlist)) {
        # read from global options
        optvalue = getOption(optname);
        
        if (is.null(optvalue)) {
            # check default values
            defaultvalue = refpathlist[[optname]];
            
            for(path in defaultvalue) {
                if (file.exists(path)) {
                    refpathlist[[optname]] = path;
                    break;
                }
            }
        } else {
            # keep value defined in the options
            #
            refpathlist[[optname]] = optvalue;
        }
        
        # inject options
        addOptions(defaultOptions, refpathlist)
    }
}


checkAndDefineBasePath <- function(basepathlist){
    for(optname in names(basepathlist)) {
        
        # read from global options
        optvalue = getOption(optname);
        
        if (!is.null(optvalue)) {
            basepathlist[[optname]] = c(optvalue, basepathlist[[optname]]);
        }

        for(path in basepathlist[[optname]]) {
            if (file.exists(path)) {
                addPathToPATH(path);
                break;
            }
        }
    }
    
    # inject options
    addOptions(defaultOptions, basepathlist)
}


checkAndDefineToolOptions <- function(tooloptionlist) {
    
    for(optname in names(tooloptionlist)) {
        
        ##optvalue = getOption(optname);
        
        # check default values
        #
        optionpath = getOption(optname, default = tooloptionlist[[optname]])
        
        toolpath = checkOneTool(optname, optionpath);
        toolname = getToolName(optname);
        
        if (!is.null(toolpath)) {
            log.info("Found ", toolpath, "/", toolname);
            
            
            # check for version mismatch
            # report mismatch
            if (toolFolderName(toolpath) != toolFolderName(optionpath)) {
                log.warning(toolpath, " doesn't match to ", optionpath, " defined in the '",optname,"' option");
            }
            
            tooloptionlist[[optname]] = toolpath;
            
        } else {
            log.warning(toolname, " not found");
            log.info("Use options('",optname,"' = '/path/to/",toolname,"') to define the location.")
            log.info("Use .Rprofile to make options persistent.");
            log.info("");
        }
    }
    
    # inject options
    #
    addOptions(defaultOptions, tooloptionlist)
}


getToolName <- function(tooloptionname) {
    substr(tooloptionname, nchar("ArrayExpressHTS.") + 1,1000);
}

checkOneTool <- function(optionname, toolpath) {
    
    toolname = getToolName(optionname);
    
    # check tool exists
    newpath = checkToolExists(toolname, toolpath);
    
    if (is.null(newpath)) {
        
        # check toolbase
        toolbase = unlist(strsplit(Sys.getenv('PATH'), ":"));
        
        if (!is.null(toolbase)) {
            
            newpath = locateInToolBase(toolname, toolpath, toolbase);
            
            if (!is.null(newpath)) {
                return(newpath);
            }
        }
        
        # try which
        newpath = locateInSystem(toolname);
    }
    
    return(newpath);
    
}

checkToolExists <- function(toolname, location) {
    if (file.exists(paste(location, "/", toolname, sep="/"))) {
        return(location);
    }
    return(NULL);
}

toolFolderName <- function(toolpath) {
    splitresult = unlist(strsplit(toolpath, "/"));
    return(splitresult[length(splitresult)]);
}

locateInToolBase <- function(toolname, toolpath, toolbase) {
    
    toolfolder = toolFolderName(toolpath)
    
    for(base in toolbase) {
        newtoolpath = paste(base, "/", toolfolder, sep="");
        
        if (file.exists(paste(newtoolpath, "/", toolname, sep="/"))) {
            return(newtoolpath);
        }
    }
    
    return(NULL);
}

locateInSystem <- function(toolname) {
    
    toolpath = suppressWarnings(system(paste("which ", toolname, sep=""),
        intern=TRUE, ignore.stderr=TRUE));
    
    if (length(toolpath) > 0) {
        #return(toolpath);
        return(substr(toolpath, 1, nchar(toolpath)-nchar(toolname)-1));
    }
    
    return(NULL);
}

createOptions <- function(options, new) {
    if (! is.null(new)) {
        options <- new.env(parent = options)
        names <- names(new)
        for (i in seq(along = new))
        assign(names[i], new[[i]], env = options)
    }
    options
}

addOptions <- function(options, new) {
    if (! is.null(new)) {
        names <- names(new)
        for (i in seq(along = new))
        assign(names[i], new[[i]], env = options)
    }
    options
}


assignPipelineOptions <- function(options){
    assignInNamespace('defaultOptions', options, ns="ArrayExpressHTS");
}


getPipelineOptions <- function() {
    return(defaultOptions);
}

getPipelineOption <- function(name, default = NULL, options = defaultOptions) {
    optval = getOption(name);
    
    if (is.null(optval)) {
        if (name %in% ls(options)) {
            optval = get(name, env = options);
        } else {
            optval = default;
        }
    }
    
    return(optval);
}

setPipelineOption <- function(name, value) {
    
    oldpath = NULL;
    
    # check the option is a path to a tool
    if (name %in% getToolOptionNames()) {
        # check the path exists
        if (!file.exists(value)) {
            log.warning("Path doesn't exist: ", value);
        }

        oldpath = getPipelineOption(name);
    }
    
    assign(name, value, env = defaultOptions)
    
    if (!is.null(oldpath)) {
        replacePathInPATH(oldpath, value);
    }
}


setPipelineOptions <- function(...) {
    list <- list(...)
    names <- names(list)
    for (i in seq(along = list)) {
        setPipelineOption(names[i], list[[i]]);
    }
}

