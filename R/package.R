# global defs 
#
.project = list();
.projects = list();
.stop.on.warnings = FALSE;
#.HTML.file = ""


log.info = function(...) {
    cat(paste(date(), getMark(), ... , "\n", sep=""));
}

log.warning = function(...) {
    cat(paste(date(), getMark(), "*** WARNING ***\n", sep=""));
    cat(paste(date(), getMark(), "\n", sep=""));
    cat(paste(date(), getMark(), "    ", ... , "\n", sep=""));
    cat(paste(date(), getMark(), "\n", sep=""));
    
    if (.stop.on.warnings) {
        stop();
    }
}

log.error = function(...) {
    cat(paste(date(), getMark(), "*** ERROR ***\n", sep=""));
    cat(paste(date(), getMark(), "\n", sep=""));
    cat(paste(date(), getMark(), "    ", ... , "\n", sep=""));
    cat(paste(date(), getMark(), "\n", sep=""));
    
    stop();
}


getMark = function() {
    " [AEHTS] ";
}

getDefaultRCloudOptions = function(){
    trace.enter("getDefaultRCloudOptions");
    on.exit({ trace.exit() })
    
    return ( list(
                    nnodes         = "automatic",
                    pool           = "32G",
                    nretries  = 4));
    
}

getMark = function() {
    " [AEHTS] ";
}

cleanupAfterExecution = function() {
    assignStopOnWarnings(FALSE);
}



ArrayExpressHTS <- function( accession,
        options = list (
                stranded          = FALSE,
                insize            = NULL,
                insizedev         = NULL,
                reference         = "genome",
                aligner           = "tophat",
                aligner_options   = NULL,
                count_feature     = "transcript",
                count_options     = "",
                count_method      = "cufflinks",
                filter            = TRUE,
                filtering_options = NULL),
        usercloud = TRUE,
        rcloudoptions = list (
                nnodes     = "automatic",
                pool       = c("4G", "8G", "16G", "32G", "64G"),
                nretries   = 4 ),
        
        dir = getwd(),
        refdir = getDefaultReferenceDir(),
        want.reports = TRUE, 
        stop.on.warnings = FALSE ) {
                             
    trace.enter("ArrayExpressHTS"); 
    on.exit({ trace.exit() })
    
    # add stop on waarnings & cleanup
    assignStopOnWarnings(stop.on.warnings);
    on.exit({ assignStopOnWarnings(FALSE); cleanupAfterExecution() }, add = TRUE);
    
    # merge user and default options
    if (!missing(options)) {
        options = mergeOptions(getDefaultProcessingOptions(), options);
        
        if (!is.null(options$filtering_options)) {
            options$filtering_options = mergeOptions(getDefaultFilteringOptions(), options$filtering_options);
        }
    } else {
        options = getDefaultProcessingOptions();
    }
    
    # merge user and default options
    if (!missing(rcloudoptions)) {
        rcloudoptions = mergeOptions(getDefaultRCloudOptions(), rcloudoptions);
    } else {
        rcloudoptions = getDefaultRCloudOptions();
    }
    
    # check refdir
    if (!file.exists(refdir)) {
        log.info(refdir, " Not Found.");
        stop();
    }
    
    refdir = normalizePath(refdir);
    
    # memory usage
    print.memory.usage("Master Memory Usage: Start")
    
    log.info("creating projects");
    
    resetProjectErrors();
    
    # create projects
    projects <- createAEprojects( accession = accession, dir = dir, refdir = refdir, 
                                  localmode = getPipelineOption("ebilocalmode"), options = options );
    
    checkProjectErrors();
    
    if (length(projects) > 0) { 
        attr(projects, 'metadata')$logfolder = paste(dir, accession, 'cluster-log', sep='/');
        
        assignProjects(projects);
        
        # optimize node usage
        if (rcloudoptions$nnodes == "automatic") {
            rcloudoptions$nnodes = length(.projects);
        }
        
        ecount = runProjects( .projects, usercloud = usercloud, rcloudoptions = rcloudoptions,
                want.reports = want.reports, filter = options$filter );
        
    } else {
        log.info("No data found to compute.");
        ecount = NULL;
    }
    
    # memory usage
    print.memory.usage("Master Memory Usage: End")

    log.info( "Computation Completed" );

    return(ecount);
}

ArrayExpressHTSFastQ <- function( accession, 
        organism = c("automatic", "Homo_sapiens", "Mus_musculus" ), 
        quality = c("automatic", "FastqQuality", "SFastqQuality"), 
        options = list(
                stranded          = FALSE,
                insize            = NULL,
                insizedev         = NULL,
                reference         = "genome",
                aligner           = "tophat",
                aligner_options   = NULL,
                count_feature     = "transcript",
                count_options     = "",
                count_method      = "cufflinks",
                filter            = TRUE,
                filtering_options = NULL),
        usercloud = TRUE,
        rcloudoptions = list(
                nnodes        = "automatic",
                pool          = c("4G", "8G", "16G", "32G", "64G"),
                nretries      = 4),
        dir = getwd(),
        refdir = getDefaultReferenceDir(),
        want.reports = TRUE,
        stop.on.warnings = FALSE ) {
                                  
    trace.enter("ArrayExpressHTSFastQ");
    on.exit({ trace.exit() })
    
    # add stop on waarnings & cleanup
    assignStopOnWarnings(stop.on.warnings);
    on.exit({ assignStopOnWarnings(FALSE); cleanupAfterExecution() }, add = TRUE);
    
    # merge user and default options
    if (!missing(options)) {
        options = mergeOptions(getDefaultProcessingOptions(), options);
        
        if (!is.null(options$filtering_options)) {
            options$filtering_options = mergeOptions(getDefaultFilteringOptions(), options$filtering_options);
        }
    } else {
        options = getDefaultProcessingOptions();
    }
    
    if (!missing(rcloudoptions)) {
        rcloudoptions = mergeOptions(getDefaultRCloudOptions(), rcloudoptions);
    } else {
        rcloudoptions = getDefaultRCloudOptions();
    }
    
    # set default quality
    if (!missing(organism)) {
        organism = "automatic";
    }
    
    # set default quality
    if (!missing(quality)) {
        quality = "automatic";
    }
    
    if (!file.exists(refdir)) {
        log.info(refdir, " Not Found.");
        stop();
    }
    
    refdir = normalizePath(refdir);
    
    
    # memory usage
    print.memory.usage("Master Memory Usage: Start")
    
    log.info( "creating projects" );
    
    resetProjectErrors();
        
    # create projects
    projects <- createFastQProjects( accession = accession, organism = organism, quality = quality, 
            dir = dir, refdir=refdir, options = options )
    
    checkProjectErrors();
    
    if (length(projects) > 0) {
    
        attr(projects, 'metadata')$logfolder = paste(dir, accession, 'cluster-log', sep='/');
        
        assignProjects(projects);
        
        # optimize node usage
        if (rcloudoptions$nnodes == "automatic") {
            rcloudoptions$nnodes = length(.projects);
        }
        
        ecount = runProjects( .projects, usercloud = usercloud, rcloudoptions = rcloudoptions, 
                want.reports = want.reports, filter = options$filter);
        
    } else {
        log.info("No data found to compute.");
        ecount = NULL;
    }   

    # memory usage
    print.memory.usage("Master Memory Usage: End")

    log.info( "Computation Completed" );

    return(ecount);
}

resetProjectErrors <- function(){
    
    # reset missing fastq files
    setPackageVariable("missing-fastq" = FALSE);
    
    # add more as necessary
    # ...
    #
    #
}


checkProjectErrors <- function(){
    
    # check missing fastq files
    #
    missingfastq = getPackageVariable("missing-fastq")
    
    if (!is.null(missingfastq) && missingfastq == TRUE) {
        log.info("ERROR One or more Fastq files are missing!")
        stop();
    }
    
    # add more as necessary
    # ...
    #
    #
}

assignProjects <- function(projects){
    assignInNamespace('.projects', projects, ns="ArrayExpressHTS");
}

assignOneProject <- function(project){
    assignInNamespace('.project', project, ns="ArrayExpressHTS");
}

assignStopOnWarnings <- function(stop.on.warnings){
    assignInNamespace('.stop.on.warnings', stop.on.warnings, ns="ArrayExpressHTS");
}

runProjects <- function( .projects, usercloud = TRUE, rcloudoptions = rcloudoptions, 
        want.reports = FALSE, filter = TRUE ) {
    
    trace.enter("runProjects");
    on.exit({ trace.exit() })
    
    # keep current folder
    currentfolder = getwd();
    on.exit({try( setwd(currentfolder), silent = TRUE )}, add = TRUE);
    
    # memory usage
    print.memory.usage("Memory: runProjects step 1")
    
    
    if( usercloud ) {
        # cluster cleanup
        on.exit({try( cleanupCluster(), silent = TRUE )}, add = TRUE);
        
        prepareCluster(rcloudoptions, logfolder = attr(.projects, 'metadata')$logfolder);
        
        # populate projects data
        log.info( "populating project data" );
        sink <- clusterApply(getCluster(), 1:length(getCluster()), ArrayExpressHTS:::setProjectData, .projects);
        
        # final coutdown
        log.info( "final countdown" );
        clusterEvalQ(getCluster(), ArrayExpressHTS:::finalCountdown());
        
        log.info( "cluster started" );
        result = try( clusterApplyLB(getCluster(), .projects, ArrayExpressHTS:::processOneProject, 
                        .projects, want.reports, .stop.on.warnings) ); 
        
        if (inherits(result, "try-error")) { 
            log.info(" Error during the processing");
        }
        
    } else {
        log.info( "setting environment" );
        
        # setup environment
        sink <- initEnvironmentVariables();
        
        for( .project in .projects ) {
            log.info( "working project ", .project$name );
            
            result = try( processOneProject(.project, .projects, want.reports, .stop.on.warnings) ); 
            
            # memory usage
            print.memory.usage("Memory: runProjects step 2")
            
            if (inherits(result, "try-error")) { 
                log.info(.project$name, " Error during the processing");
            }
        }
    }
                                        
    # memory usage
    print.memory.usage("Memory: runProjects step 3")
    
    # restore current folder                                   
    setwd(currentfolder);
    
    # collect all data into an ExpressionSet object
    #log.info( "building FPKM expression set" );
    #efpkm = to_expressionset( "fpkm", filter = filter )
        
    #log.info( "building COUNT expression set" );
    #ecount = to_expressionset( "count", filter = filter )
    
    log.info( "building expression set" );
    ecount = to_expressionset( filter = filter )

    # memory usage
    print.memory.usage("Memory: runProjects step 4")
    
    log.info( "plotting compared report" );
    
    if (want.reports) {
        if( length(.projects) > 1 )
            plot_compared_raw_report(.projects[[1]])
        else 
            log.info( "Compared report omitted for project with only one sample" )
    }
        
    # memory usage
    print.memory.usage("Memory: runProjects step 5")
        
    return(ecount);
}

# cluster variable
#
cluster001 <- NULL;

# cluster getter
#
getCluster <- function() {
    cluster001;
}

# cluster setter
#
setCluster <- function(cluster) {
    assignInNamespace('cluster001', cluster, ns="ArrayExpressHTS");
}

# prepare cluster
#
prepareCluster <- function(rcloudoptions, logfolder) {
    
    trace.enter("prepareCluster");
    on.exit({ trace.exit() })
    
    if ( !isRCloud() ) {
        log.info('RCLOUD cluster is not available. Please run from R-Cloud Workbench.');
        log.info('http://www.ebi.ac.uk/Tools/rcloud');
        stop();
        return();
    }
    
    # convert, in case a string is provided
    rcloudoptions$nnodes   = as.integer(rcloudoptions$nnodes);
    rcloudoptions$nretries = as.integer(rcloudoptions$nretries);
    
    foundCluster = FALSE;
    
    if (!is.null(getCluster())) {
        
        log.info( "cluster found" );
        log.info( "cleaning.." );
        
        cleanupCluster();
    }
    
    log.info( "creating new cluster" );
    
    finalcluster = list();
    
    for (attempt in 1:rcloudoptions$nretries) {
        
        # allocate cluster
        #
        cl0 = makeCluster(rcloudoptions$nnodes - length(finalcluster), 
                type="RCLOUD", pool = rcloudoptions$pool);
        
        finalcluster = mergeClusters(cl0, finalcluster);
        
        if (length(finalcluster) == rcloudoptions$nnodes) {
            # cluster created successfully
            #
            break;
        } else {
            log.info( attempt, " attempt to create a cluster" );
        
            if (attempt == rcloudoptions$nretries) {
            
                if (length(finalcluster) == 0) {
                    
                    log.info("Error: computing cluster could not be created");
                    stop();
                
                } else {
                    
                    log.warning(length(finalcluster), " nodes of ", 
                            rcloudoptions$nnodes," were created.");
                    break;
                }
                
            } else {
                minutes.to.wait = ceiling(runif(1, 5, 15))
                
                log.info( "Waiting ", minutes.to.wait ," minutes before next attempt" );
                
                Sys.sleep(minutes.to.wait * 60);
            }
        }
    }
    
    # add attributes
    #
    attr(finalcluster, "pool") = rcloudoptions$pool;
    
    # store it
    #
    setCluster(finalcluster);
    
    # delay in case we allocated a server 
    # that is being initialized, let it finish
    #
    
    delaytime = 10;
    
    log.info( delaytime, " seconds cluster init wait" );
    
    Sys.sleep(delaytime);
    
    log.info( "starting cluster logs" );
    
    # install output files
    #
    dir.create(logfolder, showWarnings = FALSE);
    
    log.info("log folder = ", logfolder );
    
    sink <- clusterApply(getCluster(), 1:length(getCluster()), 
            ArrayExpressHTS:::createServerLog, logfolder, 'log');
    
    clusterEvalQ(getCluster(), options('outfile0'));
    
    # display main server name
    #
    log.info( "cluster master ", get('.PrivateEnv')$getServerName() );
    
    log.info( "cluster nodes" );

    # diaply cluster and node names
    #
    print( clusterEvalQ(getCluster(), c( get('.PrivateEnv')$getServerName(), Sys.info()[4] )) )
    
    log.info( "initializing cluster" );
    
    # run .jinit()
    #
    sink <- clusterEvalQ(getCluster(), .jinit());
    
    # propagate global to the cluster options
    #
    log.info( "loading global options" );
    
    opt001 = options();
    
    sink <- clusterApply(getCluster(), 1:length(getCluster()), function(i, x){ options(x) }, opt001)
    
    log.info( "loading libraries on cluster nodes" );

    # detach the package
    #
    sink <- clusterEvalQ(getCluster(), 
        if (("package:ArrayExpressHTS" %in% search())) { detach('package:ArrayExpressHTS', unload=TRUE) } )
    
    # load the library
    #
    sink <- clusterEvalQ(getCluster(), library( "ArrayExpressHTS" ));
    
    optionsToPrint = c(
        "ArrayExpressHTS.fasta_formatter",
        "ArrayExpressHTS.cufflinks",
        "ArrayExpressHTS.samtools",
        "ArrayExpressHTS.bwa",
        "ArrayExpressHTS.mmseq",
        "ArrayExpressHTS.bam2hits",
        "ArrayExpressHTS.bowtie",
        "ArrayExpressHTS.tophat")
    
    clusterApply(getCluster(), 1:length(getCluster()), 
        function(i, opts){ for(o in opts) { message("Global ", o, " = ", 
            getOption(o));} }, optionsToPrint)
    
    log.info( "loading pipeline options" );
    
    opt002 = getPipelineOptions();
    
    sink <- clusterApply(getCluster(), 1:length(getCluster()), 
        function(i, x){ ArrayExpressHTS:::assignPipelineOptions(x) }, opt002)

    clusterApply(getCluster(), 1:length(getCluster()), 
        function(i, opts){ for(o in opts) { message("Package ", o, " = ", 
            ArrayExpressHTS:::getPipelineOption(o));} }, optionsToPrint)
    
    
    #log.info( "setting environment" );
    # setup environment
    #sink <- clusterEvalQ(getCluster(), ArrayExpressHTS:::initEnvironmentVariables())

}

cleanupCluster <- function() {
    trace.enter("cleanupCluster");
    on.exit({ trace.exit() })
    
    log.info( "cleaning up cluster" );
    
    try(stopCluster(getCluster()), silent = TRUE);
    try(setCluster(NULL), silent = TRUE)
    try(cleanupClusters(), silent = TRUE);
}

createServerLog <- function(x, logfolder, name) { 
    trace.enter("createServerLog");
    on.exit({ trace.exit() })
    
    path = paste(logfolder, paste(name, x, sep='-'), sep='/'); 
    sinkWorkerOutput(path); 
    options('outfile0' = path) 
}

setProjectData <- function(unused, projects) {
    trace.enter("setProjectData");
    on.exit({ trace.exit() })
    
    assignProjects(projects)
}

processOneProject <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProject");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); cleanupAfterExecution(); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects)
    assignOneProject(project)
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name );
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 1")
    
    #.project$qual_type = class(quality( seq [[1]] )[1])[1]
    
    if (want.reports) {
        
        seq = try( fastq_to_shortreadq() )
        if (inherits(seq, "try-error")) { 
            log.info(" Error reading data file");
            return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 2")
        
        tab = try( shortread_to_tab( seq ) );
        if (inherits(tab, "try-error")) { 
            log.info(" Error converting data file");
            return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 3")
        
        # make sure there's at least one device there
        #
        plottingresult = try( plot_raw_report( seq, tab ) )
        if (inherits(plottingresult, "try-error")) { 
            log.info(" Error plotting raw report");
            # continue even if it fails
        }
    
        # memory usage
        print.memory.usage("Memory: processOneProject step 4")
    }
    
    alignmentresult = try( align( sure.i.want.to.run.this = TRUE ) )
    if (inherits(alignmentresult, "try-error")) { 
        log.info(" Error aligning data");
        return();
    }
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 5")
    
    #estimationresult = try( call_cufflinks( run = TRUE ) )
    estimationresult = try( call_estimator( run = TRUE ) )
    if (inherits(estimationresult, "try-error")) { 
        log.info(" Error estimating expression");
        return();
    }
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 6")
    
    
    #bam = try( bam_to_list( return = FALSE ) )
    #if (inherits(bam, "try-error")) { 
    #    log.info(" Error converting bam to list");
        # continue
    #}
    
    if (want.reports) {
        plottingresult = try( plot_aligned_report( seq ) )
        
        if (inherits(plottingresult, "try-error")) { 
            log.info(" Error plotting alignment report");
            # continue
        }
    
        # memory usage
        print.memory.usage("Memory: processOneProject step 7")
    
        #rm(bam);
        rm(seq);
    
        # memory usage
        print.memory.usage("Memory: processOneProject step 8")
    }
    
    # call it here to pre-cache 
    # when we're still in parallel
    total_reads = try( number_raw_reads( ) )
    if (inherits(total_reads, "try-error")) { 
        log.info(" Error getting number of reads");
        # continue
    }
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 9")
    
    # call it here to pre-cache 
    aligned_reads = try( number_aligned_reads() )
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 10")

    gc();
    log.info( "done" );
}

isRCloud <- function() {
    trace.enter("isRCloud");
    on.exit({ trace.exit() })
    
    # exists('makeRCLOUDcluster')
    return( exists('.PrivateEnv') );
}


finalCountdown <- function(x){ 
    trace.enter("finalCountdown");
    on.exit({ trace.exit() })
    
    log.info( "final countdown" );
}

# organism = "Homo_sapiens", version="GRCh37.60", type="genome"/"transcriptome", location="/ebi/microarray/home/biocep/sequencing"; aligner="bowtie"
prepareReference <- function( organism, version = "current", type = c("genome", "transcriptome"), 
        location = getDefaultReferenceDir(), aligner = c("bwa", "bowtie", "tophat"), refresh = FALSE, run = TRUE ) {
        
    trace.enter("prepareReference");
    on.exit({ trace.exit() })
    
    if(missing(type)) {
        log.info("'type' is not defined");
        stop();
    };
    
    if(missing(aligner)) {
        log.info("'aligner' is not defined");
        stop();
    };
    
    if (!file.exists(location)) {
        log.info(location, " Not Found. Please create it.");
        stop();
    }
    
    location = normalizePath(location);
    
    version = getEnsemblReference( organism = organism, type = type, 
        version = version, location = location, refresh = refresh, run = run)
    
    indexReference( organism = organism, type = type, version = version, 
        location = location, aligner = aligner, refresh = refresh, run = run)
    
    return( version )
}

# organism = "Homo_sapiens"  version="GRCh37.60" location="/ebi/microarray/home/biocep/sequencing" aligner="bowtie"
prepareAnnotation <- function( organism, version = "current", 
        location = getDefaultReferenceDir(), refresh = FALSE, run = TRUE ) {
    
    trace.enter("prepareAnnotation");
    on.exit({ trace.exit() })
    
    if (!file.exists(location)) {
        log.info(location, " Not Found. Please create it.");
        stop();
    }
    
    location = normalizePath(location);
    
    version = getEnsemblAnnotation( organism = organism, version = version, 
        location = location, refresh = refresh, run = run )
    
    return( version )
}

