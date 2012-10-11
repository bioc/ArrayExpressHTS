# global defs 
#
.project = list();
.projects = list();
.stop.on.warnings = FALSE;
#.HTML.file = ""

.submid = NULL;

getDefaultRCloudOptions = function(){
    trace.enter("getDefaultRCloudOptions");
    on.exit({ trace.exit() })
    
    return ( list(
                    nnodes         = "automatic",
                    pool           = "32G",
                    nretries  = 4));
    
}

cleanupAfterExecution = function() {
    assignStopOnWarnings(FALSE);
}

setupSUBM_ID = function() {
    # check if submission id has been defined
    #
    #
    submid = getOption("AEHTS.SUBM_ID");
    
    assignSubmID(submid)
    
    if (!is.null(submid)) {
        #
        #
        log.info("SUBM_ID: ", submid);
    }
}

setupPSRFolder = function( accession, dir, submid, create = TRUE ) {
    
    psrdir = paste( setupBaseFolder( accession, dir, submid ), "PSR", sep="/" );
    
    if (!file.exists( psrdir ) && create) {
        dir.create( psrdir, recursive = TRUE, showWarnings = FALSE );
    }
    
    return( psrdir );
}

setupPSRLogging = function( accession, dir, submid ){
    
    pokePSRExpId(accession);
    
    psrfolder = setupPSRFolder( accession, dir, submid, create = FALSE );
    
    log.info("PSR: ", psrfolder);
    
    createDirAndBackup(psrfolder, recursive = TRUE, showWarnings = FALSE);
    
    pokePSRFolder(psrfolder);
}

makeClusterLogFolder = function( accession, dir, submid ){
    
    logfoldername = "cluster-log";
    
    logfolder = paste(dir, accession, 
            ifelse(is.null(submid), logfoldername, paste(submid, logfoldername, sep="/")), sep="/");
    
    return(logfolder);
}

sendReminder = function(url, issue) {
    from = "rcloud@ebi.ac.uk";
    to   = c("andrew@ebi.ac.uk");
    subject = "Reminder: ENA issue to be addressed";
    body = list(
            paste("ENA request to ", url, " returned unexpected result. Please address.", sep=""), 
            mime_part(issue));
    
    sapply(to, function(x){ 
                sendmail(from, x, subject, body, control=list(smtpServer="smtp.ebi.ac.uk"));
            });
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
        
        steplist = c("align", "count", "eset"),
        
        dir = getwd(),
        refdir = getDefaultReferenceDir(),
        want.reports = TRUE, 
        stop.on.warnings = FALSE ) {
                             
    trace.enter("ArrayExpressHTS");
    on.exit({ trace.exit() });
    
    
    # setup submission ID
    setupSUBM_ID();
    
    # setup folders
    setupBaseFolder(accession, dir, getSubmID());
    setupDataFolder(accession, dir);
    
    # setup PSR logging
    setupPSRLogging( accession, dir, getSubmID() );
    
    
    # add stop on waarnings & cleanup
    assignStopOnWarnings(stop.on.warnings);
    
    # set experiment status
    #
    setExpStatus(STARTED, paste(accession, " started"));
    
    # register step status
    #
    registerExpStepStarted(PARAMETER_SANITY);
    
    # add cleanup and status monitor On Exit
    #
    on.exit({ assignStopOnWarnings(FALSE); 
                cleanupAfterExecution(); 
                updateExpStatusOnExit(); }, add = TRUE);
    
    
    
    # deal with options
    #
    #
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
    #
    if (!file.exists(refdir)) {
        
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                refdir, " Not Found." );
        
        stop();
    }
    
    refdir = normalizePath(refdir);
    
    
    # memory usage
    #
    #
    print.memory.usage("Master Memory Usage: Start")
    
    
    # start creating the workset
    #
    #
    log.info("creating projects");
    
    resetProjectErrors();
    
    # create projects
    projects <- createAEprojects( accession = accession, dir = dir, refdir = refdir, 
                                  localmode = getPipelineOption("ebilocalmode"), options = options );
    
    checkProjectErrors();
    
    # register step status
    #
    registerExpStepCompleted(PARAMETER_SANITY);
    
    
    if (length(projects) > 0) { 
        attr(projects, 'metadata')$logfolder = makeClusterLogFolder(accession, dir, getSubmID());
        
        assignProjects(projects);
        
        # optimize node usage
        if (rcloudoptions$nnodes == "automatic") {
            rcloudoptions$nnodes = length(.projects);
        }
        
        ecount = runProjects( .projects, usercloud = usercloud, rcloudoptions = rcloudoptions,
                want.reports = want.reports, filter = options$filter, steplist = steplist );
        
        # update status to OK
        #
        setExpStatus(COMPUTED, "Completed");

    
    } else {
        log.info("No data found to compute.");
        ecount = NULL;

        # update status to OK
        #
        setExpStatus(NODATA, "No data found to compute.");

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
        
        steplist = c("align", "count", "eset"),
        
        dir = getwd(),
        refdir = getDefaultReferenceDir(),
        want.reports = TRUE,
        stop.on.warnings = FALSE ) {
                                  
    trace.enter("ArrayExpressHTSFastQ");
    on.exit({ trace.exit() })
    
    # setup submission ID
    setupSUBM_ID();
    
    # setup folders
    setupBaseFolder(accession, dir, getSubmID());
    setupDataFolder(accession, dir);
    
    # setup PSR logging
    setupPSRLogging( accession, dir, getSubmID() );
    
    
    # add stop on waarnings & cleanup
    #
    assignStopOnWarnings(stop.on.warnings);
    
    # set experiment status
    #
    setExpStatus(STARTED);
    
    # register step status
    #
    registerExpStepStarted(PARAMETER_SANITY);
    
    # On Exit
    #
    #
    on.exit({ assignStopOnWarnings(FALSE); 
                cleanupAfterExecution(); 
                updateExpStatusOnExit(); }, add = TRUE);
    
    
    # parse options
    #
    #
    
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
    if (missing(organism)) {
        organism = "automatic";
    }
    
    # set default quality
    if (missing(quality)) {
        quality = "automatic";
    }
    
    # check refdir exists
    #
    #
    
    if (!file.exists(refdir)) {
        
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                refdir, " Not Found." );
        
        stop();
    }
    
    refdir = normalizePath(refdir);
    
    
    
    # memory usage
    print.memory.usage("Master Memory Usage: Start")
    
    
    # create projects
    #
    #
    log.info( "creating projects" );
    
    resetProjectErrors();
        
    # create projects
    projects <- createFastQProjects( accession = accession, organism = organism, quality = quality, 
            dir = dir, refdir=refdir, options = options )
    
    checkProjectErrors();
    
    # register step status
    #
    registerExpStepCompleted(PARAMETER_SANITY);

    
    # start projects
    #
    #
    
    if (length(projects) > 0) {
    
        attr(projects, 'metadata')$logfolder = makeClusterLogFolder(accession, dir, getSubmID());
        
        assignProjects(projects);
        
        # optimize node usage
        if (rcloudoptions$nnodes == "automatic") {
            rcloudoptions$nnodes = length(.projects);
        }
        
        ecount = runProjects( .projects, usercloud = usercloud, rcloudoptions = rcloudoptions, 
                want.reports = want.reports, filter = options$filter, steplist = steplist);
        
        # update status to OK
        #
        setExpStatus(COMPUTED, "Completed");

    } else {
        log.info("No data found to compute.");
        ecount = NULL;

        # update status to OK
        #
        setExpStatus(NODATA, "No data found to compute.");

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
        
        # register step status
        #
        registerExpStepFailure(OBTAINING_RAWDATA, 
                "One or more Fastq files are missing!");
        
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

assignSubmID <- function(submid){
    assignInNamespace('.submid', submid, ns="ArrayExpressHTS");
}

getSubmID <- function(){
    get('.submid');
}

runProjects <- function( .projects, usercloud = TRUE, rcloudoptions = rcloudoptions, 
        want.reports = FALSE, filter = TRUE, steplist ) {
    
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
        
        
        # Raw Reports
        #
        #
        #
        if (want.reports) {
            
            log.info( "generating Raw Reports" );
            
            registerExpStepStarted(CLUSTER_RAW_REPORT);
            
            result = try( clusterApplyLB(getCluster(), .projects, ArrayExpressHTS:::processOneProjectRawReport, 
                            .projects, want.reports, .stop.on.warnings) ); 
            
            if (inherits(result, "try-error")) { 
                #
                #
                registerExpStepFailure(CLUSTER_RAW_REPORT, " Error during the Raw Report");
                
            } else {
                #
                #
                registerExpStepCompleted(CLUSTER_RAW_REPORT);
            }
        }
        
        if ("align" %in% steplist || "alignment" %in% steplist) {
            
            # Alignment
            #
            #
            #
            log.info( "producing Alignment" );
            
            registerExpStepStarted(CLUSTER_ALIGNMENT);
            
            result = try( clusterApplyLB(getCluster(), .projects, ArrayExpressHTS:::processOneProjectAlignment, 
                            .projects, want.reports, .stop.on.warnings) ); 
            
            if (inherits(result, "try-error")) { 
                #
                #
                registerExpStepFailure(CLUSTER_ALIGNMENT, " Error during the Alignment");
                
            } else {
                #
                #
                registerExpStepCompleted(CLUSTER_ALIGNMENT);
            }
            
        }
        
        if ("count" %in% steplist || "estimation" %in% steplist) {
            
            # Estimation
            #
            #
            #
            log.info( "producing Estimation" );
            
            registerExpStepStarted(CLUSTER_ESTIMATION);
            
            result = try( clusterApplyLB(getCluster(), .projects, ArrayExpressHTS:::processOneProjectEstimation, 
                            .projects, want.reports, .stop.on.warnings) ); 
            
            if (inherits(result, "try-error")) { 
                #
                #
                registerExpStepFailure(CLUSTER_ESTIMATION, " Error during the Estimation");
                
            } else {
                #
                #
                registerExpStepCompleted(CLUSTER_ESTIMATION);
            }
            
        }
        
        
        # Alignment Reports
        #
        #
        #
        if (want.reports) {
            
            log.info( "generating Alignment Reports" );
            
            registerExpStepStarted(CLUSTER_ALN_REPORT);
            
            result = try( clusterApplyLB(getCluster(), .projects, ArrayExpressHTS:::processOneProjectAlnReport, 
                            .projects, want.reports, .stop.on.warnings) ); 
            
            if (inherits(result, "try-error")) { 
                #
                #
                registerExpStepFailure(CLUSTER_ALN_REPORT, " Error during the Alignment Report");
                
            } else {
                #
                #
                registerExpStepCompleted(CLUSTER_ALN_REPORT);
            }
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
    
    ecount = NULL;
    
    if ("eset" %in% steplist || "robject" %in% steplist || "object" %in% steplist) {
        
        registerExpStepStarted(ASSEMBLING_ESET);
        
        
        log.info( "building expression set" );
        
        ecount = to_expressionset( filter = filter )

        registerExpStepCompleted(ASSEMBLING_ESET);

    }

    
    # memory usage
    print.memory.usage("Memory: runProjects step 4")
    
    log.info( "plotting compared report" );
    
    if (want.reports) {
        if( length(.projects) > 1 ) {
            
            # register REPORT step
            registerExpStepStarted(COMPARED_REPORT);
            
            plot_compared_raw_report(.projects[[1]])
        
            # step completed
            registerExpStepCompleted(COMPARED_REPORT);

        }  else {
            log.info( "Compared report omitted for project with only one sample" );
        }
    }
        
    registerExpStepCompleted(EXP_STATUS);
    
    # memory usage
    print.memory.usage("Memory: runProjects step 5")
        
    return(ecount);
}


setProjectData <- function(unused, projects) {
    trace.enter("setProjectData");
    on.exit({ trace.exit() })
    
    assignProjects(projects)
}


processOneProjectRawReport <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProjectRawReport");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects);
    assignOneProject(project);
    
    # store ids for PSR reporting
    #
    pokePSRRunId(.project$name);
    pokePSRFolder(.project$psrdir);
    
    # set the RUN status
    #
    setRunStatus(STARTED, " ");
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name, " Raw Report" );
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 1")
    
    #.project$qual_type = class(quality( seq [[1]] )[1])[1]
    
    if (want.reports) {
        
        registerRunStepStarted(RAW_REPORT);
        
        seq = try( fastq_to_shortreadq() )
        if (inherits(seq, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error reading data file");
            stop();
            
            #return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 2")
        
        tab = try( shortread_to_tab( seq ) );
        if (inherits(tab, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error converting data file");
            stop();
            
            #return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 3")
        
        # make sure there's at least one device there
        #
        plottingresult = try( plot_raw_report( seq, tab ) )
        if (inherits(plottingresult, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error plotting raw report");
            stop();

            # continue even if it fails
        }
        
        registerRunStepCompleted(RAW_REPORT);
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 4")
    }
}
    
processOneProjectAlignment <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProjectAlignment");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects);
    assignOneProject(project);
    
    # store ids for PSR reporting
    #
    pokePSRRunId(.project$name);
    pokePSRFolder(.project$psrdir);

    # set the RUN status
    #
    setRunStatus(STARTED, " ");
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name, " Alignment" );
    
    
    
    registerRunStepStarted(ALIGNMENT);
    
    alignmentresult = try( align( sure.i.want.to.run.this = TRUE ) )
    if (inherits(alignmentresult, "try-error")) { 
        registerRunStepFailure(ALIGNMENT, " Error aligning data");
        stop();

        #return();
    }
    
    registerRunStepCompleted(ALIGNMENT);
    
    setRunStatus(OK, " ");

    
    # memory usage
    print.memory.usage("Memory: processOneProject step 5")
}    
    
processOneProjectEstimation <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProjectEstimation");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects);
    assignOneProject(project);
    
    # store ids for PSR reporting
    #
    pokePSRRunId(.project$name);
    pokePSRFolder(.project$psrdir);
    
    # set the RUN status
    #
    setRunStatus(STARTED, " ");
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name, " Estimamtion" );
    
    registerRunStepStarted(ESTIMATION);
    
    #estimationresult = try( call_cufflinks( run = TRUE ) )
    estimationresult = try( call_estimator( run = TRUE ) )
    if (inherits(estimationresult, "try-error")) { 
        
        registerRunStepFailure(ESTIMATION, " Error estimating expression");
        stop();

        #return();
    }
    
    registerRunStepCompleted(ESTIMATION);
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 6")

    # call it here to pre-cache 
    # when we're still in parallel
    total_reads = try( number_raw_reads( ) )
    if (inherits(total_reads, "try-error")) { 
        log.info(" Error getting number of reads");
        # continue
    }
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 7")
    
    # call it here to pre-cache 
    aligned_reads = try( number_aligned_reads() )
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 8")

    # set the RUN status
    #
    setRunStatus(OK, " ");
    
}    
    
processOneProjectAlnReport <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProjectAlnReport");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects);
    assignOneProject(project);
    
    # store ids for PSR reporting
    #
    pokePSRRunId(.project$name);
    pokePSRFolder(.project$psrdir);
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name, " Alignment Report" );
    
    if (want.reports) {
        
        registerRunStepStarted(ALN_REPORT);
        
        plottingresult = try( plot_aligned_report( ) )
        
        if (inherits(plottingresult, "try-error")) { 
            
            registerRunStepFailure(ALN_REPORT, " Error plotting alignment report");
            stop();
        }
        
        registerRunStepCompleted(ALN_REPORT);
        
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 9")
        
    }
    
    # run complete
    #
    #registerRunStepCompleted(RUN_STATUS);
    
    # set the RUN status
    #
    setRunStatus(OK, " ");

    
    gc();
    log.info( "done" );
}

processOneProject <- function(project, projects, want.reports, stop.on.warnings) {
    trace.enter("processOneProject");
    on.exit({ trace.exit() })
    
    on.exit({ assignStopOnWarnings(FALSE); cleanupAfterExecution(); })
    
    assignStopOnWarnings(stop.on.warnings);
    assignProjects(projects)
    assignOneProject(project)
    
    # store ids for PSR reporting
    #
    pokePSRRunId(.project$name);
    pokePSRFolder(.project$psrdir);
    
    # set current directory to project folder 
    # so that nothing written in parallel from
    # different projects collides/clashes
    setwd(.project$projectdir);
    
    log.info( "processing ", .project$name );
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 1")
    
    #.project$qual_type = class(quality( seq [[1]] )[1])[1]
    
    if (want.reports) {
        
        registerRunStepStarted(RAW_REPORT);
        
        seq = try( fastq_to_shortreadq() )
        if (inherits(seq, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error reading data file");
            stop();
            
            #return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 2")
        
        tab = try( shortread_to_tab( seq ) );
        if (inherits(tab, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error converting data file");
            stop();

            #return();
        }
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 3")
        
        # make sure there's at least one device there
        #
        plottingresult = try( plot_raw_report( seq, tab ) )
        if (inherits(plottingresult, "try-error")) { 
            registerRunStepFailure(RAW_REPORT, " Error plotting raw report");
            stop();
            
            # continue even if it fails
        }

        #rm(bam);
        rm(seq);

            
        registerRunStepCompleted(RAW_REPORT);
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 4")
    }
    
    registerRunStepStarted(ALIGNMENT);
    
    alignmentresult = try( align( sure.i.want.to.run.this = TRUE ) )
    if (inherits(alignmentresult, "try-error")) { 
        registerRunStepFailure(ALIGNMENT, " Error aligning data");
        stop();
        
        #return();
    }
    
    registerRunStepCompleted(ALIGNMENT);
    
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 5")
    
    registerRunStepStarted(ESTIMATION);
    
    #estimationresult = try( call_cufflinks( run = TRUE ) )
    estimationresult = try( call_estimator( run = TRUE ) )
    if (inherits(estimationresult, "try-error")) { 
        
        registerRunStepFailure(ESTIMATION, " Error estimating expression");
        stop();

        #return();
    }
    
    registerRunStepCompleted(ESTIMATION);
    
    # memory usage
    print.memory.usage("Memory: processOneProject step 6")
    
    
    #bam = try( bam_to_list( return = FALSE ) )
    #if (inherits(bam, "try-error")) { 
    #    log.info(" Error converting bam to list");
        # continue
    #}
    
    if (want.reports) {
        
        registerRunStepStarted(ALN_REPORT);
        
        plottingresult = try( plot_aligned_report( ) )
        
        if (inherits(plottingresult, "try-error")) { 
            
            registerRunStepFailure(ALN_REPORT, " Error plotting alignment report");
            stop();

            # continue
        }
        
        registerRunStepCompleted(ALN_REPORT);
        
        
        # memory usage
        print.memory.usage("Memory: processOneProject step 7")
        
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
    
    # run complete
    #
    # registerRunStepCompleted(RUN_STATUS);
    
    #
    #
    setRunStatus(OK, " ")
        
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

        
setupBaseFolder = function(accession, dir, submid) {
    basedir = "";
    
    if (!is.null( submid )) {
        #
        #
        basedir = paste( dir, accession, submid, sep="/");
    } else {
        #
        #
        basedir = paste( dir, accession, sep="/");
    }
    
    if (!file.exists( basedir )) {
        dir.create( basedir, recursive = TRUE, showWarnings = FALSE );
    }
    
    return( basedir );
}

setupDataFolder = function( accession, dir ) {
    
    datadir = paste( dir, accession, "data", sep="/" );
    
    if (!file.exists( datadir )) {
        dir.create( datadir, recursive = TRUE, showWarnings = FALSE );
    }
    
    return( datadir );
}


