#
# Everything Cluster Related
#

# cluster variable
#
cluster001 <- NULL;

# cluster getter
#
getCluster <- function() {
    return(cluster001);
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
        
        # register step status
        #
        registerExpStepFailure(ALLOCATING_CLUSTER,
                "RCLOUD cluster is not available. Please run from R-Cloud Workbench.");
        
        stop();
    }
    
    # register step status
    #
    registerExpStepStarted(ALLOCATING_CLUSTER);
    
    
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
        cl0 = try({ makeCluster(rcloudoptions$nnodes - length(finalcluster), 
                            type="RCLOUD", pool = rcloudoptions$pool) });
        
        if (inherits(cl0, "try-error")) {
            
            # cannot create computing cluster
            #
            registerExpStepFailure(ALLOCATING_CLUSTER,
                    "Failed to create RCLOUD computing cluster");
                
            stop();
        }
        
        finalcluster = mergeClusters(cl0, finalcluster);
        
        if (length(finalcluster) == rcloudoptions$nnodes) {
            # cluster created successfully
            #
            break;
        } else {
            log.info( attempt, " attempt to create a cluster" );
            
            if (attempt == rcloudoptions$nretries) {
                
                if (length(finalcluster) == 0) {
                    
                    # register step status
                    #
                    registerExpStepFailure(ALLOCATING_CLUSTER,
                            "Failed to create RCLOUD computing cluster");
                    
                    stop();
                    
                } else {
                    
                    registerExpStepWarning(ALLOCATING_CLUSTER, 
                            "Only ", length(finalcluster), " nodes of ", 
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
    createDirAndBackup(logfolder, showWarnings = FALSE);
    
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
    
    # register step status
    #
    registerExpStepCompleted(ALLOCATING_CLUSTER);
    
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
    on.exit({ trace.exit() });
    
    hostname = Sys.info()[4];
    
    jobid = try({ do.call(".jcall", list( obj="uk/ac/ebi/rcloud/server/RListener", 
                                "[Ljava/lang/String;" ,"getProperty", as.character("job.id"))) });
            
    b = "]";
    
    if (!inherits(jobid, "try-error")) {
        # jobid found
        #
        logname = paste(name, jobid, hostname, sep="-")
        
        storage = "/ebi/microarray/home/biocep/service/VirtualRWorkbench/";
        
        loglog = paste("rcloud-server.", jobid, ".log", sep="");
        lsferr = paste(jobid, ".err", sep="");
        lsfout = paste(jobid, ".out", sep="");
        
        system(paste("ln -s ", storage, "/logs/", loglog, " ", logfolder, "/", loglog, sep="" ));
        system(paste("ln -s ", storage, "/lsflogs/", lsferr, " ", logfolder, "/", lsferr, sep="" ));
        system(paste("ln -s ", storage, "/lsflogs/", lsfout, " ", logfolder, "/", lsfout, sep="" ));
        
    } else {
        logname = paste(name, hostname, sep="-")
    }
    
    path = paste(logfolder, logname, sep='/'); 
    sinkWorkerOutput(path); 
    options('outfile0' = path) 

    
    

    
            
}
        
        