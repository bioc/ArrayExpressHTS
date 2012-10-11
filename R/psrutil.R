#
# PSR (Programmatic Status Report) related Routines
#
#

# package scope
#
.expid = "";
.runid = "";
.psrfolder = "";

#
# EXP
#
PARAMETER_SANITY    = "PARAMETER_SANITY"
OBTAINING_METADATA  = "OBTAINING_METADATA"
METADATA_SANITY     = "METADATA_SANITY"
REFERENCE_SANITY    = "REFERENCE_SANITY"
OBTAINING_RAWDATA   = "OBTAINING_RAWDATA"
ALLOCATING_CLUSTER  = "ALLOCATING_CLUSTER"

CLUSTER_RAW_REPORT  = "CLUSTER_RAW_REPORT"
CLUSTER_ALIGNMENT   = "CLUSTER_ALIGNMENT"
CLUSTER_ESTIMATION  = "CLUSTER_ESTIMATION"
CLUSTER_ALN_REPORT  = "CLUSTER_ALN_REPORT"

ASSEMBLING_ESET  = "ASSEMBLING_ESET"
COMPARED_REPORT  = "COMPARED_REPORT"
EXP_STATUS  = "EXP_STATUS"


#
# RUN
#
RAW_REPORT  = "RAW_REPORT"
ALIGNMENT   = "ALIGNMENT"
ESTIMATION  = "ESTIMATION"
ALN_REPORT  = "ALN_REPORT"
RUN_STATUS  = "RUN_STATUS"


# step status
#
STARTED   = "STARTED"
FAILED    = "FAILED"
WARNING   = "WARNING"
OK        = "OK"

# exp status
#
COMPUTED  = "COMPUTED"
NODATA    = "NODATA"
NOSERVICE = "NOSERVICE"

.exp.status = FAILED;
.run.status = FAILED;


#
# functions
#

pokePSRFolder = function(psrfolder) {
    assignInNamespace('.psrfolder', psrfolder, ns="ArrayExpressHTS");
}

peekPSRFolder = function() {
    return(get('.psrfolder'));
}

pokePSRRunId = function(runid) {
    assignInNamespace('.runid', runid, ns="ArrayExpressHTS");
}

peekPSRRunId = function() {
    return(get('.runid'));
}

pokePSRExpId = function(expid) {
    assignInNamespace('.expid', expid, ns="ArrayExpressHTS");
}

peekPSRExpId = function() {
    return(get('.expid'));
}

#
# EXP status
#

# setter, additionally writes on PSR
setExpStatus = function(expstatus, text = " ") {
    assignInNamespace('.exp.status', expstatus, ns="ArrayExpressHTS");

    # write it on PSR
    registerExpStepStatusText(EXP_STATUS, expstatus, text);
}

# getter
getExpStatus = function() {
    return(get('.exp.status'));
}

# update on exit
updateExpStatusOnExit = function() {
    
    if (getExpStatus() == STARTED) {
        # if status is STARTED (means we exit prematurely
                # on 'error' or 'stop') set status to FAILED
        # 
        registerExpStepStatusText(EXP_STATUS, FAILED, "Unspecified Failure")
    }
}


# COMPLETED
#
registerExpStepCompleted = function(step) {
    registerExpStepStatusText(step, OK, " ");
}

# FAILED
#
registerExpStepFailure = function(step, ...) {
    message = paste(date(), getMark(), ... , sep="");
    
    registerExpStepStatusText(step, FAILED, message);
    
    expstatus = paste(step, "_NOT_COMPLETED", sep = "");
    
    setExpStatus(expstatus, message);
    
    log.error(...);
}

# WARNING
#
registerExpStepWarning = function(step, ...) {
    registerExpStepStatusText(step, WARNING, paste(date(), getMark(), ... , sep=""));
    log.warning(...);
}

# STARTED
#
registerExpStepStarted = function(step) {
    registerExpStepStatusText(step, STARTED, " ");
}


registerExpStepStatusText = function(step, status, text = " ") {
    
    expid = peekPSRExpId();
    folder = peekPSRFolder();
    
    if (nchar(expid) > 0 && nchar(folder) > 0) {
        # write the status
        #
        fname = paste(folder, "/", expid, ".EXP.STATUS", sep="");
        appendStatusText(fname, step, status, text);
    }
    
}


#
# RUN status
#

setRunStatus = function(runstatus, text = " ") {
    # update status in memory
    assignInNamespace('.run.status', runstatus, ns="ArrayExpressHTS");

    # write it on PSR
    registerRunStepStatusText(RUN_STATUS, runstatus, text);
}

getRunStatus = function() {
    return(get('.run.status'));
}


updateRunStatusOnExit = function() {
    
    if (getRunStatus() == STARTED) {
        # if status is STARTED
        # set status to FAILED
        #
        registerRunStepStatusText(RUN_STATUS, FAILED, "Unspecified Failure");
    }
    
}

registerRunStepCompleted = function(step) {
    registerRunStepStatusText(step, OK, " ");
}

registerRunStepFailure = function(step, ...) {
    registerRunStepStatusText(step, FAILED, paste(date(), getMark(), ... , sep=""));
    log.error(...);
}

registerRunStepWarning = function(step, ...) {
    registerRunStepStatusText(step, WARNING, paste(date(), getMark(), ... , sep=""));
    log.warning(...);
}

registerRunStepStarted = function(step) {
    registerRunStepStatusText(step, STARTED, " ");
}


# low level 'writer'
#
registerRunStepStatusText = function(step, status, text = " ") {
    runid = peekPSRRunId();
    folder = peekPSRFolder();
    
    if (nchar(runid) > 0 && nchar(folder) > 0) {
        # write the status
        #
        fname = paste(folder, "/", runid, ".RUN.STATUS", sep="");
        appendStatusText(fname, step, status, text);
    }
    
}


# low level 'writer'
#
appendStatusText = function(filename, step, status, text) {
    df0 = data.frame("step" = c(step), "status" = c(status), "text"=c(text));
    write.table(df0, file = filename, append = TRUE,
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE);
}


