#
# Logs
#
# moved to a separate file
#

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
    
    #stop();
}


getMark = function() {
    " [AEHTS] ";
}
