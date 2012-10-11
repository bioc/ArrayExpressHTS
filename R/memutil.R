#
# Memory Monitor
#
# This routine prints memory usage statistics. It calls garbage 
# collector, therefore must be used with caution, should not be 
# called from a time critical, processing intensive regions like 
# processing cycles, etc.
#
print.memory.usage = function( logtext ) {
    
    if (getPipelineOption("memorymonitor") == "enabled") {
    
        # print text log
        log.info(logtext);
        
        # collect the stats
        memusage = gc();
        
        # print the stats
        print(memusage);
    
    }

}
