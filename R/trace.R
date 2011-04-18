getTraceMark = function() {
    "[TRACE] ";
}

trace.enter = function(name, ...) {
    if (getPipelineOption("trace") == "enabled") {
        step = getPackageVariable("trace.step");
        level = getPackageVariable("trace.level");
        level = level + 1;
        
        setPackageVariable("trace.level" = level);
        
        cat(paste(date(), getMark(), getTraceMark(), 
            paste(rep(step, level),collapse=""), ">", " ", name, ... , "\n", sep=""));
    }
}

trace.exit = function(name, ...) {
    if (getPipelineOption("trace") == "enabled") {
        level = getPackageVariable("trace.level");
        if (level > 0) {
            setPackageVariable("trace.level" = (level - 1));
        }
    }
}

