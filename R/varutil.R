#
#  Package Variables
#
#
# package veriages are a useful thing to have if 
# you want to cache something on a global level and 
# access whenever needed, plus it provides an easy
# trace means when you need to check all places a 
# variable is read or written to, plus it doesn't 
# overlap with local vars, minus it makes the code a 
# little ugly and hard to read.
#
# access to global environment variables is a very 
# bad thing, completely untrackable, prone to errors 
# and should preferably be avioded.
#
#

defaultVariables <- NULL

initPackageVariables <- function() {
    variables <- list(
        trace.level = 0,
        trace.step = "--"
    )
    
    defaultVariables <<- addPackageVariables(emptyenv(), variables)
}

addPackageVariables <- function(variables, new) {
    if (! is.null(new)) {
        variables <- new.env(parent = variables)
        names <- names(new)
        for (i in seq(along = new))
        assign(names[i], new[[i]], envir = variables)
    }
    variables
}

getPackageVariable <- function(name, variables = defaultVariables) {
    
    if (exists(name, envir=variables)) {
        return(get(name, envir = variables));
    } else {
        return(NULL)
    }
}

isPackageVariableDefined <- function(name, variables = defaultVariables) {
    return(exists(name, envir=variables))
}

setPackageVariable <- function(...) {
    list <- list(...)
    names <- names(list)
    for (i in seq(along = list))
    assign(names[i], list[[i]], envir = defaultVariables)
}

