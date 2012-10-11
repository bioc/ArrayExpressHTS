#
# On-load/unload Routines
#

.onLoad <- function(libname, pkgname = "ArrayExpressHTS") {
    #library.dynam("ArrayExpressHTS", pkgname, libname)

    # first: init options and defaults
    #
    if (is.null(defaultOptions)) {
        initPipelineOptions();
    }

    if (is.null(defaultVariables)) {
        initPackageVariables();
    }
    
    # second: init default environment
    #
    initEnvironmentVariables();

}

.onUnload <- function(libpath) {
    #library.dynam.unload("ArrayExpressHTS", libpath = libpath)
}
