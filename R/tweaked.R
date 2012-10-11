myHTMLplot = function (Caption = "", file = get(".HTML.file"), append = TRUE, 
        GraphDirectory = ".", GraphFileName = "", GraphSaveAs = "png", 
        GraphBorder = 1, Align = "center", Width = 500, Height = 500, 
        WidthHTML = NULL, HeightHTML = NULL, GraphPointSize = 12, 
        GraphBackGround = "white", GraphRes = 72, plotFunction = NULL, ...) 
{
    if (exists("HTMLenv", where = ".GlobalEnv", mode = "environment")) {
        GraphDirectory = get(".HTML.outdir", envir = get("HTMLenv", 
        envir = .GlobalEnv))
    }
    cat("\n", file = file, append = append, ...)
    if (GraphFileName == "") {
        nowd <- date()
        GraphFileName <- paste("GRAPH_", substring(nowd, 5, 7), 
            substring(nowd, 9, 10), "_", substring(nowd, 12, 
            13), substring(nowd, 15, 16), substring(nowd, 
            18, 19), sep = "")
    }
    GraphFileName <- paste(GraphFileName, ".", GraphSaveAs, sep = "")
    AbsGraphFileName <- file.path(GraphDirectory, GraphFileName)
    
    if (GraphSaveAs == "png") {
        if (is.null(plotFunction)) 
            dev.print(device = png, file = AbsGraphFileName, 
                width = Width, height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
        else {
            png(filename = AbsGraphFileName, width = Width, 
                height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
            plotFunction()
            dev.off()
        }
    } else {
        log.error(paste(GraphSaveAs, "is not supported"));
        stop();
    }
    
    cat(paste("<p align=", Align, "><img src='", GraphFileName, 
            "' border=", GraphBorder, 
            if (!is.null(Width)) paste(" width=", Width, sep = "")
            else "", 
            if (!is.null(HeightHTML)) paste(" height=", HeightHTML, sep = ""), 
            if (!is.null(WidthHTML)) paste(" width="), ">", sep = "", collapse = ""), 
            file = file, append = TRUE, sep = "")
    
    if (Caption != "") {
        cat(paste("<br><font class=caption>", Caption, "</font>"), 
        file = file, append = TRUE, sep = "")
    }
    cat("</p>", file = file, append = TRUE, sep = "\n")
    if (substitute(file) == ".HTML.file") 
        try(assign(".HTML.graph", TRUE, envir = get("HTMLenv", envir = .GlobalEnv)))
    invisible(return(TRUE))
}

