plot_compared_hits <- function( update ) {
    trace.enter("plot_compared_hits");
    on.exit({ trace.exit() })
    
    cols=c("darkgoldenrod","darkgoldenrod4","darkgoldenrod1")
    par( mfrow=c(1,1), xpd = FALSE )
    
    data = NULL
    for( i in names(.projects) ) {
        log.info( "At project ", i, ":" )
        
        assignOneProject(.projects[[i]])
        
        bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )
        total_reads = number_raw_reads( update=update )
        
        if( .project$pairing$type == "PE" ) {
            flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                    isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                    isFirstMateRead = TRUE )
        } else {
            flags = scanBamFlag( isUnmappedQuery = FALSE )
        }
        filename = paste( .project$aligner$out_dir, "/hits_per_read.RData", sep="" )
        if( !file.exists(filename) || update ) {
            log.info( "\treading BAM" )
            qnames = scanBam( bam_filename, param=ScanBamParam( flag=flags, 
                        what=("qname") ) )[[1]]$qname
            log.info( "\tcounting reads" )
            tab = table( qnames )
            rm(qnames)
            gc()
            save(tab, file=filename)
        } else {
            log.info( "\tfile exists. loading..." )
            log.info( "\tLoading..." )
            load( filename )
        }
        
        total_aligned = length(tab)
        unique_match = sum( tab == 1 )
        multiple_match = total_aligned-unique_match
        rm(tab)
        gc()
        
        if( !is.null(data) ) 
            data = cbind( data, 
                    matrix( c(unique_match, multiple_match, total_reads-(unique_match+multiple_match)), ncol=1 ) )
        else 
            data = matrix(  c(unique_match, multiple_match, total_reads-(unique_match+multiple_match)), ncol=1 )
        colnames(data)[ncol(data)] = .project$name
    }
    log.info( "Done with the projects, plotting" )
    xc = barplot( data, beside = FALSE, col=cols, ylab="Number of reads", las=2 )
    par( xpd=NA )
    for( i in seq(1,length(.projects)) ) {
        legend( xc[i], data[1,i]/2, paste("~",round(data[1,i]*100/sum(data[,i])), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )
        legend( xc[i], data[1,i]+data[2,i]/2, paste("~",round(data[2,i]*100/sum(data[,i])), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )
        legend( xc[i], data[1,i]+data[2,i]+data[3,i]/2, paste("~", round(data[3,i]/sum(data[,i])*100), "%", sep=""), 
                bty="n", xjust=0.5, yjust=0.5 )
        #legend( xc[i+1], data[1,i+1]/2, paste("unique\nhit\n~",round(data[1,i+1]*100/data[1,i]), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )
        #legend( xc[i+1], data[1,i+1]+(data[2,i+1]/2), paste("mult.\nhits\n~",round(data[2,i+1]*100/data[1,i]), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )
    }
    legend( "topright", legend=c( "unique hit", "mult. hit", "not aligned"), fill=cols )
    log.info( "Done!" )
}

plot_scatter_samples <- function() {
    trace.enter("plot_scatter_samples");
    on.exit({ trace.exit() })
    
    e = to_expressionset( )
    
    cb = combn( sampleNames(assayData(e)), 2 )
    e10 = log10( exprs(e) )
    if( ncol(cb) >= 6 ) {
        log.info( "*** Please note, this function can only plot 6x6 samples, any extra will not be plotted. ***" )
    }
    ncols = min( 6, ncol(cb) )
    par( mfrow=c(ceiling(ncols/2),ceiling(ncols/2)) )
    max = ceiling(max(e10))
    
    for( i in 1:ncols ) {
        plot( e10[,cb[1,i]], e10[,cb[2,i]], 
                xlab=paste( cb[1,i], "log10 FPKM"), 
                ylab=paste( cb[2,i], "log10 FPKM"), 
                ylim=c(0, max), xlim=c(0, max) )
    } 
}

# a plot like the one in the Oshlack paper to investigate whether DE is biased by transcript length
plot_DE_length <- function( e, de, pval=c(0.01,0.05) ) {
    trace.enter("plot_DE_length");
    on.exit({ trace.exit() })
    
    de = de[!is.na(de$pval),]
    par( mfrow=c(1,2) )
    for( p in pval ) {
        det = de$id[de$padj <= p]
        h1 = hist( log10(fData(e)[,]), breaks=250, plot = FALSE )
        l2 = fData(e)[featureNames(e) %in% det,]
        #h2 = sapply( h1$breaks[-1], function(x) sum(l2 < x) )
        h2 = hist( log10(l2), breaks=h1$breaks, plot = FALSE )
        filt = h2$counts > 0
        plot( h1$breaks[-1][filt], h2$counts[filt]/h1$counts[filt], 
                main=paste( "pval", p ), xlab="Transcript length", 
                ylab="Percentage of transcripts reported DE" )
    }
}

plot_DE_report <- function( e, de, update = TRUE ) {
    trace.enter("plot_DE_report");
    on.exit({ trace.exit() })
    
    report_dir = .projects[[1]]$report["compared_report_dir"]
    
    if( !file.exists(get(".HTML.file")) || update ) {
        HTML( as.title( "DE report" ), append = FALSE )
        plot_DE_length( e, de )
        myHTMLplot( GraphDirectory=report_dir, Width=2*500 )
        HTMLEndFile()
    } else {
        log.info( "File exists, not updated" )
    }
}

plot_insize_boxplot <- function() {
    trace.enter("plot_insize_boxplot");
    on.exit({ trace.exit() })
    
    col="darkgoldenrod3"
    
    isizes = list()
    for( i in names(.projects) ) {
    
        assignOneProject(.projects[[i]])
        
        bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
        isizes[[.project$name]] = scanBam( bam_filename, param=ScanBamParam( flag=flags, what=("isize") ) )[[1]]$isize
    }
    
    boxplot( isizes, outline = FALSE, xlab="Insert size", main="Distribution of insert sizes", col=col )
}

plot_avgbasequal_boxplots <- function( update = FALSE ) {
    trace.enter("plot_avgbasequal_boxplots");
    on.exit({ trace.exit() })
    
    par( mfrow=c(1,1) )
    
    col="darkgoldenrod3"
    quals = list()
    for( i in 1:length(.projects) ) {
        #progress( i, max.value=length(.projects) )
        
        assignOneProject(.projects[[i]])
        
        log.info( "Plotting avgbasequals for project ", .project$name )
        
        filename = .project$aux_files["avgQuals"]
        
        if( !file.exists( filename ) ) {
            log.info( "File ", filename, " does not exist. Calculating avgquals..." )
            seq = fastq_to_shortreadq()
            avgquals = phred_to_avgqual( seq, update=update )
        } else {
            load( filename );
        }
        if (length(avgquals) == 1) {
            quals[[.project$name]] = avgquals[[1]]
        } else {
            for (n in 1:length(avgquals)) {
                quals[[ paste(.project$name, "_", n, sep="") ]] = avgquals[[n]];
            }
        }
        gc()
    }
    par( las=2 )
    
    boxplot( quals, outline = FALSE, ylab="Average read quality", col=col )
}

# lanes are project indexes 
plot_dustyScore_distribs <- function( update ) {
    trace.enter("plot_dustyScore_distribs");
    on.exit({ trace.exit() })
    
    par( mfrow=c(1,2) )
    
    dusts = list()
    
    #cols = paste( brewer.pal(max(length(.projects),3), "Dark2"), "B4", sep="" )
    #cols = brewer.pal(max(length(.projects),3), "Dark2")
    cols = brewer.pal(max(min(length(.projects),12),3), "Set3")
    lanes = seq(1:min(length(.projects),12))
    
    for( i in lanes ) {
        #progress( i, max.value=length(lanes) )
        ##j = lanes[i]
        
        assignOneProject(.projects[[i]])
        
        filename = .project$aux_files["dusty"]
        
        if( !file.exists( filename ) || update ) {
            log.info( "File ", filename, " does not exist!" );
            alldata = fastq_to_shortreadq()
            dustS = calculate_dustyScore( alldata, update )
        } else {
            load( filename )
        }
        dusts[[.project$name]] = dustS
        if( i == 1 ) {
            plot( ecdf(dusts[[i]]), main="Read complexity", xlab="Dusty score", ylab="Proportion of reads",
                    do.points = FALSE, verticals = TRUE, lwd=3, col=cols[i] )
        } else {
            lines( ecdf(dusts[[i]]), do.points = FALSE, verticals = TRUE, lwd=3, col=cols[i] )
        }
        dusts[[.project$name]] = density(dusts[[i]])
        
    }
    
    for( i in lanes ) {
        if( i == 1 ) {
            plot( dusts[[i]], main="Read complexity", xlab="Dusty score", 
                    ylim=c(0,max(sapply(dusts, function(x) x$y))), col=cols[i], lwd=2 )
        } else {
            lines( dusts[[i]], col=cols[i], lwd=2 )
        }
    }
    legend( x="topright", sapply(.projects, function(x) x$name)[lanes], fill=cols[lanes] )
    
    par( mfrow=c(1,1) )
}

