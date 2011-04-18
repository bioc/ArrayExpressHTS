# Plotting functions called by plot_report() on tools.R
# 
# Author: filimon
###############################################################################

plot_quality_per_cycle <- function( odata ) {
    trace.enter("plot_quality_per_cycle");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "qualityPerCycle.png", sep="" )
    
    col="darkgoldenrod3"
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(odata)*480 )
    #}
    
    par( mfrow=c(1,length(odata)) )
    
    i = 1
    for( data in odata ) {
        qfreq = alphabetByCycle( quality(data) ) 
        
        offset = 64
        if( .project$qual_type == "FastqQuality" )
            offset = 33
        scores = sapply( rownames(qfreq), function(x) ASCIIconvert( x, type="CtoI") - offset )
        
        cycles = 1:unique(width(sread(data)))
        avgquals = sapply(cycles, function(x) sum(qfreq[,x]*scores)/length(data) )
        
        
        #boxplot( qfreq*scores, ylim=c(-5,40) )
        stdevs = sqrt( sapply(cycles, function(x) wtd.var( scores, qfreq[,x] ) ) )
        errbar( cycles, avgquals, avgquals-stdevs, avgquals+stdevs, 
                xlab="Cycle", ylab="Average base quality", ylim=c(-5,45), col=col )
        #title( main="1st mate" )
        par( xpd=NA )
        if( length(odata) > 1 ) {
            if( i == 1)
                legend( 30, 46, "1st mates", horiz = TRUE, bty="n" )
            else
                legend( 30, 46, "2nd mates", horiz = TRUE, bty="n" )
        }
        i = i+1
    }
    
    par(mfrow=c(1,1))
    title( main="Average base quality per cycle with 1 std dev error bar" ) 
    
    
    #if( plottofile ){
    #    dev.off()
    #}
    
}

#plot_chrhit_bars <- function( aln, index ) {
#    trace.enter("plot_chrhit_bars");
#    on.exit({ trace.exit() })
#    
#    cols = brewer.pal(3,"Spectral")
#    heights = sort(table(chromosome(aln)) )
#    names(heights) = get_chr_representation( names(heights) )
#    chrinfo = sort( chromInfo_to_dataframe(), decreasing = TRUE )
#    genlen = sum( chrinfo )
#    totreads = sum( heights )
#    data = t(matrix( c(chrinfo/genlen, heights[match(names(heights),names(chrinfo))]/totreads), ncol=2))
#    margins = par()$mai
#    margins[4] = margins[1]
#    par( mai=margins )
#    barplot( data, ylim=c(0,max(data)), beside = TRUE, names.arg=names(chrinfo), yaxt="n", las=2, col=cols[c(1,3)] )
#    axis( 2, at=NULL,labels = TRUE, col.axis=cols[1], col=cols[1], lwd=3 )
#    axis( 4, at=NULL,labels = TRUE, col.axis=cols[3], col=cols[3], lwd=3 )
#    legend( x="topleft", c("Proportion of chromosome length relative to the genome", "Proportion of reads relative to the total"), fill=cols[c(1,3)] )
#}

# TODO: make it a percentage of the length of the chromosome
plot_chromosome_distrib <- function( bam_filename ) {
    trace.enter("plot_chromosome_distrib");
    on.exit({ trace.exit() })
    
    if( .project$pairing$type == "PE" )
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
    else
        flags = scanBamFlag( isUnmappedQuery = FALSE )
    
    rname = scanBam( bam_filename, param=ScanBamParam( flag=flags, what="rname") )[[1]]$rname
    tab = table( rname )
    
    barplot( tab[tab>0], las=2, col="darkgoldenrod3", main="Read count per chromosome" )
}

searchBiomart <- function( organism ) {
    
    av_marts = listMarts()[,1]
    n_marts = c( "ensembl", "metazoa", "bacterial", "fungal", "plant", "protist" )
    marts = as.character(av_marts[pmatch(n_marts, av_marts )])
    
    organism = strsplit( tolower(organism), "_" )[[1]]
    for( mart in marts ) {
        log.info( "Searching mart ", mart, " \n" )
        ensembl = useMart( mart )
        av_datasets = listDatasets(ensembl)[,1]
        dataset = as.character(
                av_datasets[regexpr(paste("*",organism[length(organism)],"*",sep=""),
                                av_datasets)>0])
        
        if( length( dataset ) != 0 ) {
            return( c(mart, dataset) )
        } 
    }
    return( NULL )
}

# for alignment to the transcriptome only
# TODO: not tested yet!
plot_transcript_distrib <- function( bam_filename ) {
    trace.enter("plot_transcript_distrib");
    on.exit({ trace.exit() })
    
    if( .project$pairing$type == "PE" )
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
    else
        flags = scanBamFlag( isUnmappedQuery = FALSE )
    
    bam = scanBam( bam_filename, param=ScanBamParam( flag=flags, what=c("pos","rname") ) )[[1]]
    
    mart_info = searchBiomart( .project$organism )
    if( !is.null(mart_info) ) {
        ensembl = useMart( mart_info[1], dataset=mart_info[2] )
        attrib = c( "ensembl_transcript_id", "transcript_start", "transcript_end" )
        chrs = getBM( attributes="chromosome_name", mart=ensembl )$chromosome_name
        annot = data.frame()
        for( chr in chrs ) {
            annot = rbind( annot, getBM( attributes=attrib, mart=ensembl, 
                            filters="chromosome_name", values=chr ) ) 
        }
        annot = cbind( annot, length=annot$transcript_end-annot$transcript_start )
    }
    
}

plot_strand_distrib <- function( bam_filename ) {
    trace.enter("plot_strand_distrib");
    on.exit({ trace.exit() })
    
    if( .project$pairing$type == "PE" )
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
    else
        flags = scanBamFlag( isUnmappedQuery = FALSE )
    
    strands = scanBam( bam_filename, param=ScanBamParam( flag=flags, what="strand") )[[1]]$strand
    tab = table( strands )
    
    pie( tab[tab>0], col=c("darkgoldenrod3", "darkgoldenrod4"), main="Read count per strand" )
}

plot_mapq <- function( bam_filename ) {
    trace.enter("plot_mapq");
    on.exit({ trace.exit() })
    
    col="darkgoldenrod3"

    if( .project$pairing$type == "PE" ) 
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
    else 
        flags = scanBamFlag( isUnmappedQuery = FALSE )
    
    mapq = scanBam( bam_filename, param=ScanBamParam( flag=flags, what=("mapq") ) )[[1]]$mapq
    
    par( mfrow=c(1,1) )
    hist(mapq, col=col, xlab="Mapping quality", main="Mapping quality" )
}

plot_insertsize <- function( bam_filename ) {
    trace.enter("plot_insertsize");
    on.exit({ trace.exit() })
    
    flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
            isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
            isFirstMateRead = TRUE )
    isizes = scanBam( bam_filename, param=ScanBamParam( flag=flags, what=("isize") ) )[[1]]$isize
    #isizes = isizes[!is.na(isizes)]
    isizes = abs(isizes)
    col="darkgoldenrod3"
    hist(isizes, breaks=c(-0.5:(max(isizes)+0.5)), main=paste("Insert size as estimated by", 
                    .project$aligner$type), col=col )
}

table_pileup_summary <- function( pile ) {
    trace.enter("table_pileup_summary");
    on.exit({ trace.exit() })
    
    bases = c("A","T","C","G")
    df = data.frame(row.names=bases)
    for( base in c(bases, "M", "R", "W", "S", "Y", "K") ) {
        df = cbind( df, rep(0,4) )
    }
    colnames( df ) = c(bases, "M", "R", "W", "S", "Y", "K")
    for( i in 1:nrow(pile) ) {
        #browser( expr=(i == 3775))
        if( pile[i,"refbase"] %in% rownames(df) )
            df[pile[i,"refbase"],pile[i,"consbase"]] = df[pile[i,"refbase"],pile[i,"consbase"]] + 1
    }
}

# receives a list of AlignedRead or ShortRead objects
# if plot = TRUE plots on a window, otherwise plots into a file
plot_basecall_per_cycle <- function( alldata ) {
    trace.enter("plot_basecall_per_cycle");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "basecallPerCycle.png", sep="" )
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(alldata)*480 )
    #}
    
    par( mfrow=c(1,length(alldata)) )
    
    cols = c( "red", "blue", "chartreuse4", "darkgoldenrod3", "black" )
    lwd = 2
    
    allfreq = list()
    j = 1
    ymax = 0
    for( data in alldata ) {
        data = sread(data)
        log.info( "Calculating alphabet by cycle..." )
        qfreq = alphabetByCycle( data )
        log.info( "Summing frequencies..." )
        treads = sum( qfreq[,1] )
        qfreq = qfreq[c(1,2,3,4,15),]
        qperc = round(qfreq/treads*100,2)
        
        log.info( "Plotting base call per cycle..." )
        ymax = max( ymax, max(qperc) )
        plot( 1:length(qperc[1,]), qperc[1,], xlab="Cycle", ylab="% of base count", 
                ylim=c( 0, ymax ), type="l", col=cols[1], lwd=lwd )
        for( i in 2:length(qperc[,1]) ) {
            lines( 1:length(qperc[1,]), qperc[i,], col=cols[i], lwd=lwd )
        }
        
        par( xpd=NA )
        if( length(alldata) > 1 ) {
            if( j == 1) {
                legend( ncol(qfreq)/2, ymax+round(ymax*0.04), "1st mates", horiz = TRUE, bty="n" )
                title( main=paste( "Base call per cycle" ) )
            } else
                legend( ncol(qfreq)/2, ymax+round(ymax*0.04), "2nd mates", horiz = TRUE, bty="n" )
        }
        
        allfreq[[j]] = qfreq
        j=j+1
    }
    par( xpd=NA )
    #top of the left hand side plot
    #legend( -38, 42, rownames(qperc), horiz = TRUE, fill=cols, bty="n" )
    #legend( 18, 42, rownames(qperc), horiz = TRUE, fill=cols, bty="n" )
    legend( 0, ymax+round(ymax*0.04)+2, rownames(qperc), horiz = TRUE, fill=cols, bty="n" )
    #right of the right hand side plot
    #legend( 37, 26, rownames(qperc), fill=cols, bty="n" )

    par( mfrow=c(1,1) )
    
    #if( plottofile ) {
    #    dev.off()
    #}
    
    return( allfreq )
}

# the same as plot_basecall_per_cycle, but log transformed
plot_basecall_per_cycle_log <- function( allqfreq=NULL ) {
    trace.enter("plot_basecall_per_cycle_log");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "basecallPerCycleLog.png", sep="" )
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(allqfreq)*480 )
    #}
    
    par( mfrow=c(1,length(allqfreq)) )
    
    cols = c( "red", "blue", "chartreuse4", "darkgoldenrod3", "black" )
    lwd = 2
    
    j = 1
    ymax = 0
    for( qfreq in allqfreq ) {
        qlog = log10(qfreq+1)
        
        log.info( "Plotting log transformed base call per cycle..." )
        ymax = max( ymax, max(qlog) )
        plot( 1:length(qlog[1,]), qlog[1,], xlab="Cycle", ylab="log10 base count + 1", 
                ylim=c( 0, ymax ), type="l", col=cols[1], lwd=lwd )
        for( i in 2:length(qlog[,1]) ) {
            lines( 1:length(qlog[1,]), qlog[i,], col=cols[i], lwd=lwd )
        }
        
        par( xpd=NA )
        if( length(allqfreq) > 1 ) {
            if( j == 1) {
                legend( ncol(qlog)/2, ymax+ymax*0.04, "1st mates", horiz = TRUE, bty="n" )
                title( main=paste( "Base call per cycle" ) )
            } else
                legend( ncol(qlog)/2, ymax+ymax*0.04, "2nd mates", horiz = TRUE, bty="n" )
        }
        j=j+1
    }
    par( xpd=NA )
    legend( 0, 3, rownames(qlog), horiz = TRUE, fill=cols, bty="n" )
    
    par( mfrow=c(1,1) )
    #if( plottofile ) {
    #    dev.off()
    #}
}

plot_N_distrib <- function( alldata ) {
    trace.enter("plot_N_distrib");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "Ncdf.png", sep="" )
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(alldata)*480 )
    #}
    
    lwd = 3
    col="darkgoldenrod3"
    par( mfrow=c(1,2), xpd = FALSE )
    
    qfreq = matrix(0, nrow = length(alldata[[1]]), ncol = 1)
    for( data in alldata ) {
        qfreq = qfreq + alphabetFrequency( sread(data) )[,"N"]
    }
    
    log.info( "Plotting histogram of Ns..." )
    #hist( log10(qfreq[5,qfreq[5,] > 0]), main="Histogram of 'N' occurences", 
    #        xlab="log10 number of Ns", ylab="Number of reads", col="darkgoldenrod3" )
    hist( qfreq, main="Histogram of 'N' occurences", 
            xlab="Number of Ns", ylab="Number of reads", 
            col="darkgoldenrod3", breaks=max(qfreq, 1) )
    log.info( "Plotting ecdf of Ns..." )
    #plot( ecdf( qfreq[5,] ), main="CDF of 'N' occurences", xlab="Number of Ns", ylab="Proportion of reads", 
    #        do.points = FALSE, verticals = TRUE, lwd=lwd, col=col )
    plot( ecdf( qfreq ), main="ecdf of 'N' occurences", xlab="Number of Ns", ylab="Proportion of reads", 
            do.points = FALSE, verticals = TRUE, lwd=lwd, col=col )    
    
    #if( plottofile ) {
    #    dev.off()
    #}
}

plot_polyX <- function( data, ind ) {
    trace.enter("plot_polyX");
    on.exit({ trace.exit() })
    
    reads = as.character( data )
    count = polyX_size(reads)
    bases = names(count)
    cols = brewer.pal(4,"Dark2")
    names(cols) = bases
    
    xlim = c(2,width(data)[1])
    ylim = c(0,1)
    par( mfrow=c(1,3) )
    plot( x = 1, y=1, xlim=xlim, ylim=ylim, type="n",
            main=paste("All reads"), 
            xlab=paste("Length of the tail"), ylab="Cumulative proportion of the reads" )
    for( b in bases ) {
        plot( ecdf(count[[b]][count[[b]] > 1]), xlim=xlim, ylim=ylim, add = TRUE,
                do.points = FALSE, col=cols[b], verticals = TRUE, lwd=3 )
    }
    
    plot( x = 1, y=1, xlim=xlim, ylim=ylim, type="n",
            main=paste("Aligned reads"), 
            xlab=paste("Length of the tail"), ylab="Cumulative proportion of the reads" )
    for( b in bases ) {
        plot( ecdf(count[[b]][!ind & (count[[b]] > 1)]), xlim=xlim, ylim=ylim, add = TRUE,
                do.points = FALSE, col=cols[b], verticals = TRUE, lwd=3 )
    }
    
    plot( x = 1, y=1, xlim=xlim, ylim=ylim, type="n",
            main=paste("Unaligned reads"), 
            xlab=paste("Length of the tail"), ylab="Cumulative proportion of the reads" )
    legend( x="bottomright", bases, fill=cols )
    for( b in bases ) {
        plot( ecdf(count[[b]][ind & (count[[b]] > 1)]), xlim=xlim, ylim=ylim, add = TRUE,
                do.points = FALSE, col=cols[b], verticals = TRUE, lwd=3 )
    }
    
    par( mfrow=c(1,1) )
}

# the higher the dustyScore lower the complexity of the read
plot_dustyScore_distrib <- function( alldata, update = FALSE ) {
    trace.enter("plot_dustyScore_distrib");
    on.exit({ trace.exit() })
    
    imagefilename = paste( .project$report["raw_report_dir"], "/", "dusty.png", sep="" )
    
    #if( plottofile ) {
    #    png( filename=imagefilename, width=length(alldata)*480 )
    #}
    
    lwd = 2
    col="darkgoldenrod3"
    par( mfrow=c(1,2), xpd = FALSE )
    
    dustS = calculate_dustyScore( alldata, update )
    
    if( length(dustS) == length(alldata[[1]]) ) {
        many = "all"
    } else {
        many = paste( length(dustS), "random" )
    }
    
    log.info( "Plotting dusty score distribution..." )
    plot( density(dustS), main=paste("Read complexity density for", many, "reads"), 
            xlab="Dusty score", lwd=lwd, col=col )
    log.info( "Plotting dusty score ecdf..." )
    plot( ecdf(dustS), main=paste("ecdf of read complexity for", many, "reads"), 
            xlab="Dusty score", ylab="Proportion of reads",
            do.points = FALSE, verticals = TRUE, lwd=lwd, col=col )
    par( mfrow=c(1,1) )
    
    #if( plottofile ) {
    #    dev.off()
    #}
    
}

calculate_dustyScore <- function( alldata, update ) {
    filename = .project$aux_files["dusty"];
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        
        tmax = 1000000000
        #tmax = 10000
        
        if( length(alldata) > 1 ) {
            if( as.numeric(width(alldata[[1]])[1])*length(alldata)*length(alldata[[1]]) <= tmax  ) {
                data = DNAStringSet( x=paste(sread(alldata[[1]]), sread(alldata[[2]]), sep="") )
            } else {
                log.info( "The dataset it too big, calculating dusty score for a random subset of it!" )
                index = as.logical(srswor(round(tmax/(width(alldata[[1]])[1]*length(alldata))),length(alldata[[1]])))
                #index = round(runif(round(tmax/(width(alldata[[1]])[1]*length(alldata))), min=1, max=length(alldata[[1]])))
                data = DNAStringSet( x=paste(sread(alldata[[1]])[index], sread(alldata[[2]])[index], sep="") )
            }
            
        } else {
            data = sread(alldata[[1]])
        }
        log.info( "Calculating the dusty score, this may take a while..." )
        dustS = dustyScore( data )
        save( dustS, file=filename )
        log.info( "Saved to ", filename, "" )
        
    }
    
    return( dustS )
}

# receives the top table obtained from tables(), quant is the upper quantile of the nr
# of occurences in percentage (the default is the upper 5%)
plot_repeatedreads_vs_complexity <- function( data, quant=2, logscale = FALSE ) {
    trace.enter("plot_repeatedreads_vs_complexity");
    on.exit({ trace.exit() })
    
    qt = quantile( data, probs=seq(0,1,quant/100) )
    x = data[data >= qt[length(qt)-1]]
    # the higher the dustyScore lower the complexity of the read
    y = dustyScore( DNAStringSet(names(x)) )
    if( logscale )
        plot( log10(x), log10(y), xlab="Number of occurences", 
                ylab="Dusty Score (higher score = lower complexity)")
    else 
        plot( log10(x), log10(y), xlab="Number of occurences", 
                ylab="Dusty Score (higher score = lower complexity)")
}

plot_multiplealigns_hist <- function( data, update = FALSE ) {
    trace.enter("plot_multiplealigns_hist");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report_dir,"/", "multipleAligns.RData", sep="" )
    
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
        return( alntab )
    } 
    alntab = tables( id(data), n=length(data) )
    
    log.info( "Saving to file ", filename, "" )
    save( alntab, file=filename )
}

plot_multiplealigns_vs_complexity <- function( aln, data, n = 50, update = FALSE ) {
    trace.enter("plot_multiplealigns_vs_complexity");
    on.exit({ trace.exit() })
    
    x = data[1:n]
    index = match( names(x), id(aln) )
    y = dustyScore( sread(aln)[index] )
    plot( x, y, xlab="Number of alignments", ylab="Dusty score" )
}

plot_hits <- function( total_reads, bam_filename ) {
    trace.enter("plot_hits");
    on.exit({ trace.exit() })
    
    cols=c("darkgoldenrod3","darkgoldenrod4")
    par( mfrow=c(1,2), xpd = FALSE )
    
    if( .project$pairing$type == "PE" ) {
        flags = scanBamFlag( isPaired = TRUE, isProperPair = TRUE, 
                isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                isFirstMateRead = TRUE )
    } else {
        flags = scanBamFlag( isUnmappedQuery = FALSE )
    }
    
    qnames = scanBam( bam_filename, param=ScanBamParam( flag=flags, what=("qname") ) )[[1]]$qname
    tab = table( qnames )
    rm(qnames)
    
    total_aligned = length(tab)
    unique_match = sum( tab == 1 )
    multiple_match = total_aligned-unique_match
    
    data = matrix( c(total_reads, 0, unique_match, multiple_match), ncol=2 )
    colnames(data) = c("total reads", "aligned reads")
    xc = barplot( data, beside = FALSE, col=cols, ylab="Number of reads" )
    par( xpd=NA )
    legend( xc[1], total_reads/2, "raw reads 100%", bty="n", xjust=0.5, yjust=0.5 )
    legend( xc[2], unique_match/2, paste("unique hit ~",round(unique_match*100/total_reads), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )
    legend( xc[2], unique_match+(multiple_match/2), paste("multiple hits ~",round(multiple_match*100/total_reads), "%", sep=""), bty="n", xjust=0.5, yjust=0.5 )

    hist( tab, col=cols[1], breaks=c(0.5:(max(tab)+0.5)), xlab="Number of hits", 
            ylab="Number of reads", main="Histogram of aligner hits per read" )
}

# if calling on an AlignedRead object, call plot_bars_totalunique( alignData(aln)$cigar, ...
# else if calling on a bam, call plot_bars_totalunique( bam$cigar, ...
plot_bars_totalunique <- function( cigars, hitcount, seq, index ) {
    trace.enter("plot_bars_totalunique");
    on.exit({ trace.exit() })
    
    total_reads = length( seq$seqdata )
    total_aligned = length( hitcount$top )
    unique_match = sum(index)
    # TODO, this gives the incorrect number for bam files generated by tophat
    multiple_match = sum( !index )
    unique_ungapped = sum( as.character(cigars[index]) == paste( width(seq$seqdata)[1], "M", sep="" ) )
    unique_gapped = sum( as.character(cigars[index]) != paste( width(seq$seqdata)[1], "M", sep="" ) )
    
    names = data.frame( c("raw reads", "aligned", "unique match", "ungapped" ), 
            c( "", "not aligned", "multiple match", "gapped" ), 
            c( "", "", "", "" ) )
    # data lengths
    data = data.frame( c( 100, round(total_aligned*100/total_reads), round(unique_match*100/total_reads), round(unique_ungapped*100/total_reads) ), 
            c( 0, 100-round(total_aligned*100/total_reads), round(multiple_match*100/total_reads), round(unique_gapped*100/total_reads) ), 
            c( 0, 0, 0, 0 ) )
    # relative percentages
    reldata = data.frame( c( 100, round(total_aligned*100/total_reads), round(unique_match*100/total_aligned), round(unique_ungapped*100/unique_match) ), 
            c( 0, 100-round(total_aligned*100/total_reads), round(multiple_match*100/total_aligned), round(unique_gapped*100/unique_match) ), 
            c( 0, 0, 0, 0 ) )
    
    height = round(100/nrow(data))
    col=rainbow(11, start=.7,end=.1)
    plot(c(0, 100), c(0, 100), type = "n", xlab=paste( "Percentage of the total reads (", total_reads, ")", sep="" ), 
            ylab="Percentage relative to the previous",
            main = "Read count", yaxt="n")
    for( i in 1:nrow(data) ) {
        pos = 0
        for( j in 1:ncol(data) ) {
            if( data[i,j] != 0 ) {
                if( j == 1 ) {
                    rect(0, height*(i-1), data[i,j], height*(i), col=col[i])
                    text( 0+data[i,j]/2, height*(i-1)+(height*(i)-height*(i-1))/2, paste( names[i,j], " (", reldata[i,j], "%)", sep="" ) )
                } else {
                    rect(pos,height*(i-1),pos+data[i,j],height*(i), col=col[i])
                    text( pos+data[i,j]/2, height*(i-1)+(height*(i)-height*(i-1))/2, paste( names[i,j], " (", reldata[i,j], "%)", sep="" ) )
                }
                pos = pos + data[i,j]
            }
        }
    }
}

plot_readoccurence_cdf <- function( data ) {
    trace.enter("plot_readoccurence_cdf");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "readOccurenceHist.png", sep="" )
    
    lwd = 2
    col="darkgoldenrod3"
    tmax = 1000000000
    #tmax = 10000
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(data)*480 )
    #}
    
    par( mfrow=c(1,1), xpd = FALSE )
    
    if( length(data) > 1 ) {
        if( as.numeric(width(data[[1]])[1])*length(data)*length(data[[1]]) <= tmax  ) {
            tab = tables( DNAStringSet( x=paste(sread(data[[1]]), sread(data[[2]]), sep="") ), n=length(data[[1]]) )
        } else {
            log.info( "The dataset it too big, calculating read occurence for a random subset of it!" )
            #index = round(runif(round(tmax/(width(data[[1]])[1]*length(data))), min=1, max=length(data[[1]])))
            index = as.logical(srswor(round(tmax/(width(data[[1]])[1]*length(data))),length(data[[1]])))
            tab = tables( DNAStringSet( x=paste(sread(data[[1]])[index], sread(data[[2]])[index], sep="") ), n=sum(index) )
        }
    } else {
        tab = tables( sread(data[[1]]), n=length(data[[1]]) )
    }
    
    if( length(tab$top) == length(data[[1]]) ) {
        many = "all"
    } else {
        many = paste( length(tab$top), "random" )
    }
    
    plot( ecdf( log10(tab$top) ), main=paste("ecdf of times a sequence appears for", many, "reads"), 
            xlab="log10 number of occurences", 
            ylab="Proportion of sequences", do.points = FALSE, verticals = TRUE, lwd=lwd, col=col )
    
    #if( plottofile ) {
    #    dev.off()
    #}
    
}

plot_ranges <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", sep = 0.5, ...) { 
    trace.enter("plot_ranges");
    on.exit({ trace.exit() })
    
    height <- 1 
    if (is(xlim, "Ranges"))
        xlim <- c(min(start(xlim)), max(end(xlim))) 
    bins <- disjointBins(IRanges(start(x), end(x) + 1)) 
    plot.new() 
    plot.window(xlim, c(0, max(bins) * (height + sep))) 
    ybottom <- bins * (sep + height) - height 
    rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, col = col, ...) 
    title(main)
    axis(1)
}

plot_avgbasequal_distrib <- function( alldata, update = FALSE ) {
    trace.enter("plot_avgbasequal_distrib");
    on.exit({ trace.exit() })
    
    filename = paste( .project$report["raw_report_dir"], "/", "basequalDens.png", sep="" )
    
    #if( plottofile ) {
    #    png( filename=filename, width=length(alldata)*480 )
    #}
    
    par( mfrow=c(1,length(alldata)), xpd = FALSE )
    
    lwd = 2
    col="darkgoldenrod3"
    
    j = 1
    avgquals = phred_to_avgqual( alldata, update=update )    
    for( avgq in avgquals ) {
        log.info( "Calculating average base qualities..." )
        log.info( "Plotting average base qualities..." )
        
        # because the C function already does -33
        offset = 31
        if( .project$qual_type == "FastqQuality" )
            offset = 0
        avgq = avgq - offset
        den = density(avgq)
        
        plot( den, main="", xlab="Phred scale average base call quality", 
                xlim=c(-5,45), lwd=lwd, col=col )
        #abline( v=.project$filtering_options$minqual, col="red" ) 
        
        par( xpd=NA )
        if( length(alldata) > 1 ) {
            if( j == 1)
                legend( -5, max(den$y), "1st mates", horiz = TRUE, bty="n" )
            else
                legend( -5, max(den$y), "2nd mates", horiz = TRUE, bty="n" )
        }
        j = j + 1
    }
    
    par( mfrow=c(1,1) )
    title( main="Distribution of average base call read quality" )
    
    #if( plottofile ) {
    #    dev.off()
    #}
}

# if calling on an AlignedRead do plot_mapq_distribution( as.numeric(as.character(quality(alignQuality(aln)))) )
# if the sample is single-end index can be NULL
plot_mapq_distribution <- function( bam, index ) {
    trace.enter("plot_mapq_distribution");
    on.exit({ trace.exit() })
    
    #quals = as.numeric(as.character(mapq))
    
    # mapping quality should take into account both mates so we only
    # need to plot one of them
    if( .project$pairing$type == "PE" ) {
        log.info( "Sample is paired-end, will plot the quality for only one of the mates." )
        quals = bam$mapq[index[["mates"]]]
    } else
        quals = bam$mapq
    filt = (quals != 255)
    plot( density(quals[filt]), main=paste( "Mapping quality distribution (~", 
                    round(sum(!filt)/length(quals)*100), "% had no info)", sep="" ), 
            xlab="Phred scale mapping quality" )
}
