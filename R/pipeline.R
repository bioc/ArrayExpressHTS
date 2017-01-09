#
#  Alignment, Estimation, Report & Util Routines
#

align <- function( update = FALSE, sure.i.want.to.run.this = FALSE ) {
    trace.enter("align");
    on.exit({ trace.exit() })
    
    if( (.project$aligner$type != "custom") && (update || !file.exists(paste(.project$aligner$out_dir, "/", 
                    .project$aligner$files[1], sep="")) ) ) {
        fail = tryCatch( do.call(paste( "align_", .project$aligner$type, sep="" ), 
                        list( run=sure.i.want.to.run.this )),
                error = function(e) { 
                    #
                    #
                    registerRunStepFailure(ALIGNMENT, "Error during alignment :", e);
                    stop();
                })
    } else {
        log.info( "Alignment already exists. Skipping..." )
        fail = FALSE
    }
    if( !fail ) {
        call_samtools( sure.i.want.to.run.this, update )
    }
    
    recordAlignedProperties();
    
}

align_tophat <- function( run = FALSE ) { 
    trace.enter("align_tophat");
    on.exit({ trace.exit() })
    
    indexes_dir = .project$aligner$indexes_dir
    output_dir = .project$aligner$out_dir
    
    cmds = list()
    if( .project$pairing$type == "SR" ) {
        fq_files = paste( .project$fq_files, collapse="," )
    } else if( .project$pairing$type == "PE" ) {
        fq_files = paste( paste( .project$fq_files[1], collapse="," ),
                paste( .project$fq_files[2], collapse="," ) )
    }
    options(scipen=10)
    
    # override all options with the user defined ones
    if( .project$aligner$override_options ) {
        options = .project$aligner$options
    } else {
        options = paste( "-p", .project$aligner$threads )
        for( i in 1:length(.project$aligner$options) ) {
            if( !is.null(.project$aligner$options[[i]] ) ) {
                if( !is.null(names(.project$aligner$options)) &&
                        (names(.project$aligner$options)[i] == "--mate-inner-dist") 
                        && !is.null(.project$pairing$insize) )
                    options = paste( options, names(.project$aligner$options)[i], .project$pairing$insize )
                else if( (names(.project$aligner$options)[i] == "--mate-std-dev") 
                        && !is.null(.project$pairing$insizedev) )
                    options = paste( options, names(.project$aligner$options)[i], .project$pairing$insizedev )
                else
                    options = paste( options, names(.project$aligner$options)[i], .project$aligner$options[[i]] )
            }
        }
        
        # default is phred33-quals
        #
        #--solexa-quals                          
        #--solexa1.3-quals                          (same as phred64-quals)
        #--phred64-quals                            (same as solexa1.3-quals)
        
        if( .project$qual_type == "FastqQuality" ) {
            # the tophat uses phred+33 by default,
            # which is fine for FastqQuality type
            #
            #options = paste( options, "--solexa-quals" )
        } else {
            options = paste( options, "--phred64-quals" )
        }
    }
    if( .project$count$method == "mmseq" ) {
        options = paste( options, " -G ", .project$annot$folder, "/", .project$organism, 
                ".", .project$reference$version, ".gtf",
                " --no-novel-juncs --min-isoform-fraction 0.0 --min-anchor-length 3", sep="" )
    }
    
    tophatcmd = paste(getPipelineOption("ArrayExpressHTS.tophat"), "/tophat", sep="");

    cmds = c( cmds, paste( tophatcmd, " ", options, " -o ", output_dir, " ", 
                    .project$reference$file, " ", fq_files, sep="" ) )
    
    i = 1
    for( cmd in cmds ) {
        if( i == 1 ) {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep="") )
        } else {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep=""), append = TRUE )
        }
        i = i+1
    }
    
    fail = run_cmds( cmds, run )
    return( fail )
}

align_bowtie <- function( run = FALSE ) {
    trace.enter("align_bowtie");
    on.exit({ trace.exit() })
    
    indexes_dir = .project$aligner$indexes_dir
    output_dir = .project$aligner$out_dir
    fq_files = .project$fq_files
    
    cmds = list()
    
    if( .project$pairing$type == "SR" ) {
        fq_files = paste( .project$fq_files, collapse="," )
    } else if( .project$pairing$type == "PE" ) {
        fq_files = paste( "-1", paste( .project$fq_files[1], collapse="," ),
                "-2", paste( .project$fq_files[2], collapse="," ) )
    }
    
    # override all options with the user defined ones
    if( .project$aligner$override_options ) {
        options = .project$aligner$options
    } else {
        options = paste( "-p", .project$aligner$threads )        
        
        
        #--phred33-quals    input quals are Phred+33 (default)
        #--phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
        #--solexa-quals     input quals are from GA Pipeline ver. < 1.3
        #--solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
        #--integer-quals    qualities are given as space-separated integers (not ASCII)
        
        if( .project$qual_type == "FastqQuality" ) {
            #
            #
            #options = paste( options, "--solexa-quals" )
            options = paste( options, "--phred33-quals" )
        } else {
            #options = paste( options, "--solexa1.3-quals" )
            options = paste( options, "--phred64-quals" )
        }
        
        # for working with mmseq
        if( .project$reference$type == "transcriptome" ) {
            options = paste( options, "--fullref" )
        }
        
        for( i in 1:length(.project$aligner$options) ) {
            if( !is.null(.project$aligner$options[[i]] ) ) {
                options = paste( options, names(.project$aligner$options)[i], .project$aligner$options[[i]] )
            }
        }
    }
    
    bowtiecmd = paste(getPipelineOption("ArrayExpressHTS.bowtie"), "/bowtie", sep="");
    
    cmds = c(cmds, paste( bowtiecmd, " ", options, 
                    .project$reference$file, fq_files,
                    paste( output_dir, "/", .project$aligner$files[1], sep="") ) ) # "accepted_hits.sam"
    
    i = 1
    for( cmd in cmds ) {
        if( i == 1 ) {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep="") )
        } else {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep=""), append = TRUE )
        }
        i = i+1
    }
    fail = run_cmds( cmds, run )
    return( fail )
}

align_bwa <- function( run = FALSE ) {
    trace.enter("align_bwa");
    on.exit({ trace.exit() })
    
    reference_file = paste( .project$reference$file, ".fa", sep="" )
    output_dir = .project$aligner$out_dir;
    
    cmds = list()
    fq_files = .project$fq_files
    
    # override all options with the user defined ones
    if( .project$aligner$override_options ) {
        options = .project$aligner$options[1]
    } else {
        options = paste( "-t", .project$aligner$threads )    
        for( i in 1:length(.project$aligner$options) ) {
            if( !is.null(.project$aligner$options[[i]] ) && 
                    !is.null(names(.project$aligner$options)) &&
                    names(.project$aligner$options)[i] != "hits" ) {
                options = paste( options, names(.project$aligner$options)[i], .project$aligner$options[[i]] )
            } 
        }
    }
    
    sai = paste( "fq", 1:length(fq_files), ".sai", sep="" )
    
    bwacmd = paste(getPipelineOption("ArrayExpressHTS.bwa"), "/bwa", sep="");
    
    cmds = c( cmds, paste( bwacmd, " aln ", options, " ", reference_file, " ", .project$fq_files, " > ", output_dir, "/", sai, sep="" ) )
    
    options = "-n 40"
    if( .project$pairing$type == "SR" ) {
        cmds = c( cmds, paste( bwacmd, " samse ", options, " ", 
                        reference_file, " ", output_dir, "/", sai, " ", .project$fq_files, " > ", 
                        output_dir, "/", .project$aligner$files[1], sep="" ) )
    } else {
        # -a is not used, only when there are not enough alignments to infer the insert size
        #if( !is.null(.project$pairing$insize) )
        #    options = paste( "-a", .project$pairing$insize )
        cmds = c( cmds, paste( bwacmd, " sampe ", options, " ", reference_file, " ", paste( output_dir, "/", sai, sep="", collapse=" " ), " ",
                        paste( fq_files, collapse=" "), " > ", output_dir, "/", .project$aligner$files[1], sep="" ) )
    }
    
    
    i = 1
    for( cmd in cmds ) {
        if( i == 1 ) {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep="") )
        } else {
            write( cmd, file=paste(output_dir, "/", "run_options.txt",sep=""), append = TRUE )
        }
        i = i+1
    }
    
    fail = run_cmds( cmds, run )
    return( fail )
}

filterbam = function(filename, newfile) {
    argc = 3;
    argv = c("fltbam", filename, newfile);
    result = 123;
    
    .C('fltbam', as.integer(argc), as.character(argv), as.integer(result), PACKAGE = "ArrayExpressHTS")[[3]]
}

call_samtools <- function( run, update = FALSE ) {
    trace.enter("call_samtools");
    on.exit({ trace.exit() })
    
    reference = paste( .project$reference$file, ".fa", sep="" )
    output_dir = .project$aligner$out_dir
    cmds = list()
    
    samfile = "accepted_hits.sam"
    bamfile = "accepted_hits.bam"
    sortedbam = sub( ".bam", ".sorted.bam", bamfile )
    
    samtoolscmd = paste(getPipelineOption("ArrayExpressHTS.samtools"), "/samtools", sep="");
    
    
    # produce .sam from .bam if .sam is not present
    if( !file.exists( paste( output_dir, "/", samfile, sep="") ) || update ) {
        if ( file.exists( paste( output_dir, "/", bamfile, sep="") )) {
            cmds = c( cmds, paste( samtoolscmd, " view ", output_dir, "/", bamfile, " > ", 
                            output_dir, "/", samfile, sep="" ) )
                            
        }
    }
    
    # produce .bam from .sam if .bam is not pres
    if( !file.exists( paste( output_dir, "/", bamfile, sep="" )) || update ) {
        referencefai = paste(reference,".fai",sep="")
        
        if( !file.exists(referencefai) || update ) {
            cmds = c( cmds, paste( samtoolscmd, " faidx", reference ) )
        }
        
        cmds = c( cmds, paste( samtoolscmd, " import ", referencefai, " ", output_dir, "/", samfile, " ", 
                        output_dir, "/", bamfile, sep="" ) )
    }
    
    if( !file.exists( paste( output_dir, "/", sortedbam, sep="" ) ) || update ) {
        
        # needs to be sorted to be read into R
        cmds = c( cmds, paste( samtoolscmd, " sort ", output_dir, "/", bamfile, " ", output_dir, "/", "accepted_hits.sorted", sep="" ) )
        # needed by mmseq
        cmds = c( cmds, paste( samtoolscmd, " sort -n ", output_dir, "/", bamfile, " ", output_dir, "/", "accepted_hits.sortednames", sep="" ) )
        # needs to be indexed to be read into R
        cmds = c( cmds, paste( samtoolscmd, " index ", output_dir, "/", sortedbam, sep="" ) )
    }
    
    filename = paste( output_dir, "/", bamfile, sep="" )
    newfile = sub( ".bam", ".filt.bam", filename )
    sorted = sub( ".bam", ".sorted", newfile)
    final = paste( sorted, ".bam", sep="" )
    
    if( .project$pairing == "PE" && .project$aligner$type %in% c("bowtie","tophat") ) { # filter
        
        # needs to be filtered because Bowtie/Tophat don't return only the best 
        # stratum of alignment quality when using Paired end reads
        if( !file.exists(newfile) ) {
            
            # execute commands accumulated so far
            #
            fail = run_cmds( cmds, run )
            
            # clean command list
            #
            cmds = list();
            
            # run bam mfiltering
            #
            rc0 = filterbam(filename, newfile);
            
            # log the result
            #
            log.info("FLTBAM ", filename, " RC = ", rc0);
        }
        
        # needs to be sorted to be read into R
        if( !file.exists(sorted) ) {
            cmds = c( cmds, paste( samtoolscmd, " sort ", newfile, " ", sorted, sep="" ) ) 
        }
        
        # needs to be indexed to be read into R
        
        bai = paste( final, ".bai",sep="")
        if( !file.exists(bai) ) {
            cmds = c( cmds, paste( samtoolscmd, " index ", final, sep="" ) )
        } 
    } else {
        cmds = c( cmds, paste( "ln -s ", output_dir, "/", sortedbam, " ", final, sep="" ) )
    }
    
    if( !file.exists( paste( output_dir, "/", "raw.pileup", sep="" ) ) || update ) { # pileup -c -f
        
        # new workflow
        #
        cmds = c( cmds, paste( samtoolscmd, " mpileup -f ", reference, " ", output_dir, "/", sortedbam, " > ",
                        output_dir, "/", "raw.pileup", sep="" ) )
    }
    
    fail = run_cmds( cmds, run )
    
    return(fail)
}

# reads in text files as generated by SAMtools pileup
pileup_to_dataframe <- function( update = FALSE ) {
    trace.enter("pileup_to_dataframe");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["pileup"]
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading... " )
        load( filename )
    } else {
        pileup_file = paste( .project$aligner$out_dir, "/", "raw.pileup", sep="" )
        pile = read.table( pileup_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE )
        colnames(pile) = c( "rname", "pos", "refbase", "consbase", "consq", "snpq", "RMS", "nreads", "readbases", "qualbases" )
        save( pile, file=filename )
        log.info( filename, " saved" )
    }
    
    return( pile )
}


map_quality_scale0 <- function( num ) {
    trace.enter("map_quality_scale0");
    on.exit({ trace.exit() })
    
    if (num == 33) {
        # sanger, phred+33
        return("FastqQuality");
    } else if (num == 64) {
        # solexa, illumina 1.3, illumina 1.5, phred+64
        return("SFastqQuality");
    } else {
        # default
        return("FastqQuality");
    }
}

# get fastq quality
get_fastq_quality0 <- function( fname ) {
    trace.enter("get_fastq_quality0");
    on.exit({ trace.exit() })
    
    readmax = getPipelineOption("fastqreadmax"); 
    result = 0;
    result0 = .C("checkQuality", as.character(fname), as.integer(readmax), as.integer(result))[[3]];
    
    return(map_quality_scale0(result0));
}


# get fastq quality
get_fastq_quality_scale <- function( fname ) {
    trace.enter("get_fastq_quality_scale");
    on.exit({ trace.exit() })
    
    splitarr = unlist(strsplit(fname, '/'));
    name = splitarr[length(splitarr)]
    
    outfname = paste(tempdir(), "/", name, ".tmp", sep = "");
    
    writeLines( readLines(fname, n = 8), file(outfname) )
    
    seq = readFastq(outfname)
    
    return(class(quality(seq))[1]);
}


# read fastq files 
# the fastq format defines:
# 1st line: @id:lane:tile:x:y#multiplexIndex/paired-endNumber(1|2)
# 2nd line: sequence
# 3rd: +
# 4th: quality
fastq_to_shortreadq <- function( update = FALSE ) {
    trace.enter("fastq_to_shortreadq");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["fqAsShortReadQ"]
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading... " )
        varnames = load( filename )
        log.info( paste( "loaded", paste(varnames, collapse=", ") ) )
        
        # update project data
        #
        proj0 = .project;
        proj0$qual_type <- class(quality(seqdata[[1]]))[1];
        
        assignOneProject(proj0);
        
        
    } else {
        if( .project$pairing$type == "PE" ) {
            seqdata = list()
            for( i in 1:2 ) {
                log.info( paste("Reading fastq file ", dir(.project$datadir, pattern=paste(.project$name,"*_", i,sep="")), "...", sep="") )
                seqdata[[i]] = readFastq( .project$datadir, pattern=paste(.project$name,"*_", i,sep=""), qualityType="Auto" )
            }
            log.info( "Reordering the mates..." )
            ord = IRanges::match( subseq(id(seqdata[[1]]), 1, width(id(seqdata[[1]]))-1), 
                    subseq(id(seqdata[[2]]), 1, width(id(seqdata[[2]]))-1))
            seqdata[[2]] = seqdata[[2]][ord]
            
            # if any base quality < 59 (ASCII ';') then FastqQuality
            # else SFastqQuality
            qual = class(quality(seqdata[[1]]))[1] 
        } else {
            seqdata = list()
            log.info( paste("Reading fastq file ", dir(.project$datadir, pattern=.project$name), "...", sep="") )
            seqdata[[1]] = readFastq( .project$datadir, pattern=.project$name, qualityType="Auto" )
            qual = class(quality(seqdata[[1]]))[1]
        }
        
        # update project data
        #
        proj0 = .project;
        proj0$qual_type <- qual;
        assignOneProject(proj0);
        
        save( seqdata, file=filename )
        log.info( filename, " saved" )
    }
    
    return(seqdata)
}

shortread_to_tab <- function( seqdata, update = FALSE ) {
    trace.enter("shortread_to_tab");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["shortReadQtab"]
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading... " )
        varnames = load( filename )
        log.info( paste( "loaded", paste(varnames, collapse=", ") ) )
    } else {
        if( .project$pairing$type == "PE" ) {
            log.info( "Joining mates..." )
            reads = DNAStringSet( paste( sread(seqdata[[1]]), sread(seqdata[[2]]), sep="" ) )
            log.info( "Calculating table..." )
            tab = tables(reads, n=length(reads) )
        } else {
            reads = sread(seqdata[[1]])
            tab = tables(reads, n=length(reads) )
        }
        
        save( tab, file=filename )
        log.info( filename, " saved" )
    }
    
    return(tab)
}

# does not support stranded reads
shortread_to_dataframe <- function( data, update = FALSE, return = TRUE ) {
    trace.enter("shortread_to_dataframe");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["fqAsDataFrame"]
    
    if( file.exists(filename) && !update ) {
        log.info( paste( "File ", filename, " exists." ) )
        if( return ) {
            log.info( "Loading... " )
            varnames = load( filename )
            log.info( paste( "loaded", paste(varnames,collpase=", ") ) )
        } else {
            seqframe = NULL
        }
    } else {
        seqframe = data.frame( id=as.character( id(data) ), seq=as.character( sread(data) ), quality=as.character( quality(data)@quality ) )
        save( seqframe, file=filename )
        log.info( filename, " saved" )
    }
    
    return( seqframe )
}

# parses a read ID as comes out of Solexa into a data.frame
stringids_to_dataframe <- function( ids, update = FALSE, return = TRUE ) {
    trace.enter("stringids_to_dataframe");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["idsAsDataFrame"]
    
    if( file.exists(filename) && !update ) {
        log.info( paste( "File ", filename, " exists." ) )
        if( return ) {
            log.info( paste( "Loading..." ) )
            load( filename )
        } else {
            idtab = NULL
        }
    } else {
        lines = strsplit( ids, ":" )
        end = sapply( lines, function(x) x[5] )
        endline = strsplit( end, "/" )
        idtab = data.frame( instrumentName=sapply( lines, function(x) x[1] ),
                flowcellLane=sapply( lines, function(x) x[2] ),
                tile=sapply( lines, function(x) x[3] ),
                clusterX=sapply( lines, function(x) x[4] ),
                clusterY=sapply( endline, function(x) x[1] ),
                sample=sapply( endline, function(x) x[2] ), stringsAsFactors = FALSE )
        save( idtab, file=filename )
        log.info( "Saved to ", filename )
    } 
    return( idtab )
}

# build an AlignedRead object from a list obtained from a bam file
bamlist_to_alignedread <- function( bam, update = FALSE, return = TRUE ) {
    trace.enter("bamlist_to_alignedread");
    on.exit({ trace.exit() })
    
    tags = .project$aligner$bamtags
    filename = .project$aux_files["alignedRead"]
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        if( .project$pairing == "SR" ) {  # single read
            
            df = data.frame( bam[["flag"]], bam[["cigar"]] )
            tagn = vector()
            for( tag in tags ) {
                if( length(bam[["tag"]][[tag]]) == length(bam[["flag"]]) ) {
                    log.info( paste( "tag", tag, "is present in the bam" ) )
                    df = cbind( df, bam[["tag"]][[tag]] )
                    tagn = c( tagn, tag )
                }
            }
            names( df ) = c( "flag", "cigar", tagn )
            metadata = data.frame( labelDescription=c( "flag", "cigar", tagn ), row.names=c("flag", "cigar", tagn) )
            aligndata = AlignedDataFrame( df, metadata )
            if( .project$reference$type == "genome" ) {
                log.info( "Checking chromosome name representation..." )
                chrnames = get_chr_representation( bam[["rname"]] )
            } else {
                chrnames = as.character(bam[["rname"]])
            }
            aln = AlignedRead( bam[["seq"]], BStringSet( bam[["qname"]] ), do.call( .project$qual_type, 
                            list( BStringSet( bam[["qual"]] ) ) ), chrnames, bam[["pos"]], bam[["strand"]], 
                    FastqQuality( as.character(bam[["mapq"]]) ), aligndata )
        } else { # paired reads
            log.warning( "Warning: not implemented for paired reads!" )
            aln = NULL
        }
        
        log.info( "Saving...")
        save( aln, file=filename )
        log.info( "Saved to file ", filename )
    }
    
    if( return )
        return(aln)
    return( NULL )
}

# instead of putting in a list like seq$seqdata, seq$seqtab these are separated for aln and alntab because
# of the size they take
# if calling on Aligned read do xstringset_to_hitcount( id(aln) )
xstringset_to_hitcount <- function( data, update = FALSE, return = TRUE ) {
    trace.enter("xstringset_to_hitcount");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["hitCount"]
    
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists. Loading..." )
        load( filename )
    } else {
        hitcount = tables( data, n=length(data) )
        save( hitcount, file=filename )
        log.info( "Saved to file ", filename )
    }
    
    if( return )
        return(hitcount)
    return( NULL )
}

bam_to_occurrencetab <- function( bam, update = FALSE, return = TRUE ) {
    trace.enter("bam_to_occurrencetab");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files$seqOccurrenceTab
    
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists. Loading..." )
        load( filename )
    } else {
        octab = tables( bam$seq, n=length(bam$seq) )
        save( octab, file=filename )
        log.info( "Saved to file ", filename )
    }
    
    if( return )
        return(octab)
    return( NULL )
}

# from a bam file to a list
bam_to_list <- function( update = FALSE, return = TRUE ) {
    trace.enter("bam_to_list");
    on.exit({ trace.exit() })
    
    bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )
    rdata_filename = .project$aux_files["bamNamesAsList"]
    if( file.exists( rdata_filename ) && !update ) {
        log.info( "File ", rdata_filename, " exists." )
        if( return ) {
            log.info( "Loading..." )
            load( rdata_filename )
        }
    } else {
        log.info( paste( "Reading bam file ", bam_filename, "..." ) )
        tags = .project$aligner$bamtags
        bam = scanBam( bam_filename, param=bamParam( tags=tags ) )[[1]]
        save( bam, file=rdata_filename )
        log.info( "Saved data in ", bam_filename, " into ", rdata_filename )
    }
    
    # this can be deleted later, so far it's only being used in the MouseCrosses project
    filename = paste(.project$aligner$out_dir, "/", "alignedQnames.RData", sep="" )
    if( !file.exists( filename ) ) {
        qn = bam$qname
        save( qn, file=filename )
    }
    
    if( return )
        return( bam )
    return( NULL )
}

bamParam <- function( what=scanBamWhat(), tags=character(), firstmate = FALSE ) {
    trace.enter("bamParam");
    on.exit({ trace.exit() })
    
    if( .project$pairing$type == "SR" ) 
        param = ScanBamParam( flag=scanBamFlag( isUnmappedQuery = FALSE ), tag=tags, what=what ) 
    else {
        if( firstmate )
            param = ScanBamParam( flag=scanBamFlag( isFirstMateRead = TRUE, isUnmappedQuery = FALSE, 
                            isProperPair = TRUE, isPaired = TRUE, hasUnmappedMate = FALSE ), tag=tags, what=what )
        else
            param = ScanBamParam( flag=scanBamFlag( isUnmappedQuery = FALSE, isProperPair = TRUE, 
                            isPaired = TRUE, hasUnmappedMate = FALSE ), tag=tags, what=what )
    }
    
    return( param )
}

call_tmm <- function( run = FALSE, e ) {
    trace.enter("call_tmm");
    on.exit({ trace.exit() })
    
    nfac = calcNormFactors(exprs(e))
    exprs(e) = exprs(e) * nfac
    
    return(e)
}

call_none <- function( run = FALSE, e = NULL ) {
    trace.enter("call_none");
    on.exit({ trace.exit() })
    
    log.info( "doing nothing" )
    
    return(e)
}

call_normalisator <- function( eset, run = FALSE ) {
    trace.enter("call_normalisator");
    on.exit({ trace.exit() })
    
    res = tryCatch( do.call(paste( "call_", .project$count$normalisation, sep="" ), 
                    list( run=run, e=eset )),
            error = function(e) { 
                #
                #
                registerExpStepFailure(ASSEMBLING_ESET, "Error during normalisation: ", e);
            } )
    return(res)
}

call_estimator <- function( update = FALSE, run = FALSE ) {
    trace.enter("call_estimator");
    on.exit({ trace.exit() })
    
    filename = paste( .project$aligner$out_dir, "/", 
            .project$count$method, "_", .project$count$files["transcript"], sep="" )
    if( update || !file.exists(filename) ) {
        fail = tryCatch( do.call(paste( "call_", .project$count$method, sep="" ), 
                        list( run=run )),
                error=function(e) { 
                    #
                    #
                    registerRunStepFailure(ESTIMATION, "Error during estimation :", e);
                    stop();

                });
    } else {
        log.info( "Estimates ", filename, " already exists. Skipping..." )
        fail = FALSE
    }
}

call_cufflinks <- function( gtf = TRUE, run = FALSE, update = FALSE ) {
    trace.enter("call_cufflinks");
    on.exit({ trace.exit() })
    
    cmds = list()
    samf = paste( .project$aligner$out_dir, "/", .project$aligner$files[1], sep="" )
    ssamf = sub( "sam", "sorted.sam", samf )
    usamf = sub( "sam", "unsorted.sam", samf )
    fail = FALSE
    
    cufflinkscmd = paste(getPipelineOption("ArrayExpressHTS.cufflinks"), "/cufflinks", sep="");

    
    filename = paste( .project$aligner$out_dir, "/", 
            .project$count$method, "_", .project$count$files["transcript"], sep="" )
    filename2 = paste( .project$aligner$out_dir, "/", 
            .project$count$method, "_", .project$count$files["gene"], sep="" )
    filename3 = paste( .project$aligner$out_dir, "/", 
            .project$count$method, "_", sub( "expr", "gtf", .project$count$files["transcript"]), sep="" )
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists." )
    } else { 
        if( .project$aligner$type != "tophat" && .project$aligner$type != "custom" ) {
            
            registerExpStepWarning(PARAMETER_SANITY, 
                    "Cufflinks should be run on the output of TopHat only.");
            
            log.warning( "Cufflinks should be run on the output of TopHat only.")
            cmds = c( cmds, paste( "sort -k 3,3 -k 4,4n", samf, ">", ssamf ) )
            cmds = c( cmds, paste("mv", samf, usamf) )
            cmds = c( cmds, paste("mv", ssamf, samf) )
        }
        
        options = paste( "-p", .project$aligner$threads )
        if( !is.null( .project$aligner$options[["--max-intron"]] ) ) {
            options = paste( options, "-I", .project$aligner$options[["--max-intron"]] )
        }
        if( gtf ) {
            options = paste( options, " -G ", .project$annot$folder, "/", .project$organism, 
                    ".", .project$reference$version, ".gtf", sep="" )
        }
        
        #options = paste( options, " -o ", .project$aligner$out_dir, sep="" );
        
        cmds = c( cmds, paste( cufflinkscmd, " ", options, " ", .project$aligner$out_dir, "/", .project$aligner$files[1], sep="" ) )
        
        # new workflow
        #
        cmds = c( cmds, paste( "mv isoforms.fpkm_tracking", filename ))
        cmds = c( cmds, paste( "mv genes.fpkm_tracking", filename2 ))        
        
        cmds = c( cmds, paste( "mv transcripts.gtf", filename3 ))
        
        fail = run_cmds( cmds, run )
    }
    return( fail )
}

call_mmseq <- function( run = FALSE, update = FALSE ) {
    trace.enter("call_mmseq");
    on.exit({ trace.exit() })
    
    cmds = list()
    fail = FALSE
    
    bam2hitscmd = paste(getPipelineOption("ArrayExpressHTS.bam2hits"), "/bam2hits", sep="");
    mmseqcmd = paste(getPipelineOption("ArrayExpressHTS.mmseq"), "/mmseq", sep="");

    
    filename = paste( .project$aligner$out_dir, "/", "hits_file", sep="" ) 
    if( file.exists(filename) && !update && file.info(filename)$size != 0 ) {
        log.info( "File ", filename, " exists." )
    } else {
           cmds = c( cmds, paste( bam2hitscmd, " -m \"(E\\S+).*gene:(E\\S+)\" 1 2 ",
                        .project$reference$file, ".fa", " ", 
                        .project$aligner$out_dir, "/", .project$aligner$files[3], " > ", filename, sep="") )
    }
    
    mmseqfilename = paste( .project$aligner$out_dir, "/", "mmseq_transcripts.expr", sep="" )
    if( file.exists(mmseqfilename) && !update ) {
        log.info( "File ", mmseqfilename, " exists." )
        fail = 0
    } else {
        mmseqfilename = paste( .project$aligner$out_dir, "/", "mmseq_out", sep="" )
        
        cmds = c( cmds, paste(mmseqcmd, " ", filename, " ", mmseqfilename, sep="") )
        cmds = c( cmds, paste("mv ", .project$aligner$out_dir, "/mmseq_out.mmseq ", 
                        .project$aligner$out_dir, "/mmseq_transcripts.expr", sep="") )
        cmds = c( cmds, paste("mv ", .project$aligner$out_dir, "/mmseq_out.gene.mmseq ", 
                        .project$aligner$out_dir, "/mmseq_genes.expr", sep="") )
        fail = run_cmds( cmds, run )
    }
    return( fail )
}

call_count <- function( run = FALSE, update = FALSE ) {
    trace.enter("call_count");
    on.exit({ trace.exit() })
    fail = FALSE
    
    countfile =  paste( .project$aligner$out_dir, "/", "count_transcripts.expr", sep="" )
    if( file.exists(countfile) && !update ) {
        log.info( "File ", countfile, " exists." )
        fail = 0
    } else {
        filename = paste( .project$aligner$out_dir, "/", .projects[[1]]$aligner$files[2], sep="" )
        log.info( "Reading bam header" )
        header = scanBamHeader( filename )[[1]]$targets
        
        log.info( "Reading bam file" )
        bam = scanBam( filename, param=bamParam( what=c("qname", "rname"), 
                        firstmate = TRUE ) )[[1]]
        
        dups = duplicated(bam$qname) | rev(duplicated(rev(bam$qname)))
        bam$qname = bam$qname[!dups]
        bam$rname = bam$rname[!dups]
        
        
        trans = bam$rname
        count = table(trans)
        tab = data.frame( trans_id=names(count), count=as.numeric(count), 
                length=header[match(names(count), names(header))] )
        rownames(tab) = NULL
        
        filenames = paste( .project$aligner$out_dir, "/", .project$count$method, 
                "_", .project$count$files, sep="" )
        names(filenames) = names(.project$count$files)
        write.table( tab, file=filenames["transcript"], quote = FALSE, )
        log.info( "Saved to ", filenames["transcript"] )
        
        bam$rname = as.character(bam$rname)
        
        trans = bam$rname;
        
        count = table(trans)
        
        tab = data.frame( gene_id=names(count), count=as.numeric(count), 
                length=header[match(names(count), names(header))] )
        rownames(tab) = NULL
        write.table( tab, file=filenames["gene"], quote = FALSE, )
        log.info( "Saved to ", filenames["gene"] )
        
    }
    return(fail)
}


cufflinks_to_dataframe <- function() {
    trace.enter("cufflinks_to_dataframe");
    on.exit({ trace.exit() })
    
    # WARNING!! for old versions edit the file first: add "length" at the end of the 1st line    
    filename = paste( .project$aligner$out_dir, "/", "transcripts.expr", sep="" )
    ctab = read.table( filename, header = TRUE, stringsAsFactors = FALSE )
    return(ctab)
}

cufflinks_to_iranges <- function() {
    trace.enter("cufflinks_to_iranges");
    on.exit({ trace.exit() })

    message("cufflinks_to_iranges is deprecated in favor of cufflinks_to_granges");
    message("cufflinks_to_iranges returns NULL");
    
    #ctab = cufflinks_to_dataframe()
    #ctab = split( ctab, ctab$chr ) 
    #ctabIR = lapply( ctab, function(x) RangedData( IRanges(start = x$left, end=x$right, names=x$trans_id ), strand=rep(1,length(x$trans_id)) ) )
    #ctabIR = do.call( RangedDataList, ctabIR )
    #return( ctabIR )

    return(NULL)
}

cufflinks_to_granges <- function() {
    trace.enter("cufflinks_to_granges");
    on.exit({ trace.exit() })
    
    cuff = cufflinks_to_dataframe()
    cuffgr = GRanges(seqnames = Rle(cuff$chr), 
            ranges = IRanges(cuff$start, width = (cuff$end-cuff$start), names = cuff$transcript_id ) )
}

mmseq_to_dataframe <- function() {
    trace.enter("mmseq_to_dataframe");
    on.exit({ trace.exit() })
    
    mm = read.table( paste( .project$aligner$out_dir, "/", "mmseq_out.mmseq", sep=""), header = TRUE, stringsAsFactors = FALSE )
    return( mm )
}

cufflinks_to_countmatrix <- function() {
    trace.enter("cufflinks_to_countmatrix");
    on.exit({ trace.exit() })
    
    ctab = cufflinks_to_dataframe()
    countm = matrix( ctab$RPKM, dimnames=list( ctab$trans_id, .project$name ) )
    return( countm )
}

# output a matrix of counts per feature
count_reads_per_feature <- function( alnIR, annot, annotIR, mr.algo=c("discard"), update = FALSE, return = TRUE ) {
    trace.enter("count_reads_per_feature");
    on.exit({ trace.exit() })
    
    feature = .project$count$feature;
    
    filename = sub( ".RData", paste( feature, ".RData", sep=""), .project$aux_files["countsPerFeature"] )
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        counts = findOverlaps(alnIR, annotIR)
        if( feature=="isoform" ) { # call cufflinks
            cmd = paste( "cufflinks", .project$count_options, .project$aligner$files[1] )
            #system( cmd ) # uncomment to run
        } else if( feature == "exon" ){
            counts = sapply( counts, as.matrix )
            tab = sapply( counts, function(x) table( x[,"subject"] ) )
            tab = sapply( tab, sort )
            allc = sapply( names(tab), function(x) {
                        featnames = names(annotIR[[x]])
                        matrix( tab[[x]], dimnames=list( featnames[as.integer(names(tab[[x]]))], "counts" ) )
                    })
        } else if( feature == "gene" ) {
            
        } else {
            log.info( "Not implemented!" )
            return( NULL )
        }
        counts = list(indexes=counts,perExon=allc)
        save( counts, file=filename )
        log.info( filename, " saved" )
    }
    
    if( return )
        return( counts )
    return( NULL )
}

transcriptome_to_dataframe <- function() {
    trace.enter("transcriptome_to_dataframe");
    on.exit({ trace.exit() })
    
    reference_file = paste( .project$aligner$indexes_dir, "/", .project$organism, ".", 
            .project$reference$version, ".cdna.chromosome.fa", sep="" )
    lengths = fasta.seqlengths( reference_file )
    tnames = names( lengths )
    
    # this should be used by default for cdna files downloaded from ENSEMBL
    tnames = strsplit( tnames, " " )
    
    chrs = strsplit( sapply(tnames, function(x) x[3] ), ":" )
    transtab = data.frame( id=sapply(tnames, function(x) x[1]), 
            chr=sapply(chrs, function(x) x[3]), 
            start=sapply(chrs, function(x) x[4]),
            end=sapply(chrs, function(x) x[5]),
            strand=sapply(chrs, function(x) x[6]),
            gene=sapply(tnames, function(x) strsplit(x[4], ":")[[1]][2]) )
    
    return( transtab )
}

get_annotation_RData_name <- function( origfile, type, filtered ) {
    sub( ".RData", paste(".", type, ".", filtered, ".RData",sep=""), origfile )
}

# either from biomart or from a gff file, returns a list of data.frames by chromosome
get_annotation <- function( .project, type="any", split = FALSE, update = FALSE, returnresult = TRUE, filter = TRUE ) {
    trace.enter("get_annotation");
    on.exit({ trace.exit() })
    
    ANYTYPE = "any";
    GTFTYPE = "gtf";
    TEXTTYPE = "txt";
    BIOMARTTYPE = "biomaRt";
    
    if (filter) {
        filtered = "filtered";
    } else {
        filtered = "unfiltered";
    }
    
    for( t in c("biomaRt", "txt", "gtf") ) {
        rfilename = get_annotation_RData_name(.project$aux_files["annot"], t, filtered);
        if( file.exists(rfilename) && !update ) {
            log.info( "File ", rfilename, " exists." )
            log.info( "Loading..." )
            load( rfilename )
            return( annot )
        } 
    }
    
    column_names = c("chr","strand", "geneid", "gene_start", "gene_end",
            "transcript", "transcript_start", "transcript_end", "biotype",
            "id", "start", "end", "feature" )
    
    annot = NULL;
    
    if( type == TEXTTYPE || type == ANYTYPE ) { # as produced by the FetchEnsembl perl script
        
        filename = dir( .project$annot$folder, pattern="txt$", full.names = TRUE )
        
        if( length(filename) != 0 ) {
            
            log.info( "Using file ", filename[1], " for annotation." )
            
            annot = read.table( filename, header = FALSE, sep="\t", fill = TRUE, stringsAsFactors = FALSE,
            col.names=column_names )
            
            type = TEXTTYPE;
            
        } else {
            log.info( "Could not find a txt file in the annotation folder." )
        }
    
    } else if( type == GTFTYPE || type == ANYTYPE ) { # as downloaded from the Ensembl ftp, it's very slow
    
        filename = dir( .project$annot$folder, pattern="gtf", full.names = TRUE )
        
        if( length(filename) != 0 ) {
            
            log.info( "Using file ", filename[1], " for annotation." )
            annot = gtf_to_dataframe( filename[1], names=column_names, filter = filter )
            chrs = unique( annot$chr )
            
            if (filter) {
                chrs = filter_chr( chrs )
            }
            
            annot = annot[annot$chr %in% chrs,]
            
            if( split ) {
                annot = split( annot, annot$chr )
            }
        
            type = GTFTYPE;
        
        } else {
            log.info( "Could not find a gtf file in the annotation folder." )
            download_annotation()
        }
    
    } else if( type == BIOMARTTYPE || type == ANYTYPE ) { # only works with mouse, human and drosophila pseudoobscura
        
        attrib = c( "chromosome_name", "strand", "ensembl_gene_id", "start_position", "end_position",
                "ensembl_transcript_id", "transcript_start", "transcript_end", "gene_biotype",
                "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end" )
        mart = "ensembl"
        
        dataset = NULL;
        
        # static info, don't change
        if( length(grep("musculus", .project$organism)) ) {
            dataset = "mmusculus_gene_ensembl"
        } else if( length(grep("sapiens", .project$organism)) ) {
            dataset = "hsapiens_gene_ensembl"
        } else if( length(grep("pseudoobscura", .project$organism)) ) {
            dataset = "dpseudoobscura_eg_gene"
            mart = "metazoa_mart_5"
            attrib = c( "chromosome_name", "strand", "ensembl_gene_id", "start_position", "end_position",
                    "ensembl_transcript_id", "transcript_start", "transcript_end", "biotype",
                    "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end" )
        
        } else {
            log.info( "Failed: biomaRt dataset name not set!" )
            dataset = NULL;
        }
        
        if (!is.null(dataset)) {
            
            ensembl = useMart( mart, dataset=dataset )
            
            chrs = getBM( attributes="chromosome_name", mart=ensembl )$chromosome_name
            log.info( paste("Available chromosomes:", paste(chrs,collapse=", ") ) )
            
            if (filter) {
                chrs = filter_chr( chrs )
            }
            
            log.info( paste("Filtering chromosomes. Will download annotation for:", paste(chrs,collapse=", ") ))
            log.info( "Downloading annotation for:" )
            #annot = list()
            annot = data.frame()
            for( chr in chrs ) {
                log.info( paste( chr, " ", sep="" ) )
                annot = rbind( annot, getBM( attributes=attrib, mart=ensembl, filters="chromosome_name", values=chr ) )
                #annot[[chr]] = annot[[chr]][order(annot[[chr]]$start_position, annot[[chr]]$strand),]
                #rownames(annot[[chr]]) = NULL
            }
            annot = cbind( annot, rep("exon",length(annot$strand)))
            colnames(annot) = column_names
            annot$chr = get_chr_representation( annot$chr )
            log.info( " " )
            
            type = BIOMARTTYPE;
        }
        
        #names( annot ) = get_chr_representation( chrs )
    } else {
        log.info( "Unkown annotatio type: ", type );
        log.info( "Supported types are: ", GTFTYPE, " ", TEXTTYPE, " ", BIOMARTTYPE, " ", ANYTYPE );
    }
    
    if (!is.null(annot)) {
        
        rfilename = get_annotation_RData_name(.project$aux_files["annot"], type, filtered);
        
        log.info( "Saving to file ", rfilename, "" )
        save( annot, file=rfilename )
        
        if( split ) {
            annot = split( annot, annot$chr )
        }
    
        if (returnresult) {
            return(annot)
        } else {
            return(NULL);
        }
    
    } else {
        return(NULL);
    }
    
}

download_annotation <- function( project, run = FALSE ) {
    trace.enter("download_annotation");
    on.exit({ trace.exit() })
    
    version = project$reference$version
    anfile = paste( project$annot$folder, "/",
            project$organism, ".", version, ".gtf", sep="")
    #log.info( "Searching for", anfile )
    if( file.exists(project$annot$folder) &&
            file.exists( anfile ) ) {
        log.info( "Annotation: ", anfile )
        fail = FALSE
    } else {
        log.info('Could not find ', project$annot$folder);
        ref_subver = paste( "release-", substr( version, nchar(version)-1, 
                        nchar(version) ), sep="")
        
        cmds = paste( "mkdir ", project$refdir, "/annotation/", sep="" )
        cmds = c( cmds, paste( "mkdir ", project$annot$folder) ) 
        run_cmds( cmds, run )
        cmds = paste( "curl -O ftp://ftp.ensembl.org//pub/", ref_subver, 
                "/gtf//", tolower(project$organism), "/", project$organism, ".", 
                version, ".gtf.gz > ", project$organism, ".", 
                version, ".gtf.gz", sep="" ) 
        cmds = c(cmds, paste( "gunzip ", project$organism, ".", 
                        version, ".gtf.gz", sep="" ) )
        cmds = c(cmds, paste( "mkdir ", project$annot$folder) )
        cmds = c(cmds, paste( "mv ", project$organism, ".", version, ".gtf ",
                        project$annot$folder, sep="" ) )
        fail = run_cmds( cmds, run )
    }
    return( fail )
}

# from a data.frame into IRanges
annot_to_iranges <- function( annot, feature=NULL, extra="strand", update = FALSE, return = TRUE ) {
    trace.enter("annot_to_iranges");
    on.exit({ trace.exit() })
    
    if( is.null(feature) ) {
        feature = .project$count$feature;
    }
    
    filename = sub( ".RData", paste(feature,".RData",sep=""), .project$aux_files["annotAsIRange"] )
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        chrs = names( annot )
        #chrs = filter_chr( chrs )
        #annot = annot[annot$chr %in% chrs,]
        
        if( feature=="exon" ) {
            columns = c("strand", "ensembl_exon_id","exon_chrom_start","exon_chrom_end")
        } else if( feature=="gene" ) {
            columns = c("strand", "ensembl_gene_id","start_position","end_position")
        } else if( feature=="rRNA" ) {
            columns = c("strand", "ensembl_gene_id","start_position","end_position")
            annot = lapply( annot, function(x) x[x$gene_biotype == "rRNA",] )
        } else {
            log.info( "Not implemented yet!" )
            return( NULL )
        }
        
        for( chr in names(annot) ) {
            annot[[chr]] = unique( annot[[chr]][,columns] )
            names(annot[[chr]]) = c("strand", "id", "start", "end")
            #annot = split( annot, annot$chromosome_name )
            #annot = lapply( annot, function(x) split(x, x$strand) )
            #annotIR = lapply( annot, function(x) lapply(x, function(y) IRanges(start = y[,"start"], end=y[,"end"], names=y[,"id"]) ) )
        }
        
        #annotIR = lapply( annot, function(y) RangedData( IRanges(start = y[,"start"], end=y[,"end"], names=y[,"id"]), y[,"strand"]) )
        annotIR = lapply(1:length(annot), function(i) { y <- annot[[i]]; 
                                GRanges( rep(names(annot[i]), length(y[,"start"])),  
                                IRanges(start = y[,"start"], end= y[,"end"], names=y[,"id"]), 
                                y[,"strand"]) } )
        #annotIR = do.call( RangedDataList, annotIR )
        annotIR = do.call( GRangesList, annotIR )
        save( annotIR, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    if( return )
        return( annotIR )
    return( NULL )
}

# make an XStringSet with only the unique sequences in the ShortReadQ or AlignedRead object
shortread_to_unique <- function( data, update = FALSE, return = TRUE ) {
    trace.enter("shortread_to_unique");
    on.exit({ trace.exit() })
    
    types= c( "sread", "id" )
    
    filename = sub( ".RData", paste( class(data)[1], ".RData", sep="" ), .project$aux_files["uniqueReads"] )
    useqs = list()
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        if( return ) {
            log.info( "Loading..." )
            load( filename )
        } else {
            useqs = NULL
        }
        log.info( "" )
    } else {
        useqs = lapply( types, function(type) unique( do.call( type, list(data) ) ) )
        save( useqs, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    return( useqs )
}

make_filt_indexes <- function( annot, update = FALSE, return = TRUE ) {
    trace.enter("make_filt_indexes");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["indexes"]
    
    bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )    
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        index = list()
        
        if( .project$pairing$type == "PE" ) {
            bam = scanBam( bam_filename, param=ScanBamParam( flag=scanBamFlag(isPaired = TRUE, 
                                    isProperPair = TRUE, isUnmappedQuery = FALSE, 
                                    hasUnmappedMate = FALSE ), what=("flag") ) )[[1]]
        }
        # for this to run correctly:
        #     one has to be sure that the bam import included only reads with the mate mapped (done if using bamParam())
        #    in case of paired-end reads the bam file has mates together
        log.info( "Getting mate indexes..." )
        if( .project$pairing$type == "PE" ) {
            # bwa's SAM files have:
            #     1 line per read with optional alignments in the XA tag
            #    1 line per mate in a pair
            # tophat/bowtie alignments have:
            #     n lines per read with all possible alignments
            #    1 line per mate in a pair
            # all aligners strip the /1 /2 terminations from the names of the reads
            index[["mates"]] = flag_is_secondmate( bam$flag )    
        } else {
            index[["mates"]] = rep(0, length(bam$qname))
        }
        # duplicated sequences, leaving the 1st occurrence
        log.info( "Getting duplicated indexes..." )
        index[["dupfilt"]] = duplicated(bam$seq)
        # all sequences that appear more than once
        index[["dup"]] = index[["dupfilt"]] & duplicated(rev(bam$seq))
        
        log.info( "Getting multiple-hit reads..." )
        # reads that align multiple times to the reference
        if( .project$aligner$type == "tophat" ) {
            index[["multhit"]] = duplicated(bam$qname) & rev(duplicated(rev(bam$qname)))  
        } else if( .project$aligner$type == "bwa" ) {
            index[["multhit"]] = !is.na(bam$tag$XA)
        } else {
            log.info( "I don't know this aligner" )
        }
        
        # reads above or equal to the min quality
        avgquals = phred_to_avgqual( bam$qual )
        index[["minqual"]] = (avgquals >= .project$filtering_options[["minqual"]])
        # reads above or equal to the min map quality (watch out, some aligners might assign qual 0 to multiple match reads
        index[["minmapq"]] = (bam$mapq >= .project$filtering_options[["minmapq"]] )
        # reads with less or equal to maxN number of Ns
        filt = nFilter( .project$filtering_options[["maxN"]] )
        index[["maxN"]] = filt( bam$seq )
        
        # reads that have less or equal than the max size of a polyX tail
        count = polyX_size( bam$seq )
        countab = data.frame( count[[1]], count[[2]], count[[3]], count[[4]] )
        if( .project$filtering_options$maxpol < 1 ) {
            maxpol = round(.project$filtering_options$maxpol * width( bam$seq[1] ))
        } else {
            maxpol = .project$filtering_options$maxpol
        }
        countab = countab > maxpol
        countab = rowSums( countab )
        index[["notpoly"]] = !countab
        
        # reads that map to the chromosome of interest (watch out for multiple reads aligning to unwanted regions)
        index[["chrOI"]] = !(bam$rname %in% .project$filtering_options[["chr_ignore"]])
        
        index[["gapped"]] = grep( "N", bam$cigar )
        # reads that map to ribosomal DNA (the same warning as for chrOI)
        #ribIR = annot_to_iranges( annot, feature="rRNA" )
        #bamIR = cigarToIRangesListByRName( bam$cigar, bam$rname, bam$pos )
        #names( bamIR ) = get_chr_representation( names(bamIR) )
        #hits = findOverlaps( bamIR, ribIR )
        # TODO: how to map back to the bam file? cigarToIRangesListByRName will split reads with gaps into 2 or more intervals!
        # GenomicRnages doesn't seem to do this, but I don't know how to represent the annotation then
        # a temporary solution would be to ignore gapped reads (control with the gapped filter option in .project)
        #sapply(hits, function(x) if( length(x@matchMatrix[,1]) > 0  ) length(unique(x@matchMatrix[,1])) else 0 )
        #index[["ribosomal"]] 
        
        save( index, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    if( return )
        return( index )
    return( NULL )
}

# global defs
#
freq00101 = list()

plot_raw_report <- function( seqdata, tab, update = FALSE ) {
    trace.enter("plot_raw_report");
    on.exit({ trace.exit() })
    
    width = 480;
    height = 480;
    
    #size = 400
    algn = "left" 
    
    report_dir = .project$report["raw_report_dir"]
    filename = "plotRawReport"
    if( !file.exists(paste(report_dir, "/", filename, ".html", sep="")) || update ) {
        
        HTMLInitFile( outdir=report_dir, filename=filename, Title="Quality report on raw reads" )
        HTML( "<h1>Run info</h1>" )
        HTML( paste("Project name:", .project$name) )
        pos = gregexpr( ":", id(seqdata[[1]])[1] )[[1]]
        HTML( paste("Lane number(s):", substr( id(seqdata[[1]])[1], pos[1]+1, pos[2]-1 ) ) )
        HTML( paste("Run length:", width(seqdata[[1]])[1]) )
        HTML( paste("Paired:", .project$pairing$type) )
        HTML( paste("Stranded:", .project[["stranded"]]) )
        HTML( paste("Number of reads:", length(seqdata[[1]])) )
        HTML( paste("Quality scale:", .project$qual_type) )
        HTML( "<h1>Raw reads quality</h1>" )
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 1")
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="basecall_per_cycle", 
                plotFunction = function() {
                    # the name is odd 
                    # not to clash with proper names
                    fq = plot_basecall_per_cycle( seqdata );
                    assignInNamespace('freq00101', fq, ns="ArrayExpressHTS");
                });
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 2")
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( freq00101 )*width, Height=height, 
                Align=algn, GraphFileName="basecall_per_cycle_log", 
                plotFunction = function() {
                    plot_basecall_per_cycle_log( allqfreq = freq00101 );
                });
        
        
        # clean the freq00101
        assignInNamespace('freq00101', list(), ns="ArrayExpressHTS");
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 3")
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="N_distrib", 
                plotFunction = function() {
                    plot_N_distrib( seqdata );
                });
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 4")
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="dustyScore_distrib", 
                plotFunction= function(){
                    plot_dustyScore_distrib( seqdata );
                });

        # memory usage
        print.memory.usage("Memory: plot_raw_report step 5")
        
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="quality_per_cycle", 
                plotFunction= function(){
                    plot_quality_per_cycle( seqdata );
                });
        
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 6")

        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="avgbasequal_distrib", 
                plotFunction= function(){
                    plot_avgbasequal_distrib( seqdata );
                });
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 7")
        
        myHTMLplot( GraphDirectory=report_dir, Width=length( seqdata )*width, Height=height, 
                Align=algn, GraphFileName="readoccurence_cdf", 
                plotFunction=function(){
                    plot_readoccurence_cdf( seqdata )
                });
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 8")
        
        HTMLEndFile()
        
        # memory usage
        print.memory.usage("Memory: plot_raw_report step 9")
    
    } else {
        log.info( "File exists, not updated" )
    }
}

# returns the number of raw reads in the fastq file
number_raw_reads <- function( update = FALSE ) {
    trace.enter( "number_raw_reads" );
    on.exit({ trace.exit() })
    
    filename = paste(.project$projectdir,"/number_raw_reads.RData", sep="")
    if( file.exists(filename) & !update ) {
        log.info( "Loading ", filename )
        load( filename )
    } else {
        log.info( "Getting umber of reads from ", .project$fq_files[1] );
        nr = as.integer( system( paste( "wc -l <", .project$fq_files[1] ), intern = TRUE ) ) / 4
        save( nr, file=filename )
    }
    return(nr);
}

# count how many reads were aligned
number_aligned_reads <- function( update = FALSE ) {
    trace.enter("number_aligned_reads");
    on.exit({ trace.exit() })
    
    filename = paste(.project$aligner$out_dir, "/number_aligned_reads.txt", sep="")
    if( file.exists(filename) && !update ) {
        log.info( "File with number of aligned reads  already exists. Loading..." );
        areads = read.table( filename )[1,1]
    } else {
        log.info( "Counting aligned reads for ", .project$name )
        hits_file = paste( .project$aligner$out_dir, "/hits_file", sep="" )
        if( file.exists(hits_file) ) { # mmseq, already filtered
            log.info( "Looking at the hits file" );
            cmds = paste( "grep \">\" ", hits_file, " | wc -l >> ", filename, sep="" )
            log.info( cmds );
            system( cmds )
            areads = read.table( filename )[1,1]
        } else { # other 
            log.info( "Looking at the BAM file" );
            bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )
            if( file.exists(bam_filename) ) {
                bam = scanBam( bam_filename, param=bamParam( what="qname", firstmate = TRUE ) )[[1]]
                areads = length( unique(bam$qname) )
                rm(bam)
                save( areads, file=filename )
            } else {
                cat( "\t", bam_filename, " not found" )
            }
        }
    }
    
    return(areads)
}


# looks at the fastq file for the read length 
length_raw_reads <- function() {
    rlength = nchar( read.table(.project$fq_files[1], nrows=2, stringsAsFactors = FALSE)[2,1] )
    return(rlength)
}

#plot_aligned_report <- function( seqdata, index, hitcount, update = FALSE ) {
plot_aligned_report <- function( update = FALSE ) {
    trace.enter("plot_aligned_report");
    on.exit({ trace.exit() })
    
    width = 480;
    height = 480;
    
    algn = "left" 
    
    report_dir = .project$report["aligned_report_dir"]
    
    filename = "plotAlignedReport"
    
    if( !file.exists(paste(report_dir, "/", filename, ".html", sep="")) || update ) {
        
        log.info( "Creating report" )
        
        HTMLInitFile( outdir=report_dir, filename=filename, Title="Quality report on aligned reads" )
        HTML( "<h1>Run info</h1>" )
        HTML( paste("Project name:", .project$name) )
        HTML( paste("Read length:", getReadLength0(.project$fq_files[[1]]) )) ## width(seqdata[[ 1 ]])[1])
        HTML( paste("Paired:", .project$pairing$type) )
        HTML( paste("Stranded:", .project[["stranded"]]) )
        HTML( paste("Quality scale:", .project$qual_type) )
        HTML( paste("Aligner:", .project$aligner$type) )
        HTML( paste("Reference:", .project$organism, .project$reference$type, 
                        .project$reference$version ) )
        
        bam_filename = paste(.project$aligner$out_dir, "/", .project$aligner$files[2], sep="" )
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 1")
        
        
        # TODO: barplot of raw reads and proportions of aligned, filtered, multi-hits, gapped reads
        # needs: bam$qname
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                Align=algn, GraphFileName="hits", 
                plotFunction=function(){
                    plot_hits( number_raw_reads(), bam_filename ); ## length(seqdata[[ 1 ]])
                });
      
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 2")
      
      
        # for PE only, insert size density 
        # needs: bam$isize
        if( .project$pairing$type == "PE" ) {
            myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                    Align=algn, GraphFileName="insertsize", 
                    plotFunction=function(){
                        plot_insertsize( bam_filename );
                    });
            gc()
        }
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 3")
        
        
        # alignment quality distribution
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                Align=algn, GraphFileName="mapquality", 
                plotFunction=function(){
                    plot_mapq( bam_filename );
                });
        gc()
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 4")
        
        
        # barplot of distribution among chromosomes
        # needs: bam$rname
        #reflen = reflen_to_dataframe()
        if( .project$reference$type == "genome" ) {
            myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                    Align=algn, GraphFileName="chrbarplot", 
                    plotFunction= function(){
                        plot_chromosome_distrib( bam_filename );    
                        
                    });
        }
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 5")
        
        
        # pie chart of distribution between strands
        # needs: bam$flag
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                Align=algn, GraphFileName="strand", 
                plotFunction=function(){
                    plot_strand_distrib( bam_filename );
                });    
        
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 6")
        
        # boxplots along transcripts
        # needs: annotation 
        if( .project$reference$type == "transcriptome" ) {
            myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
                    Align=algn, GraphFileName="alongtranscript", 
                    plotFunction=function(){
                        plot_transcript_distrib( bam_filename );
                    });    
        }
        
        # memory usage
        print.memory.usage("Memory: plot_aligned_report step 7")
        
        
        #            log.info( "To be implemented" )
        #        } else {
        #            if( .project$pairing$type == "PE" ) {
        #                aln = readGAlignments( bam_filename, param=bamParam(firstmate = TRUE ) )
        #            } else {
        #                aln = readGAlignments( bam_filename, param=bamParam() )
        #                annot = get_annotation( type="gtf", split = FALSE )
        #                subannot <- unique( annot[,c("chromosome_name", "strand", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end")] )
        #                
        #                
        #                cnames = factor( paste( "EG:", subannot$chromosome_name, sep="" ), levels(rname(aln)) )
        #                ar <- GRanges(
        #                        seqnames = cnames,
        #                        ranges = IRanges( start=subannot$exon_chrom_start, end=subannot$exon_chrom_end), 
        #                        strand = Rle( subannot$strand ), 
        #                        exon=subannot$ensembl_exon_id, 
        #                        gene=subannot$ensembl_transcript_id )
        #    
        #                uar = ar
        #                strand(uar) = Rle( rep("*", length(uar) ) )
        #                ov = countOverlaps( uar, uar ) == 1
        #            
        #                if( project[["stranded"]] )
        #                    ar = ar[ov]
        #                else 
        #                    ar = uar[ov]
        #                saln = subsetByOverlaps( aln, ar )
        #                cover <- coverage(saln)
        #                islands <- slice(cover, 1)
        #                islands = islands[["EG:4_group3"]]
        #                steps = round(width(islands) * 0.05)
        #                vals = c()
        #                ind = 0
        #                for( i in 1:20 ) {
        #                    sapply( )
        #                }
        #            }
        #        }
        
        #HTML( paste("Number of reads:", length(seqdata[[1]])) )
        #
        #myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
        #    Align=algn, plotFunction=function(){
        #    plot_bars_totalunique( bam$cigar, hitcount, seq, !index$multhit )
        #} )
        #
        #myHTMLplot( GraphDirectory=report_dir, plotFunction=function(){
        #    plot_mapq_distribution( bam$mapq )
        #} )
        
        # polyX of unaligned reads
        #HTML( "<h1>Not aligned</h1>" )
        #notaln = !(subseq(id(seq$seqdata), start=1, end=-3) %in% BStringSet(bam$qname) )
        #
        #myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height,
        #    Align=algn, GraphFileName="polyX", 
        #    plotFunction=function(){
        #        plot_polyX( bam_filename );
        #    });    
        
        HTMLEndFile()
    } else {
        log.info( "File exists, not updated" )
    }
}

plot_compared_raw_report <- function( .project, update = FALSE ) {
    trace.enter("plot_compared_raw_report");
    on.exit({ trace.exit() })
    
    width = 480;
    height = 480;
    
    algn = "left" 
    
    report_dir = .project$report["compared_report_dir"]
    
    filename = "plotComparedRawReport"
    
    if( !file.exists(paste(report_dir, "/", filename, ".html", sep="")) || update ) {
        
        HTMLInitFile( outdir=report_dir, filename=filename, Title="Compared report on raw reads" )
        
        HTML( as.title( "Comparison plots" ), append = FALSE )
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 1")
        
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height, 
                Align=algn, GraphFileName="avgbasequal_boxplots", 
                plotFunction = function(){
                    plot_avgbasequal_boxplots( update );
                })
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 2")
        
        
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height, 
                Align=algn, GraphFileName="dustyScore_distribs", 
                plotFunction = function() {
                    plot_dustyScore_distribs( update );
                })
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 3")
        
        
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height, 
                Align=algn, GraphFileName="scatter",
                plotFunction=function(){
                    plot_scatter_samples();
                })
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 4")
        
        if( .project$pairing$type == "PE" ) {    
            myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height, 
                    Align=algn, GraphFileName="insizes", 
                    plotFunction=function(){
                        plot_insize_boxplot()
                    })
        }
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 5")
        
        myHTMLplot( GraphDirectory=report_dir, Width=2*width, Height=height, 
                Align=algn, GraphFileName="hits", 
                plotFunction = function() {
                    plot_compared_hits( update );
                })
        
        # memory usage
        print.memory.usage("Memory: plot_compared_raw_report step 6")

        HTMLEndFile()
    } else {
        log.info( "File exists, not updated" )
    }
}
