#  Things that are usually not called directly but through tools.R
# 
# Author: filimon
#
##############################################################################

createDirAndBackup <- function( dirname, showWarnings = TRUE, recursive = FALSE ) {
    
    if (file.exists(dirname)) {
        
        cnt = 0;
        backup = "";
        
        while( cnt < 100 ) {
            cnt = cnt + 1;
            backup = paste(dirname, "-backup-", gsub(" ", "-", date()), "-", cnt, sep="");
            
            if (!file.exists(backup)) {
                break;
            }
        } 
        
        log.info("Backing up ", dirname, " to ", backup);
        file.rename(dirname, backup);
    }
    
    log.info("Creating ", dirname);
    
    dir.create(dirname, showWarnings = showWarnings);
}


create_feature_data <- function( gtf, tids ) {
    trace.enter("create_feature_data");
    on.exit({ trace.exit() })
    
    log.info( "Assembling exons.. " );
    
    # take all exons 
    gtf_exon = gtf[grep("exon", gtf$feature), ]
    
    log.info( "Computing transcript lengths.. " );
    
    # take unique names of transcripts
    #
    gtf_exon_unique = unique(gtf_exon$transcript_id)
    
    # create of map of transcript start indexes
    # add additionally map length as the last index
    #
    gtf_exon_unique_map = c(match(gtf_exon_unique, gtf_exon$transcript_id), length(gtf_exon$transcript_id))
    
    # compute transcript lengths, sum exon lengths
    #
    gtf_exon_unique_len = sapply(1:(length(gtf_exon_unique_map)-1), 
            function(i) { i0 = gtf_exon_unique_map[i]; i1 = gtf_exon_unique_map[i+1] - 1; 
                sum(gtf_exon$end[i0:i1] - gtf_exon$start[i0:i1] + 1); })
    
    # get lengths for expressed transcripts
    #
    tids_len = gtf_exon_unique_len[match(tids, gtf_exon_unique)];
    
    log.info( "Gathering attributes.. " );
    
    # select records that that contain protein id
    #
    gtf_prot = gtf[grep("ENS*", gtf$protein_id),]
    
    # select protein ids for expressed transcripts
    #
    tids_protein_id = gtf_prot$protein_id[match(tids, gtf_prot$transcript_id)];
    
    # select gene id, gene name, transcript name for 
    # expressed transcripts
    #
    tids_map = match(tids, gtf$transcript_id)
    
    tids_gene_id     = gtf$gene_id[tids_map];
    tids_gene_name   = gtf$gene_name[tids_map];
    tids_chr         = gtf$chr[tids_map];
    tids_transcript_name = gtf$transcript_name[tids_map];
    
    log.info( "Creating Feature Data.. " );
    
    df = data.frame(chr=tids_chr, gene_id=tids_gene_id, gene_name=tids_gene_name, 
            transcript_name=tids_transcript_name, length=tids_len, protein_id=tids_protein_id)
    
    metaData <- data.frame(labelDescription=c("Chromosome", "Gene Id", 
                    "Gene Name", "Transcript Name", "Transcript Length", "Protein Id"));
    
    fd = new("AnnotatedDataFrame", data=df, varMetadata=metaData);
    return( fd ) 
}

exprs_table_to_dataframe <- function() {
    trace.enter("exprs_table_to_dataframe");
    on.exit({ trace.exit() })
    
    tids = vector()
    ctab = list()
    i=1
    for( project in .projects ) {
        log.info( "At ", project$name ) 
        
        cfilename = paste( project$aligner$out_dir, "/", 
                project$count$method, "_", project$count$files[project$count$feature], sep="" )
        
        if( file.exists(cfilename) ) {
            
            ctab[[project$name]] = read.table( cfilename, header = TRUE, stringsAsFactors = FALSE )
            
            if( project$count$method == "mmseq" ) {
                # 
                #
                colnames(ctab[[project$name]])[c(1,3)] = c( "id", "estimate" )
            } else if( project$count$method == "cufflinks" ) {
                if( project$count$feature == "transcript" ) {
                    
                    # fetch transcript info
                    #
                    
                    # new workflow
                    #
                    colnames(ctab[[project$name]])[c(1,8,9,10)] = c( "id", "length", "coverage", "estimate" )
                    
                } else {
                    
                    # fetch genes info
                    #
                    
                    # new workflow
                    #
                    colnames(ctab[[project$name]])[c(1,8,9,10)] = c( "id", "length", "coverage", "estimate" )
                }
                    
            } else {
                colnames(ctab[[project$name]]) = c( "id", "estimate", "length" )
            }
            tids = unique( c(tids, ctab[[project$name]][,"id"]) )
        } else {
            registerExpStepWarning(ASSEMBLING_ESET, "Estimation data missing, check logs");
            
            log.warning(project$name, " warning: ", cfilename, " not found")
            
            ctab[[project$name]] = NULL
        }
        i = i+1
    }
    
    tids = sort(tids)
    data = matrix( 0, nrow=length(tids), ncol=length(.projects), 
            dimnames=list(feature_id=tids,sample=names(.projects)) )
    if( "length" %in% colnames(ctab[[1]]) ) {
        lengths = matrix( NA, nrow=length(tids), ncol=1, 
                dimnames=list(feature_id=tids,"length") )
        for( pn in names(.projects) ) {
            lengths[match(ctab[[pn]][,"id"],rownames(data)),1] = ctab[[pn]][,"length"]
        }
    } else 
        lengths = NULL
    
    for( pn in names(.projects) ) {
        log.info( "Reading ", pn )
        
        if (!is.null(ctab[[pn]])) {
            data[match(ctab[[pn]][,"id"],rownames(data)),pn] = ctab[[pn]][,"estimate"]
#            if( type == "fpkm" ) {
#                data[match(ctab[[pn]][,"trans_id"],rownames(data)),pn] = ctab[[pn]][,"FPKM"]
#            } else if( type == "count" ) {
#                data[match(ctab[[pn]][,"trans_id"],rownames(data)),pn] = ctab[[pn]][,"coverage"] * ctab[[pn]][,"length"]
#            }
        } else {
            data[,pn] = NA;
        }
    }
    
    return(list(data=data,tids=tids,length=lengths))
}

verify_fqfiles <- function( files ) {
    trace.enter("verify_fqfiles");
    on.exit({ trace.exit() })
    
    if( length(files) == 0 ) {
        log.info( "ERROR No FastQ files found! " );
        setPackageVariable("missing-fastq" = TRUE);
        #stop();
        #return( NULL )
    }
    fl = files[grep( "_1\\.", files )]
    fl = c( fl, files[grep( "_2\\.", files )] )
    
    if( length(fl) == 0 ) {
        fl = files
    }
    
    for(f in fl) {
        log.info("FASTQ: ", f);
    }
    
    return( fl )
}

fasta_to_list <- function() {
    trace.enter("fasta_to_list");
    on.exit({ trace.exit() })
    
    ref = readDNAStringSet( paste(.project$reference$file, ".fa", sep=""), format="fasta" );
    
    names(ref) = sapply( ref, function(x)  { i = regexpr("ENS:", x$desc)[1]; substr( x$desc, i, i+18) } )
    return( ref )
}

reftype_to_filename <- function( type ) {
    trace.enter("reftype_to_filename");
    on.exit({ trace.exit() })
    
    if( type == "genome" || type == "dna" ) {
        reference = "dna.chromosome"
    } else if( type == "transcriptome" ) {
        reference = "cdna.chromosome"
    } else {
        log.info( "Unknown reference type" );
        return( NULL )
    }
    return( reference )
}

sam_addXS <- function() {
    trace.enter("sam_addXS");
    on.exit({ trace.exit() })
    
    infile = paste( .project$aligner$out_dir, "/", .project$aligner$files[1], sep="" )
    outfile = sub( ".sam", ".XS.sam", infile )
    res = vector( mode="integer", length=1 )
    result = .C( "addXS", as.character(infile), as.character(outfile), as.integer(res), PACKAGE = "ArrayExpressHTS" )[[3]]
}

# receives a DNAStringSet or a character vector
polyX_size <- function( reads, update = FALSE ) {
    trace.enter("polyX_size");
    on.exit({ trace.exit() })
    
    filename = .project$aux_files["polyCount"]
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." );
        log.info( "Loading..." );
        load( filename )
    } else {
        count = list()
        bases = c("A", "T", "C", "G")
        reads = as.character(reads)
        for( b in bases ) {
            log.info( "Counting ", b, "s..." )
            count[[b]] = vector( mode="integer", length=length(reads) )
            count[[b]] = .Call( "count_polyL", as.character(b), as.integer(length(reads)), as.character(reads), PACKAGE = "ArrayExpressHTS" )
        }
        save( count, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    return( count )
}


# make useful indexes of unique and not aligned reads in aln and seq objects
make_indexes_old <- function( seq, aln, alntab, update = FALSE, return = TRUE ) {
    trace.enter("make_indexes_old");
    on.exit({ trace.exit() })
    
    
    filename = .project$aux_files["indexes"]
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        index = list()
        index[["seq"]] = list()
        index[["aln"]] = list()
        
        if( .project$aligner$type == "tophat" )
            index[["aln"]][["um"]] = (as.character(id(aln)) %in% names(alntab$top)[alntab$top == 1])
        else if( .project[["aligner"]][["type"]] == "bwa" )
            index[["aln"]][["um"]] = (as.character(id(aln)) %in% names(alntab$top)[alntab$top == 1])
        else
            log.info( "I don't know how to calculate indexes for this type of aligner..." )
        
        index[["seq"]][["notaln"]] = !(subseq(id(seq$seqdata), start=1, end=-3) %in% id(aln))
        
        save( index, file=filename )
        log.info( "Saved to ", filename, "" )
    }    
    
    if( return )
        return( index )
    return( NULL )
}

# transform the names into the representation wished for the project (e.g. from '1' to 'chr1' )
get_chr_representation <- function( chrnames ) {
    trace.enter("get_chr_representation");
    on.exit({ trace.exit() })
    
    
    # if there is at least one name that does not contain the project defined representation string then change
    if( length( grep( paste( "^", .project$reference$chrRepresentation, sep="" ), chrnames) ) != length(chrnames) ) {
        chrnames = sapply( chrnames, function(x) paste( .project$reference$chrRepresentation, 
                            ifelse( x != "MT", x, "M" ), sep="" ) )
        rownames( chrnames ) = NULL
        cnames = unique( chrnames )
        log.info( "Changed chromosome names to: ", paste(cnames, collapse=", "), "" )
    } else {
        log.info( "Nothing changed" )
    }
    return( as.character(chrnames) )
}

# filter unkown contigs (NT_*, GL and HS correspond to mmusculus and human)
filter_chr <- function( chrnames ) {
    trace.enter("filter_chr");
    on.exit({ trace.exit() })
    
    return( chrnames[grep("^NT_|GL|HS|Unknown_", chrnames, invert = TRUE)] )
}

# transform a decimal number into the corresponding 8-bit binary number (to get the bitwise flag from SAM)
flag_meaning <- function( dec ) {
    trace.enter("flag_meaning");
    on.exit({ trace.exit() })
    
    
    meaning = c( "paired", "in proper pair", "unmapped", "mate unmapped", "reverse strand", "mate on reverse strand", 
            "first of pair", "second of pair", "not primary alignement", "failed quality", "artifact" )
    unmeaning = c( "not paired", "not in proper pair", "mapped", "mate mapped", "forward strand", "mate on forward strand", 
            "second of pair", "first of pair", "primary alignement", "passed quality", "not artifact" )
    mapped = TRUE
    
    bin = flagint_to_bit( dec )
    
    log.info( "The read is:" )
    for( i in 1:length(bin) ) {
        if( rev(bin)[i] == 1 ) {
            log.info( "\t",meaning[i],"" )
            if( i == 3 )
                mapped = FALSE
        } else {
            log.info( "\t", unmeaning[i], "" )
        }
    }
    log.info( bin )
    return( mapped )
}

flagint_to_bit <- function( dec ) {
    trace.enter("flagint_to_bit");
    on.exit({ trace.exit() })
    
    bin = vector()
    
    while( dec > 0 ) {
        bin = c( dec%%2, bin )
        dec = floor( dec/2 )
    }
    
    zeros = 11 - length(bin)
    for( i in 1:zeros )
        bin = c( 0, bin )
    
    return( bin )
}

# before using this, which is very slow, consider using scanBamFlag() when using scanBam() (see bam_to_list())
flag_is_mapped_and_proper <- function( decs ) {
    trace.enter("flag_is_mapped_and_proper");
    on.exit({ trace.exit() })
    
    k = 1
    
    proper = vector( mode="logical", length=length(decs) )
    for( dec in decs ) {
        progress( k*100/length(decs) )
        bin = flagint_to_bit( dec )
        if( (bin[3] | bin[4]) == 0 )
            proper[k] = TRUE
        k = k+1
    }
    
    return(proper)
}

flag_is_mapped <- function( decs ) {
    trace.enter("flag_is_mapped");
    on.exit({ trace.exit() })
    
    k = 1
    mapped = vector()
    for( dec in decs ) {
        progress( k )
        mapped[k] = TRUE
        bin = vector()
        
        while( dec > 0 ) {
            bin = c( dec%%2, bin )
            dec = floor( dec/2 )
        }
        
        zeros = 11 - length(bin)
        for( i in 1:zeros )
            bin = c( 0, bin )
        
        for( i in 1:length(bin) ) {
            if( i == 3 && rev(bin)[i] == 1 )
                mapped[k] = FALSE
        }
        k = k + 1
    }
    return( mapped )
}

flag_is_firstmate <- function( flags ) {
    trace.enter("flag_is_firstmate");
    on.exit({ trace.exit() })
    
    #res = sapply( decs, function(x) { bin = binary(x)$binary; any(c(bin, rep(0,11-length(bin))) & c(0,0,0,0,0,0,1,0,0,0,1)) } )
    #return(res)
    .Call( "is_firstmate", as.integer(length(flags)), as.integer(flags), PACKAGE = "ArrayExpressHTS" )
}

flag_is_secondmate <- function( flags ) {
    trace.enter("flag_is_secondmate");
    on.exit({ trace.exit() })
    
    .Call( "is_secondmate", as.integer(length(flags)), as.integer(flags), PACKAGE = "ArrayExpressHTS" )
}

ASCIIconvert <- function(i, type=c("CtoI", "ItoC") ) {
    #trace.enter("ASCIIconvert");
    #on.exit({ trace.exit() })
    
    
    ASCII <- c( "","\001","\002","\003","\004","\005","\006","\007", # 000-007
            "\010","\011","\012","\013","\014","\015","\016","\017", # 010-017
            "\020","\021","\022","\023","\024","\025","\026","\027", # 020-027
            "\030","\031","\032","\033","\034","\035","\036","\037", # 030-037
            "\040","\041","\042","\043","\044","\045","\046","\047", # 040-047
            "\050","\051","\052","\053","\054","\055","\056","\057", # 050-057
            "\060","\061","\062","\063","\064","\065","\066","\067", # 060-067
            "\070","\071","\072","\073","\074","\075","\076","\077", # 070-077
            "\100","\101","\102","\103","\104","\105","\106","\107", # 100-107
            "\110","\111","\112","\113","\114","\115","\116","\117", # 110-117
            "\120","\121","\122","\123","\124","\125","\126","\127", # 120-127
            "\130","\131","\132","\133","\134","\135","\136","\137", # 130-137
            "\140","\141","\142","\143","\144","\145","\146","\147", # 140-147
            "\150","\151","\152","\153","\154","\155","\156","\157", # 150-157
            "\160","\161","\162","\163","\164","\165","\166","\167", # 160-167
            "\170","\171","\172","\173","\174","\175","\176","\177", # 170-177
            "\200","\201","\202","\203","\204","\205","\206","\207", # 200-207
            "\210","\211","\212","\213","\214","\215","\216","\217", # 210-217
            "\220","\221","\222","\223","\224","\225","\226","\227", # 220-227
            "\230","\231","\232","\233","\234","\235","\236","\237", # 230-237
            "\240","\241","\242","\243","\244","\245","\246","\247", # 240-247
            "\250","\251","\252","\253","\254","\255","\256","\257", # 250-257
            "\260","\261","\262","\263","\264","\265","\266","\267", # 260-267
            "\270","\271","\272","\273","\274","\275","\276","\277", # 270-277
            "\300","\301","\302","\303","\304","\305","\306","\307", # 300-307
            "\310","\311","\312","\313","\314","\315","\316","\317", # 310-317
            "\320","\321","\322","\323","\324","\325","\326","\327", # 320-327
            "\330","\331","\332","\333","\334","\335","\336","\337", # 330-337
            "\340","\341","\342","\343","\344","\345","\346","\347", # 340-347
            "\350","\351","\352","\353","\354","\355","\356","\357", # 350-357
            "\360","\361","\362","\363","\364","\365","\366","\367", # 360-367
            "\370","\371","\372","\373","\374","\375","\376","\377" ) # 370-377
    
    if( type == "ItoC" )
        return( ASCII[i %% 256 + 1] )
    else if( type == "CtoI" ) 
        return(match(i, ASCII) - 1)
}        

# to move files with the wrong name to files with names that can be found 
# in the project. e.g. usage: move_files( "d_melanogaster", "dmelanogaster", "." )
move_files <- function( pattern, replacement, dir, sure.i.want.to.run.this = FALSE ) {
    trace.enter("move_files");
    on.exit({ trace.exit() })
    
    files = dir( dir, pattern=pattern, full.names = TRUE )
    mod_files = sub( pattern, replacement, files ) 
    for( i in 1:length(files) ) {
        cmd = paste( "mv", files[i], mod_files[i] )
        log.info( cmd, "" )
        if( sure.i.want.to.run.this )
            system( cmd )
    }
}

# get back a dusty score to help identify low complexity reads
# for the meaning of the Dusty score see: https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2009-February/000170.html
stringset_to_dustyscore <- function( strings ) {
    trace.enter("stringset_to_dustyscore");
    on.exit({ trace.exit() })
    
    return( dustyScore(strings) )
}

# reads a fastq file into a data.frame
fastq_to_dataframe <- function( project, update = FALSE, return = TRUE ) {
    trace.enter("fastq_to_dataframe");
    on.exit({ trace.exit() })
    
    filename = project$aux_files["fqAsDataFrame"]
    
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists." )
        if( return ) {
            log.info( "Loading... " )
            varnames = load( filename )
            log.info( "loaded ", paste(varnames,collpase=", "), "" )
        } else 
            return( NULL )
    } else {
        log.info( "Reading file " )
        i = 1
        for( file in project$fq_files ) {
            log.info( file, "" )
            lines = read.table( file, stringsAsFactors = FALSE, comment.char="" )
            len = dim(lines)[1]
            if( i == 1 ) {
                seqframe = data.frame( id=lines[seq(1,len,4),1], seq=lines[seq(2,len,4),1], strand=lines[seq(3,len,4),1],
                        quality=lines[seq(4,len,4),1], stringsAsFactors = FALSE )
            } else {
                seqframe = rbind( seqframe, data.frame( id=lines[seq(1,len,4),1], seq=lines[seq(2,len,4),1], strand=lines[seq(3,len,4),1],
                                quality=lines[seq(4,len,4),1], stringsAsFactors = FALSE ) )
                log.info( "rbinding..." )
            }
            i = i + 1
        }
        log.info( "done" )
        save( seqframe, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    return( seqframe )
}

fastqquality_to_integer <- function( qual ) {
    trace.enter("fastqquality_to_integer");
    on.exit({ trace.exit() })
    
    return( as(qual,"numeric") )
}

filter_shortread_by_sequence <- function( data, seq ) {    
    trace.enter("filter_shortread_by_sequence");
    on.exit({ trace.exit() })
    
    seqfilter <- srFilter( function(x) { grep(seq, sread(x)) }, name="seqfilter" )
    return( data[seqfilter(data)] )
}

# receives a list of "PhredQuality"s  
# should only be applied on qualities from aligned reads as it is saved on the aligner output dir
phred_to_avgqual <- function( quals, update = FALSE ) {
    trace.enter("phred_to_avgqual");
    on.exit({ trace.exit() })
    
    
    filename = .project$aux_files["avgQuals"]
    
    if( file.exists(filename) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading..." )
        load( filename )
    } else {
        # in pure R it would be like this, very very slow:
        # avgquals = sapply( as.character(quals), function(x) mean((as.integer(charToRaw(x))-33)) )
        # or even one could use a filter
        # qualfilter <- srFilter( function(x) { (alphabetScore(quality(x)) / width(quality(x))) < qual })
        # filtdata = data[qualfilter(data)]
        
        avgquals = list()
        for( i in 1:length(quals) ) {
            qual = quality(quality(quals[[i]]))
            
            # instead I call a C function which is much faster
            avgquals[[i]] = vector( mode="double", length=length(qual) )
            avgquals[[i]] = .Call( "phred_to_average_qual", as.integer(length(qual)), as.character(qual), PACKAGE = "ArrayExpressHTS" )
        }
        save( avgquals, file=filename )
        log.info( "Saved to ", filename, "" )
    }
    
    return( avgquals )
} 

gtf_to_dataframe <- function( filename, names=c( "gene_id", "transcript_id", "exon_number", "gene_name", "transcript_name", "protein_id" ), filter = TRUE ) {
    trace.enter("gtf_to_dataframe");
    on.exit({ trace.exit() })
    
    log.info( "loading gtf.." );
    
    #gtf = read.table( filename[1], fill = TRUE,
    #       col.names = c( "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes" ),
    #        sep="\t", strip.white = TRUE, stringsAsFactors = FALSE )
    gtf = read.table( filename[1], fill = TRUE, stringsAsFactors = FALSE )
    
    chrs = unique( gtf$chr )
    
    if (filter) {
        log.info( "filtering gtf.." );
        
        chrs = filter_chr( chrs )
        gtf = gtf[gtf$chr %in% chrs,]
    }
    
    log.info( "removing unneeded columns.." );
    
    # remove the score column
    gtf = gtf[,-6]
    
    log.info( "preparing attribute data.." );
    
    splitat = strsplit(gtf$attributes, ";")
    attribs = list()
    
    log.info( "allocating columns.." );
    
    splitatcnt = length(splitat);
    
    for( name in names ) {
        attribs[[name]] = rep( NA, splitatcnt )
    }
    
    log.info( "parsing attributes.." );
    
    count = 1
    
    perc = floor(splitatcnt / 20);
    curp = perc;
    
    for( line in splitat ) {

        pairs = strsplit(sub("^ ", "", line), " ")
        for( pair in pairs )
            attribs[[pair[1]]][count] = pair[2]
        
        if (count == perc) {
            log.info(curp / perc, "%% ",count," out of ",splitatcnt," complete.." );
            curp = curp + perc;
        }
        
        count = count + 1
    }
    
    # remove the attributes column
    gtf = gtf[,-8]
    bnames = names(gtf)
    
    log.info( "attaching columns.." );
    
    for( at in attribs ) {
        gtf = cbind( gtf, at )
    }
    names(gtf) = c(bnames, names)
    
    log.info( "done" );
    
    return(gtf)
    
    #perc = floor(splitatcnt / 20);
    #curp = perc;
    
    
    #for( line in splitat ) 
    #for( i in 1:splitatcnt ) {
    #    #pairs = strsplit(sub("^ ", "", line), " ")
    #    
    #    pairs = strsplit(line, " ")
    #    
    #    #for( pair in pairs )
    #    #attribs[[pair[1]]][count] = pair[2]
    #    if (i == curp) {
    #        log.info(ceiling(curp * 100 / splitatcnt), "%% ",count," out of ",splitatcnt," complete.." );
    #        curp = curp + perc;
    #    }
    #    
    #}
    
    
}

run_cmds <- function( cmds, run ) {
    trace.enter("run_cmds");
    on.exit({ trace.exit() })
    
    if( !run )
        log.info( "These commands were not run (you can copy&paste on the console instead):" )
    else 
        log.info( "Running:" )
    for( cmd in cmds ) {
        
        log.info( "\t", cmd, "" )
        
        if( run ){
            fail = system( cmd )
            if( fail ) {
                log.info( cmd, " Failed! Stopping..." )
                return( fail )
            }
        }
    }
    log.info( "" )
    return(0)
}

# returns a list of sizes per chromosome/transcript in the reference fasta
reflen_to_dataframe <- function( type=.project$reference$type ) {
    trace.enter("reflen_to_dataframe");
    on.exit({ trace.exit() })
    
    reference = reftype_to_filename( type )
    filename = paste( .project$aligner$indexes_dir, "/", .project$organism, ".", 
            .project$reference$version, ".", reference, ".fa", sep="" )
    chrinfo = fasta.seqlengths( filename )
    return( chrinfo )
}

#bam_to_alignedread <- function( project, update = FALSE, return = TRUE ) {
#    trace.enter("bam_to_alignedread");
#    on.exit({ trace.exit() })
#    filename = project[["aux_files"]]["alignedRead"]
#    if( file.exists(filename) && !update ) {
#        log.info( "File", filename, "exists." )
#        if( return ) {
#            log.info( "Loading..." )
#            load( filename )
#        } else 
#            aln = NULL
#    } else {
#        aln = .readAligned_bam( project$aligner_out_dir, pattern=paste( project$aligner_files[2], "$", sep="") )
#        save( aln, file=filename )
#        log.info( "Saved to file", filename, "" )
#    }
#    
#    return( aln )
#}

