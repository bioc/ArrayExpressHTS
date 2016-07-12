readIDF <- function( idffname ) {
    trace.enter("readIDF");
    on.exit({ trace.exit() })
    
    idffile = scan(idffname, character(), sep="\n")
    
    idf.data = list()
    
    for(line in idffile) { 
        strings = unlist(strsplit(line,"\t"))
        
        key = strings[1];
        key = gsub('"', '', key)
        
        if (length(strings) > 1) {
            values = strings[2:length(strings)]
        } else { 
            values=NA
        }
        
        idf.data[[key]]=values
    }
    
    return (idf.data);
}

attachExperimentData = function(eset, idf.data) {
    trace.enter("attachExperimentData");
    on.exit({ trace.exit() })
    
    ## making key matches #
    ## (Person Last Name, Person First Name, Person Mid Initials), Person Email, Person Phone, Person Address, Person Affiliation, 
    Person_Name = c(idf.data$"Person First Name", idf.data$"Person Last Name",idf.data$"Person Mid Initials");
    Personal_contact=c(idf.data$"Person Email", idf.data$"Person Phone", idf.data$"Person Address");
    
    ## making experimentData object #        
    SubmitterIndex=which(idf.data$"Person Roles"=="submitter");
    
    experimentData = new("MIAME", 
            name = as.character(paste(idf.data$"Person Last Name"[SubmitterIndex],", ",idf.data$"Person First Name"[SubmitterIndex], sep="")), #performer
            lab = as.character(idf.data$"Person Affiliation"[SubmitterIndex]) , #Person Affiliation 
            contact = as.character(idf.data$"Person Email"[SubmitterIndex]), # Person Email(Person Phone, Person Address)
            title = as.character(idf.data$"Investigation Title") , #description #Investigation Title
    
            ##abstract= "",    #not provided in the idf.data
            ##url    = "",
    
            other = list(
                    accession = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
                    identifier = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
                    ##Experimental Factor Type
                    experimentalFactor = c(idf.data$"Experimental Factor Type"), 
                    ##Experimental Design
                    type = c(idf.data$"Experimental Design"),
                    measurementType = experimentData(eset)@other$measurementType #from processed data.zip depending on user answer about QT type
                    )
            );
    
    experimentData(eset) = experimentData;
    
    return(eset);
}


read_gtf <- function( filename ) {
    
    if( file.exists(filename) ) {
        
        log.info( "Reading processed annotation ", filename )
        
        # column names of processed annotation file
        # 
        col.names = c( "chr", "source", "feature", "start", "end", "score", "strand", "phase", "gene_id", 
            "transcript_id", "exon_number", "gene_name", "transcript_name", "protein_id", "gene_biotype"  )
        
        gtf = read.table( filename, fill = TRUE,
            col.names = col.names, sep="\t", strip.white = TRUE, stringsAsFactors = FALSE )
        
        return (gtf);
    
    } else {
        log.info( "Could not find processed annotation file ", filename )
    }
    
    return (NULL);
}

standardise <- function( e, by_length = TRUE, by_libsize = TRUE ) {
    trace.enter("standardise");
    on.exit({ trace.exit() })
    
    
}

# the functions in this file assume that .projects has been set
to_expressionset <- function( update = FALSE, filter = TRUE ) {
    trace.enter("to_expressionset");
    on.exit({ trace.exit() })
    
    if(.projects[[1]]$count$standardise)
        name = "std"
    else 
        name = "notstd"
    
    if(.projects[[1]]$count$normalisation == "none")
        name = paste( name, "_notnorm", sep="" )
    else 
        name = paste( name, "_", .projects[[1]]$count$normalisation, sep="" )
    
    filename = paste( .projects[[1]]$basedir, "/eset_", name, ".RData", sep="" )
    
    if( file.exists( filename ) && !update ) {
        log.info( "File ", filename, " exists." )
        log.info( "Loading... " )
        varnames = load( filename )
        return( get(varnames[[1]]) );      
    } else {        
        res = exprs_table_to_dataframe()
        data = res$data
        tids = res$tids
        length = res$length 
        eset = new("ExpressionSet", exprs = data)
        
        # attach experiment data
        #
        idffname = dir( .projects[[1]]$datadir, ".idf", full.names = TRUE )[1]
        
        if (file.exists(idffname)) {
            eset = attachExperimentData(eset, readIDF(idffname));
            #experimentData(eset) = readIDF(idffname);
        }
        
        # attach phenoData
        #
        sdrffname = dir( .projects[[1]]$datadir, ".sdrf", full.names = TRUE )[1]
        
        if (file.exists(sdrffname)) {
            sdrf = readSDRF(sdrffname);
            
            resultset = attr(.projects, 'metadata')$resultset;
            if (!is.null(resultset)) {
                attr(sdrf, 'resultset') = resultset;
            }
           
            phenoData(eset) = sdrf;
        }
        
        # attach .projects
        # reference organism and version are 
        # stored in projects per each project
        #
        attr(eset, 'projects') = .projects;
        
        # read annotation
        #
        
        gtf = attr(.projects, 'metadata')$gtf;
        
        if (is.null(gtf)) {
            log.info( "Reading annotation.. " );
            
            gtf = read_gtf( paste( .projects[[1]]$annot$folder, "/", 
                .projects[[1]]$organism, ".", .projects[[1]]$reference$version, ".gtf.out", sep="" ) );
            
            attr(.projects, 'metadata')$gtf = gtf;
        }
        
        # get the transcript lengths
        if( !is.null(length) ) {
            log.info( "Length information available from the count method... " );
            df = data.frame(transcript_name=tids, length=length )
            metaData <- data.frame(labelDescription=c("Transcript Name", 
                            "Transcript Length"));
        
            featureData(eset) = new("AnnotatedDataFrame", data=df, varMetadata=metaData);
        } else if (!is.null(gtf)) {
            log.info( "Creating feature data from gtf.. " );
            featureData(eset) = create_feature_data( gtf, tids )
        } else {
            log.info( "Creating feature data from cDNA file... " );
            fa = .projects[[1]]$reference$file
            fa = paste(sub( ".dna.chromosome", ".cdna.chromosome", fa, fixed = TRUE), ".fa", sep="" )
            lengths = fasta.seqlengths(fa)
            pos = regexpr( " ", names(lengths)[1] )
            names(lengths) = substr(names(lengths),1,pos)
            
            df = data.frame(transcript_name=tids, length=lengths[match(tids,names(lengths))] )
            metaData <- data.frame(labelDescription=c("Transcript Name", 
                            "Transcript Length"));
            featureData(eset) = new("AnnotatedDataFrame", data=df, varMetadata=metaData);
            
        }

        if( .projects[[1]]$count$method == "cufflinks" ) { # TODO solve and remove this error
            log.info( "*** Please note that cufflinks estimates are returned as FPKMs ***" );
        } else {
            if( .projects[[1]]$count$method == "mmseq" ) {
                if( .projects[[1]]$count$standardise || .projects[[1]]$count$normalisation != "none" ) {
                    eset = un_std(eset)
                }
            }
            call_normalisator( eset, run = TRUE )
            if( .projects[[1]]$count$standardise ) {
                eset = std(eset)
            }
        }
        save( eset, file=filename );
        
        log.info( filename, " saved" );
    
        return(eset);
    }
}

std <- function( e ) {
    trace.enter("std");
    on.exit({ trace.exit() })
    
    areads = number_aligned_reads()
    new_e = exprs(e)*(1e3*1e6)/(fData(e)$length*areads)
    exprs(e) = new_e
    
    return(e)
}

# remove standardization
un_std <- function( e ) {
    trace.enter("un_std");
    on.exit({ trace.exit() })
    
    areads = number_aligned_reads()
    new_e = (exprs(e)*fData(e)$length*areads)/(1e3*1e6)
    exprs(e) = new_e
    
    return(e)
}

to_dataframe <- function() {
    trace.enter("to_dataframe");
    on.exit({ trace.exit() })
    
    ctab = list()
    for( project in .projects ) {
        filename = paste( project$aligner$out_dir, "/", project$count$files[project$count$feature], sep="" )
        ctab[[project$name]] = read.table( filename, header = TRUE, stringsAsFactors = FALSE )[,c("trans_id", "FPKM")]
        tids = unique( c(tids, ctab[[project$name]][,"trans_id"]))
    }
    tids = sort(tids)
    data = matrix( 0, nrow=length(tids), ncol=length(.projects), dimnames=list(feature_id=tids,sample=names(.projects)) )
    for( pn in names(.projects) ) {
        data[match(ctab[[pn]][,"trans_id"],rownames(data)),pn] = ctab[[pn]][,"FPKM"]
    }
    return( data )
}

get_DE <- function( e, factor ) {
    trace.enter("get_DE");
    on.exit({ trace.exit() })
    
    res = list()
    cb = combn( unique(pData(e)[[factor]]), 2 )
    conds <- c( pData(e)[[factor]] )
    for( i in 1:ncol(cb) ) {
        cds <- newCountDataSet( round(exprs(e)), conds, phenoData=phenoData(e), featureData=featureData(e) )
        cds <- estimateSizeFactors( cds )
        # if there is only a monoplicate :) for each condition, use pooling
        if( sum(conds == cb[1,i]) + sum(conds == cb[2,i]) == 2 ) {
            cds <- estimateVarianceFunctions( cds, pool = TRUE )
        } else {
            cds <- estimateVarianceFunctions( cds )
        }
        res[[paste(cb[1,i], cb[2,i], sep="x")]] <- nbinomTest( cds, cb[1,i], cb[2,i] )
    }
    return( res )
}

# table of the number of reads aligned to genome, transcriptome and their overlap
table_transcriptome_vs_genome <- function( qnames, nreads ) {
    trace.enter("table_transcriptome_vs_genome");
    on.exit({ trace.exit() })
    
    #nreads = 22881837
    df = data.frame(row.names=names(qnames))
    overlaps = list()
    for( i in 1:length(qnames) ) {
        qn1 = qnames[[i]][[1]]
        qn2 = qnames[[i]][[2]]
        log.info( "Calculating overlaps" )
        overlaps[[1]] = qn1 %in% qn2
        overlaps[[2]] = qn2 %in% qn1
        df[i,1] = round(length(unique(qn1)) / nreads * 100)
        df[i,2] = round(length(unique(qn2)) / nreads * 100)
        df[i,3] = round(sum(overlaps[[1]]) / nreads * 100)
        df[i,4] = round(sum(overlaps[[1]]) / length(unique(qn2)) * 100 )
        df[i,5] = round(sum(overlaps[[2]]) / length(unique(qn1)) * 100 )
    }
    colnames( df ) = c( "T (% of total reads)", "G (% of total reads)", "T&G (% of the total)", "T|G", "G|T" )
    return( df )
}

# compares transcript fasta files in terms of the biggest set
table_transcript_overlap <- function() {
    trace.enter("table_transcript_overlap");
    on.exit({ trace.exit() })
    
    fasta_files = paste( sapply( .projects, function(x) x$reference$file ), ".fa", sep="" )
    
    lengths = lapply( fasta_files, function(x) fasta.seqlengths(x) )
    names( lengths ) = sapply( .projects, function(x) x$organism )
    
    tnames = lapply( lengths, function(x) names(x) )
    names( tnames ) = names( lengths )
    
    transtab = list()
    for( i in 1:length(names(tnames)) ) {
        # ad hoc solution for David's transcript names
        if( length(grep("\\|", tnames[[i]][1])) ) {
            tnames[[i]] = strsplit( tnames[[i]], "\\|" )
            transtab[[i]] = data.frame( id=sapply(tnames[[i]], function(x) x[4]), 
                    chr=sapply(tnames[[i]], function(x) x[1]), 
                    start=sapply(tnames[[i]], function(x) x[2]),
                    end=sapply(tnames[[i]], function(x) x[3]),
                    strand=rep(NA, length(tnames[[i]])),
                    gene=rep(NA, length(tnames[[i]])) )
        } else { # this should be used by default for cdna files downloaded from ENSEMBL
            tnames[[i]] = strsplit( tnames[[i]], " " )
            
            chrs = strsplit( sapply(tnames[[i]], function(x) x[3] ), ":" )
            transtab[[i]] = data.frame( id=sapply(tnames[[i]], function(x) x[1]), 
                    chr=sapply(chrs, function(x) x[3]), 
                    start=sapply(chrs, function(x) x[4]),
                    end=sapply(chrs, function(x) x[5]),
                    strand=sapply(chrs, function(x) x[6]),
                    gene=sapply(tnames[[i]], function(x) strsplit(x[4], ":")[[1]][2]) )
        }
    }
    names(transtab) = names( tnames )
    
    cbn = combn( names(transtab), 2 )
    ov = data.frame( row.names="Transcripts")
    for( i in 1:ncol(cbn) ) {
        # in terms of the biggest set!
        ov[1,i] = sum(transtab[[cbn[1,i]]]$id %in% transtab[[cbn[2,i]]]$id) / 
                max(c(length(transtab[[cbn[1,i]]]$id),length(transtab[[cbn[2,i]]]$id))) *100
        colnames(ov)[i] = paste( cbn[,i], collapse=" vs " )
    }
    
    return(ov)
}
