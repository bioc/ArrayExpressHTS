#
# Everything to do with preparation of Reference and Annotation
#
#


fixReferenceDoubleY <- function( reffname ) {
    
    log.info("Fixing Reference ", reffname);
    
    # check if file is there
    #
    if (!file.exists(reffname)) {
        log.info("Reference Not Found ", reffname);
        return();
    }
    
    lines001 = system(paste("grep -rn '>Y' ", reffname, sep=""), intern=TRUE);
    
    # check if fix is required
    #
    
    if (length(lines001) > 1) {
        
        # backup
        #
        log.info("Backing Up To ", reffname, ".backup" );
        
        system(paste("cp ", reffname, " ", reffname, ".backup", sep=""));
        
        # fix the reference
        #
        
        index001 = index001 = regexpr(":", lines001)
        
        start1 = substr(lines001[1],1, index001[1]-1);
        start2 = substr(lines001[2],1, index001[2]-1);
        
        log.info("Removing Double Y Range ", start1, " - ", as.integer(start2) - 1);
        
        cmd001 = paste("sed -i ", start1, ",", as.integer(start2) - 1, "d ", reffname, sep="");
        
        system(cmd001);
        
        log.info("Reference Fixed");
        
    } else {
        
        log.info("Reference is OK");
    }
}


getMatchedText <- function( pattern, string ) {
    trace.enter("getMatchedText");
    on.exit({ trace.exit() })
    
    #match = regexpr("\\w*.\\w*", filename); 
    match = regexpr(pattern, string); 
    substr(string, match, match + attr(match, "match.length") - 1);
}

getEnsemblReference <- function( organism, type, version, location, run = FALSE, refresh = FALSE ) {
    trace.enter("getEnsemblReference");
    on.exit({ trace.exit() })
    
    #owd = getwd()
    cmds = c()
    if( type == "genome" || type == "dna" ) {
        ty = "dna"
    } else { 
        ty = "cdna"
    }
    
    reference_dir = paste(location, "/reference_genomes", sep="");
    
    if( version != "current" ) {
        release_ver = getMatchedText("\\d*$", version);
        ref_subver = 
        # pub/release-62/fasta/
        ref_subverpath = paste( "release-", release_ver, "/fasta/", sep="")
        
    } else {
        #"ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"        
        ref_subverpath = paste( "current_fasta/", sep="")
    }
    
    ensembl_path = paste("ftp://ftp.ensembl.org/pub/", 
            ref_subverpath, tolower(organism), "/", ty, sep="");
    
    # create dir for all references
    if( !file.exists(reference_dir) ) {
        cmds = c( cmds, paste( "mkdir ", reference_dir, sep="" ) )
    }
    
    dirlist = paste(organism, ".dirlist.txt", sep="");
    
    if( type == "genome" ) {
        cmds = c( cmds, paste( "curl -s -l ", ensembl_path, "/ | grep toplevel.fa > ", dirlist, sep="" ) )
    } else {
        cmds = c( cmds, paste( "curl -s -l ", ensembl_path, "/ | grep all.fa > ", dirlist, sep="" ) )
    }
    
    # run commands:
    # list cdna/dna files in the cdna/dna folder for the organism
    #
    
    run_cmds( cmds, run )
    
    cmds = c()
    
    ftpfile = read.table( dirlist, stringsAsFactors = FALSE )[1,1]
    
    version = substr( ftpfile, nchar(organism)+2, regexpr(ty,ftpfile)-2 )
    
    organism_ver = paste(organism, ".", version, sep="")
    
    filename = paste( organism_ver, ".", ty,".chromosome.fa", sep="" )
    
    log.info( "Looking for ", filename, "..." )
    
    # check if the file exists
    #
    
    if( file.exists(paste( reference_dir, "/", organism_ver, "/", filename, sep="") ) && !refresh ) {
        # do nothing
        #
        
        log.info( "File ", filename, " exists." )
    } else {
        # if the file is missing or refresh is 
        # requested, download the file and unpack it
        
        # download reference
        cmds = c( cmds, paste( "curl -O ", ensembl_path, "/", ftpfile, sep="" ) )
        
        fail = run_cmds( cmds, run )
        
        if( fail ) {
            log.error( "Failed to execute commands ", cmds);
            
            stop();
        }
        
        cmds = c()
        
        # create organism.version dir
        reference_folder = paste( reference_dir, "/", organism_ver, sep="" )
        
        if( !file.exists( reference_folder )) {
            cmds = c( cmds, paste( "mkdir ", reference_folder, sep="" ) )
        }
        
        cmds = c( cmds, paste( "gunzip ", organism_ver, "\\.*", ty, "*fa.gz", sep="" ) )
        cmds = c( cmds, paste( "mv ", organism_ver, "\\.*", ty, "*fa ", reference_folder, "/", 
                        organism_ver, ".", ty, ".chromosome.fa ", sep="" ) );
        
        run_cmds( cmds, run )
        
        cmds = c()
        log.info( "Reference saved to ", reference_dir, "/", organism_ver )
        
        # Homo_sapiens have 2 PAR regions in Y and X chromosomes 
        # that are very similar. Ensembl didn't find anything more 
        # clever to do other than to remove the PARs from Y chromosome.
        # This removal resulted in that "Y" was split into 2 parts, 
        # which is fine to aligners, but cufflinks doesn't "understand" it.
        #
        #
        # The first part is the leading telomere, which has all Ns and 
        # is not used in alignment. We will cut it off from the reference. 
        #
        # The second part is the actual Y chromosome, what's left of it.
        #
        if (organism == "Homo_sapiens") {
            
            # remove the empty Y region from the reference
            #
            fixReferenceDoubleY( paste(reference_folder, "/", filename, sep="") );
        }
        
    }
    
    return( version )
}


getEnsemblAnnotation <- function( organism, version, location, run = FALSE, refresh = FALSE ) {
    
    trace.enter("getEnsemblAnnotation");
    on.exit({ trace.exit() })
    
    reference_dir = paste(location, "/annotation", sep="");
    
    
    if( version != "current" ) {
        release_ver = getMatchedText("\\d*$", version);
        ref_subver = 
        # pub/release-62/gtf/
        ref_subverpath = paste( "release-", release_ver, "/gtf/", sep="")
        
    } else {
        #"ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/"
        ref_subverpath = paste( "current_gtf/", sep="")
    }
    
    ensembl_path = paste("ftp://ftp.ensembl.org/pub/", ref_subverpath, tolower(organism), sep="");
    
    cmds = c();
    
    # create dir for all references
    if( !file.exists(reference_dir) ) {
        cmds = c( cmds, paste( "mkdir ", reference_dir, sep="" ) )
    }
    
    dirlist = paste(organism, ".dirlist.txt", sep="");
    
    # if version is current, the actual filename will be unknown
    # we're finding out the actual filename here
    #
    cmds = c( cmds, paste( "curl -s -l ", ensembl_path, "/ | grep gtf.gz > ", dirlist, sep="" ) )
    
    run_cmds( cmds, run )
    
    cmds = c()
    
    # read the file
    #
    ftp_filename = read.table( dirlist, stringsAsFactors = FALSE )[1,1]
    
    # update the version
    #
    version = substr(ftp_filename, nchar(organism)+2, nchar(ftp_filename)-nchar(".gtf.gz"))
    
    # organism.version
    #
    organism_ver = paste(organism, ".", version, sep="")
    
    # strip off the ".gz" thingy
    #
    gtf_filename = substr(ftp_filename,1, nchar(ftp_filename)-3);
    
    gtf_fullfilename = paste(reference_dir, "/", organism_ver, "/", gtf_filename, sep="");
    
    log.info( "Looking for ", gtf_filename, " ..." )
    
    if( file.exists(gtf_fullfilename) && file.exists(paste( gtf_fullfilename, ".out", sep="")) && !refresh ) {
        # do nothing
        #
        
        log.info( "File ", gtf_filename, " exists." )
    } else {
        # if the file is missing or refresh is 
        # requested, download the file and unpack it
        
        # download reference
        cmds = c( cmds, paste( "curl -O ", ensembl_path, "/", ftp_filename, sep="" ) )
        
        fail = run_cmds( cmds, run )
        
        if( fail ) {
            log.error( "Failed to execute commands ", cmds);
            
            stop();
        }
        
        cmds = c()
        
        # create organism.version dir
        reference_folder = paste( reference_dir, "/", organism_ver, sep="" )
        
        if( !file.exists(reference_folder) ) {
            cmds = c( cmds, paste( "mkdir ", reference_folder, sep="" ) )
        }
        
        if (file.exists(gtf_filename)) {
            cmds = c( cmds, paste( "rm ", gtf_filename, sep="" ) )
        }
        
        cmds = c( cmds, paste( "gunzip ", ftp_filename, sep="" ) );
        cmds = c( cmds, paste( "mv ", gtf_filename, " ", reference_folder, "/", sep="") );
        
        run_cmds( cmds, run )
        
        cmds = c();
        
        cmds = c( cmds, paste( "perl ", .path.package(package = "ArrayExpressHTS"), "/script/gtf2txt.pl ", 
                        reference_folder, "/", gtf_filename, sep="") );
        
        run_cmds( cmds, run )
        
        log.info( "Annotation saved to ", reference_dir, "/", organism_ver )
    }
    
    return( version )
}

checkReference <- function( fname ) {
    
    fai = read.table( fname, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, 
            row.names = NULL, blank.lines.skip = TRUE, sep="\t" );
    
    if (length( unique(fai[[1]]) ) != length( fai[[1]] )) {
        log.info(fname, " Inconsistency Detected");
    } else {
        log.info(fname, " OK");
    }
}    


fixReference <- function( fname ) {
    # ".fa.fai"
    
    fai = read.table( fname, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, 
            row.names = NULL, blank.lines.skip = TRUE, sep="\t" );
    
    if (length( unique(fai[[1]]) ) != length( fai[[1]] )) {
        
        fnamebackup = paste(fname, ".backup", sep="")
        
        while (file.exists(fnamebackup)) {
            fnamebackup = paste(fnamebackup, "~", sep = "");
        }
        
        file.copy(fname, fnamebackup, overwrite = TRUE);
        
        fixedfai = fai[ match( unique(fai[[1]]), fai[[1]] ), ]
        
        write.table(fixedfai, fname, quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t" );
        
        
        log.info(fname, "Fixed");
    } else {
        log.info(fname, "OK");
    }
    
}

check_indexes <- function( run = FALSE ) {
    trace.enter("check_indexes");
    on.exit({ trace.exit() })
    
    if( !file.exists(.project$reference$file) )
    indexReference( .project$organism, .project$reference$type, .project$reference$version, 
            .project$refdir, .project$aligner$type, run = TRUE ) 
    else
    log.info( "\tindexes OK" )
}


indexReference <- function( organism, type, version, location, aligner, refresh = FALSE, run = TRUE ) {
    trace.enter("indexReference");
    on.exit({ trace.exit() })
    
    owd = getwd()
    cmds = c()
    if( type == "genome" || type == "dna" ) {
        ty = "dna"
    } else { 
        ty = "cdna"
    }
    ref_dir = paste(organism, ".", version, sep="")
    gen_dir = "reference_genomes"
    filename = paste( ref_dir, ".", ty,".chromosome.fa", sep="" )
    
    if( aligner %in% c("bowtie", "tophat") ) {
        index_file = paste( location, "/tophat_indexes/", ref_dir, "/", filename, sep="" )
        log.info( "Looking for index files ", index_file, " ..." )
        
        if( file.exists(index_file) && !refresh ) {
            log.info( "File ", index_file, " exists." )
        } else {
            file = paste(location, "/tophat_indexes", sep="")
            if( !file.exists(file) )
            cmds = c( cmds, paste("mkdir ", file, sep="") )
            file = paste(location, "/tophat_indexes/", ref_dir, sep="")
            if( !file.exists(file) )
            cmds = c( cmds, paste("mkdir ", file, sep="" ) )
            run_cmds( cmds, run )
            cmds = c()
            
            if( run ) {
                file = paste( location, "/tophat_indexes/", ref_dir, sep="" )
                setwd( file )
                log.info( "Set dir ", file )
            }
            
            cmds = c( cmds, paste("bowtie-build ", location, "/", gen_dir, "/", ref_dir, "/", 
                            ref_dir, ".", ty, ".chromosome.fa ", ref_dir, ".", ty, ".chromosome", sep="") )
            cmds = c( cmds, paste( "ln ", location, "/", gen_dir, "/", ref_dir, "/", filename, " .", sep="" ) )
            run_cmds( cmds, run )
            cmds = c()
            
            if( run ) {
                file = location
                setwd( file )
                log.info( "Set dir ", file )
            }
            cmds = c( cmds, paste( "ln -s tophat_indexes bowtie_indexes" ) )
            run_cmds( cmds, run )
            cmds = c()
            
            log.info( "Indexes saved to ", location )
        }
    } else if( aligner == "bwa" ) {
        index_file = paste( location, "/bwa_indexes/", ref_dir, "/", filename, sep="" )
        log.info( "Looking for index files ", index_file, "..." )
        
        if( file.exists(index_file) && !refresh ) {
            log.info( "File ", index_file, " exists." )
        } else {
            if( run ) {
                file = location
                setwd( file )
                log.info( "Set dir ", file )
            }
            file = "bwa_indexes"
            if( !file.exists(file) )
            cmds = c( cmds, paste("mkdir", file) )
            file = paste( location, "/bwa_indexes/", ref_dir, sep="" )
            if( !file.exists(file) )
            cmds = c( cmds, paste( "mkdir ", file, sep="" ) )
            run_cmds( cmds, run )
            cmds = c()
            if( run ) {
                file = paste(location, "/bwa_indexes/", ref_dir, sep="")
                setwd( file )
                log.info( "Set dir ", file )
            }
            cmds = c( cmds, paste( "bwa index -a bwtsw ", location, "/", gen_dir, "/", 
                            ref_dir, "/", ref_dir, ".", ty, ".chromosome.fa", sep="" ) )
            cmds = c( cmds, paste( "mv ", location, "/", ref_dir, "/", filename, ".* .", sep="" ) )
            cmds = c( cmds, paste( "ln ", location,  "/", gen_dir, "/", ref_dir, "/", filename, " .", sep="" ) )
            run_cmds( cmds, run )
            cmds = c()
            
            log.info( "Indexes saved to ", location )
        }
    } else {
        log.info( "Aligner not supported" )
    }
    
    setwd( owd )
}


# organism = "Homo_sapiens", version="GRCh37.60", type="genome"/"transcriptome", location="/ebi/microarray/home/biocep/sequencing"; aligner="bowtie"
prepareReference <- function( organism, version = "current", type = c("genome", "transcriptome"), 
        location = getDefaultReferenceDir(), aligner = c("bwa", "bowtie", "tophat"), refresh = FALSE, run = TRUE ) {
    
    trace.enter("prepareReference");
    on.exit({ trace.exit() })
    
    if(missing(type)) {
        log.error("'type' is not defined");
        
        stop();
    };
    
    if(missing(aligner)) {
        log.error("'aligner' is not defined");
        stop();
    };
    
    if (!file.exists(location)) {
        log.error(location, " Not Found. Please create it.");
        stop();
    }
    
    location = normalizePath(location);
    
    version = getEnsemblReference( organism = organism, type = type, 
            version = version, location = location, refresh = refresh, run = run)
    
    indexReference( organism = organism, type = type, version = version, 
            location = location, aligner = aligner, refresh = refresh, run = run)
    
    return( version )
}

# organism = "Homo_sapiens"  version="GRCh37.60" location="/ebi/microarray/home/biocep/sequencing" aligner="bowtie"
prepareAnnotation <- function( organism, version = "current", 
        location = getDefaultReferenceDir(), refresh = FALSE, run = TRUE ) {
    
    trace.enter("prepareAnnotation");
    on.exit({ trace.exit() })
    
    if (!file.exists(location)) {
        log.error(location, " Not Found. Please create it.");
        stop();
    }
    
    location = normalizePath(location);
    
    version = getEnsemblAnnotation( organism = organism, version = version, 
            location = location, refresh = refresh, run = run )
    
    return( version )
}

