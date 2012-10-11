# Change the parameters
#
# Author: filimon, Andrew Tikhonov
#

getReadLength0 = function( fastq ) {
    result = 0;
    .C("getReadLength", as.character(fastq), as.integer(result), PACKAGE = "ArrayExpressHTS")[[2]];
}

obtainENAFastQDataFiles <- function( accession, resultset, datadir, localmode = TRUE, refresh = FALSE ) {
    trace.enter("obtainENAFastQDataFiles");
    on.exit({ trace.exit() })
    
    log.info(accession,' mapped to ',length(resultset$enaexpid), ' ENA experiments');
    
    log.info(paste(names(resultset$enaexpid), collapse=' '));
    
    fastQs = sapply(resultset$enaexpid, function(x){ x$enafastq; });
    
    for (i in 1:length(fastQs)) { 
        if (is.null(fastQs[[i]])) { 
            
            # register step status
            #
            registerExpStepWarning(METADATA_SANITY, 
                    'No runs available for ', names(fastQs[i]));
            
        }
    }
    
    # clean data folder
    incompleteFastQs = dir(datadir, pattern="fastq.gz");
    
    for(fname00 in incompleteFastQs) {
        log.info("cleaning ", fname00);
        
        suppressWarnings( file.remove( paste(datadir,"/",fname00, sep="") ) );
    }
    
    
    for(fastq in unlist(fastQs)) {
        
        splitarr = unlist(strsplit(fastq, '/'));
        
        fname = splitarr[length(splitarr)];
        
        fnameunzipped = getMatchedText("\\w*.\\w*", fname);
        
        #fnames = dir(datadir, pattern=fnameunzipped, full.names = FALSE);
        
        if (fnameunzipped %in% dir(datadir)) {
            log.info("Found ", fnameunzipped );
            
        } else {
            log.info("Getting ", fnameunzipped);
            
            if (localmode) {
                #'/nfs/era-pub/vol1/fastq/SRR017/SRR017181'
                
                fastqpath = paste(getENAFastqPathBase(), fastq, sep='/');
                
                destfile = paste(datadir, fname, sep='/');
                
                log.info( "copying .. " );
                
                cmd = paste('cp ', fastqpath, ' ', destfile, sep='');
                
                result = system(cmd);
                
                if (result != 0) {
                    # register step status
                    #
                    registerExpStepFailure(OBTAINING_RAWDATA, 
                            'Failed copying ', fastqpath, ' to ', destfile);
                    
                    stop();
                }
                
            } else {
                url = paste(getENAFastqURLBase(), fastq, sep='/');
                
                destfile = paste(datadir, fname, sep='/');
                
                log.info('downloading ', url, ' to ', destfile);
                
                result0 = try(download.file(url = url, destfile = destfile, quiet = TRUE), silent = FALSE);
                
                if (inherits(result0, 'try-error')) {
                    # register step status
                    #
                    registerExpStepFailure(OBTAINING_RAWDATA, 
                            'Failed downloading ', url, ' to ', destfile);
                    
                    stop();
                }
                
            }
            
            log.info( "unzipping .. " );
            
            cmd = paste('gunzip ', destfile, sep='');
            
            result = system(cmd);
            
            if (result != 0) {
                # register step status
                #
                registerExpStepFailure(OBTAINING_RAWDATA, 
                        'Failed unzipping ', destfile);
                
                stop();
            }
        }
    }
    
    return(TRUE);
}


getExperimentDescriptors <- function( accession, datadir, localmode = TRUE ) {
    trace.enter("getExperimentDescriptors");
    on.exit({ trace.exit() })
    
    # download IDF
    #
    #
    idflocaion = paste(datadir,paste(accession,'.idf.txt',sep=''), sep='/');
    
    
    if (localmode) {
        idfsource = getAEIDFFilename(accession);
    } else {
        idfsource = getAEIDFURL(accession);
    }
    
    if (!is.null(idfsource)) {
        
        # remove the file
        unlink(idflocaion);
        
        if (localmode) {
            log.info('Copying IDF to ', idflocaion);
            file.copy(idfsource, idflocaion, overwrite = TRUE);
            
            # fix permissions
            Sys.chmod(idflocaion, "664");
        } else {
            log.info('Downloading IDF ', idfsource);
            download.file(idfsource, idflocaion)
        }
    } else {
        # register step status
        #
        registerExpStepWarning(OBTAINING_METADATA, 
                'No IDF descriptor found.');
        
    }
    
    
    #if (file.exists(idflocaion)) {
    #    log.info('Found IDF ', idflocaion);
    #} else {
    #}
    
    # download SDRF
    # need to check 2 possible file types .sdrf and .seq.sdrf.txt
    
    sdrflocation = paste(datadir,paste(accession,'.sdrf.txt',sep=''), sep='/');

    
    if (localmode) {
        sdrfsource = getAESDRFFilename(accession);
    } else {
        sdrfsource = getAESDRFURL(accession);
    }
    
    if (!is.null(sdrfsource)) {
        
        # remove the file
        unlink(sdrflocation);
        
        if (localmode) {
            log.info('Copying SDRF to ', sdrflocation);
            file.copy(sdrfsource, sdrflocation, overwrite = TRUE);
            
            # fix permissions
            Sys.chmod(idflocaion, "664");
        } else {
            log.info('Downloading SDRF ', sdrfsource);
            download.file(sdrfsource, sdrflocation)
        }
    } else {
        # register step status
        #
        registerExpStepFailure(OBTAINING_METADATA, 
                'No SDRF descriptor found');
        
        stop();
    }
    
    
    #if (file.exists(sdrflocation)) {
    #    log.info('Found SDRF ', sdrflocation);
    #    
    #} else {
    #}

    return (list('sdrffname' = sdrflocation, 'idffname' = idflocaion))
}

initDefaultProject <- function() {
    trace.enter("initDefaultProject");
    on.exit({ trace.exit() })
    
    p0 = list();
    p0$reference = list();
    p0$annot = list();
    p0$count = list();
    p0$pairing = list();
    p0$aligner = list();
    p0;
}

createDefaultProject <- function(projname, organism, project0) {

    trace.enter("createDefaultProject");
    on.exit({ trace.exit() })
    
    log.info("Creating project: ", projname);
    
    project = project0;
    
    project$name = projname;
    
    project$projectdir = paste(project$basedir, projname, sep="/")
    
    project$organism = organism;
    
    project$reference$version = getCurrentRefVersion( organism, project0$refdir );
    
    #project$projectdir = projname;
    
    dir.create(project$projectdir, showWarnings = FALSE)
    
    aligner_outdir = paste(project$aligner$type, "_out", sep="");
    
    project$aligner$out_dir = paste(project$projectdir, aligner_outdir, sep="/");
    
    dir.create( project$aligner$out_dir, showWarnings = FALSE );
    
    project = fixedOptions( project );
    
    project$report = c(
        raw_report_dir      = paste( project$projectdir, "/report", sep="" ),
        aligned_report_dir  = paste(project$aligner$out_dir, "/report", sep=""),
        compared_report_dir = paste( project$basedir, "/compare_report", sep="" ) );
    
    for( rdir in project$report ) {
        dir.create( rdir, showWarnings = FALSE );
    }
    
    return(project);
    
}

getFirstAvailableLayout <- function(layoutinfo) {
    # get first available layout info
    #
    for(item in layoutinfo) { 
        if(!is.null(item)) { 
            return(item); 
        }
    }
    
    return(NULL);
}


# will create the file hierarchy in the current directory by default
createAEprojects <- function( accession, options = getDefaultProcessingOptions(),
        dir = getwd(), refdir = getDefaultReferenceDir(), localmode = TRUE ) {
    
    trace.enter("createAEprojects");
    on.exit({ trace.exit() });
    
    projects = list();
    
    # check if reference folder is there
    if(!file.exists(refdir)) {
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "Reference folder ",refdir," does not exist.");
        
        stop();
    }
    
    # normalize paths
    dir = normalizePath(dir);
    refdir = normalizePath(refdir);
    
    # register step status
    #
    registerExpStepStarted(OBTAINING_METADATA);
    registerExpStepStarted(METADATA_SANITY);
    
    
    # after checking reference folder
    # load & cache all supported organisms
    #
    
    scanSupportedOrganisms(refdir);
    
    resultset = NULL;
    
    # define folders
    basedir = setupBaseFolder(accession, dir, getSubmID());
    psrdir  = setupPSRFolder(accession, dir, getSubmID());
    datadir = setupDataFolder(accession, dir);
    
    descriptors = getExperimentDescriptors( accession, datadir = datadir, localmode = localmode );
    resultset   = mapAEtoENAviaHTTP( accession, descriptors$sdrffname );
    
    if (is.null(resultset$enaexpid)) {
        # register step status
        #
        registerExpStepFailure(METADATA_SANITY, 
                'No experiment data found.');

        stop();
    }
    
    # register step status
    #
    registerExpStepStarted(OBTAINING_RAWDATA);
    
    result = obtainENAFastQDataFiles(accession, resultset, datadir=datadir, localmode=localmode);
    
    if (!result) {
        # register step status
        #
        registerExpStepFailure(OBTAINING_RAWDATA, 
                'Failed to obtain experiment raw data.');
        
        stop();
    }
    
    # register step status
    #
    registerExpStepCompleted(OBTAINING_RAWDATA);
    
    
    expidcnt = length(resultset$enaexpid);
    expidnames = names(resultset$enaexpid);
    
    # get layout info
    #
    layoutinfo = sapply(resultset$enaexpid, function(x){ x$layout$info })
    
    availablelayoutinfo = getFirstAvailableLayout(layoutinfo);
    
    projectnames = c();
    
    for (i in 1:expidcnt) {
        # expidnames[i]
        expid = resultset$enaexpid[[i]];
        
        if (is.null(expid$enafastq)) {
            
            # register step status
            #
            registerExpStepWarning(METADATA_SANITY, expidnames[i], " No ENA data found ");
            
        } else {
            log.info(expidnames[i]);
            log.info('   RUNID:', expid$enarunid);
            log.info('   LAYOUT:', expid$layout$name);
            
            if (expid$layout$name == "PAIRED") {
                log.info('      NOMINAL_LENGTH:', expid$layout$info['NOMINAL_LENGTH']);
                log.info('      ORIENTATION:', expid$layout$info['ORIENTATION']);
                log.info('      NOMINAL_SDEV:', expid$layout$info['NOMINAL_SDEV']);
            }
            
            log.info('   SAMPLE:', expid$samplename, ' TAXON:', expid$sampletaxon);
            
            # reset organism
            #
            organism = NULL;
            
            #organismfound = FALSE;
            #organismdefined = FALSE;
            
            
            ena.organism = expid$samplename;
            
            # check organism from ENA
            #
            if (!is.null(ena.organism) && nchar(ena.organism) > 0) {
                
                log.info(expidnames[i], " Organism defined in ENA: ", ena.organism);
                
                ena.organism = sub( " ", "_", ena.organism)

                if (isOrganismSupported(ena.organism, refdir)) {
                    log.info("Using organism: ", ena.organism);
                    organism = ena.organism;
                
                } else {
                    
                    # register step status
                    #
                    registerExpStepWarning(METADATA_SANITY, 
                            expidnames[i], " Organism specified in ENA ", 
                            ena.organism," is not supported");
                    
                }
                
            } else {
                
                # register step status
                #
                registerExpStepWarning(METADATA_SANITY, 
                        expidnames[i], " No organism found in ENA");
                
            }
            
            # if ENA has no organism defined, 
            # try getting it from SDRF 
            #
            if (is.null(organism)) {
                
                # get organism from SDRF
                #
                indexes = grep(expidnames[i], resultset$sdrfexpid);
                
                if (length(indexes) > 0 && nchar(resultset$sdrforganism[ indexes[[1]] ]) > 0) {
                    
                    sdrf.organism = resultset$sdrforganism[ indexes[[1]] ];
                    
                    log.info(expidnames[i], " Organism defined in SDRF: ", sdrf.organism);
                    
                    sdrf.organism = sub( " ", "_", sdrf.organism)
                    
                    if (isOrganismSupported(sdrf.organism, refdir)) {
                        
                        log.info("Using organism: ", sdrf.organism);
                        organism = sdrf.organism;
                    
                    } else {
                        
                        # register step status
                        #
                        registerExpStepWarning(METADATA_SANITY, 
                                expidnames[i], " Organism in SDRF ", 
                                sdrf.organism," is not supported");
                        
                    }
                    
                } else {
                    # register step status
                    #
                    registerExpStepWarning(METADATA_SANITY, 
                            expidnames[i], " No organism found in SDRF");
                    
                }
            }
            
            # check if organism was not fond
            #
            if (is.null(organism)) {
                
                # register step status
                #
                log.info(expidnames[i], " Organism not defined.");
                log.info("Prepare reference and anotation data for the organism.");
                log.info("Otherwise make sure proper organism is defined in SDRF.");
                log.info("");
                log.info("Supported organisms:" );
                log.info(paste("'", names(getSupportedOrganisms( refdir )), "' ", sep=''));
                
                registerExpStepFailure(REFERENCE_SANITY, 
                        expidnames[i], " Organism not defined.");
                
                stop();
            }
            
            # get organism version
            #version = getCurrentRefVersion(organism);
            
            # init meta project
            #
            project0 = initDefaultProject();
            
            # init default project
            #
            project0$psrdir = psrdir;
            project0$basedir = basedir;
            project0$datadir = datadir;
            project0$refdir = refdir;
            project0$organism = organism;
            project0$enaexpid = expidnames[i];
            
            project0 = setUserOptions( project0, options );
            
            for (fastq in expid$enafastq) {
                
                projname = unlist(strsplit(fastq, '/'))[3];
                
                #log.info(projname);
                
                if (!(projname %in% projectnames)) {
                    
                    projects[[ projname ]] = createDefaultProject(projname, organism, project0);
                    
                    # quality scores from ENA are 
                    # always Phred+33, which is "FastqQuality"
                    #
                    #
                    projects[[ projname ]]$qual_type = "FastqQuality";
                    
                    projectnames = c(projectnames, projname);
                    
                    # check if paired parameters need to
                    # be "retro-fitted" back into the project
                    #
                    
                    if (expid$layout$name == "PAIRED") {
                        
                        # get fastq filename
                        #
                        fqfilename = projects[[projname]]$fq_files[1];
                        
                        # get read length
                        #
                        readlength = getReadLength0(fqfilename);
                        
                        log.info( fqfilename, " Determined Read Length: ", readlength);
                        
                        # check if the layout info is available
                        #
                        if (!is.null(expid$layout$info)) {
                            
                            explayoutinfo = expid$layout$info
                        
                        } else {
                            
                            # register step status
                            #
                            registerExpStepWarning(METADATA_SANITY, 
                                    projname, " NO LAYOUT INFO ");
                            
                            
                            if (!is.null(availablelayoutinfo)) {
                                
                                explayoutinfo = availablelayoutinfo;
                                
                            } else {
                                # TODO: 
                                # read this from defaults 
                                #
                                explayoutinfo = c("500", "25");
                                names(explayoutinfo) = c("NOMINAL_LENGTH", "NOMINAL_SDEV");
                            }
                        }
                            
                        expNominalLength = as.integer( explayoutinfo['NOMINAL_LENGTH'] )
                        expNominalSdev = as.integer( explayoutinfo['NOMINAL_SDEV'] );
                        
                        # if user has not defined 
                        # the insert size using options
                        #
                        if(is.null(options$insize)) {
                            # update the insert size
                            #
                            
                            if (!is.null(expNominalLength) && !is.na(expNominalLength)) {
                                projects[[projname]]$pairing$insize = expNominalLength - readlength * 2;
                            } else {
                                projects[[projname]]$pairing$insize = NULL;
                            }
                        }
                        
                        # if size deviation 
                        # 
                        if( is.null(options$insizedev) ) {
                            if (!is.null(expNominalSdev) && !is.na(expNominalSdev)) {
                                projects[[projname]]$pairing$insizedev = expNominalSdev;
                            } else {
                                projects[[projname]]$pairing$insizedev = NULL;
                            }
                        }
                        
                    }
                            
                    # pairing type should have been 
                    # setup to this moment, check it
                    #
                    if ((projects[[projname]]$pairing$type == "PE" && expid$layout$name != "PAIRED") || 
                    (projects[[projname]]$pairing$type == "SR" && expid$layout$name != "SINGLE")) {
                        
                        # register step status
                        #
                        registerExpStepFailure(METADATA_SANITY, 
                                "LAYOUT & DATA Mismatch, Data Layout: ", 
                                project0$pairing$type, " Experiment Layout: ", expid$layout$name);
                        
                        stop();
                    }
                
                } else {
                    log.info(projname, " already defined");
                }
            }
        }
    }
    
    metadata = list();
    
    metadata$psrdir = psrdir;
    metadata$basedir = basedir;
    metadata$datadir = datadir;
    metadata$refdir = refdir;
    metadata$resultset = resultset;
    
    attr(projects, 'metadata') = metadata;
    
    # register step status
    #
    registerExpStepCompleted(OBTAINING_METADATA);
    registerExpStepCompleted(METADATA_SANITY);
    
    
    return(projects);
}


createFastQProjects <- function( accession, organism, quality, 
        options = getDefaultProcessingOptions(), dir = getwd(), refdir = getDefaultReferenceDir() ) {
    trace.enter("createFastQProjects");
    on.exit({ trace.exit() })
    
    # define projects
    projects = list()
    
    # register step status
    #
    registerExpStepStarted(OBTAINING_METADATA);
    registerExpStepStarted(METADATA_SANITY);
    
    # normalize paths
    dir = normalizePath(dir);

    
    # define folders
    basedir = setupBaseFolder(accession, dir, getSubmID());
    psrdir  = setupPSRFolder(accession, dir, getSubmID());
    datadir = setupDataFolder(accession, dir);
    
    # check if reference folder is there
    if( !file.exists(refdir) ) {
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "Reference folder ",refdir," does not exist.");
        
        stop();
    }
    
    # normalize reference folder
    refdir = normalizePath(refdir);
    
    # create projects
    #
    projectnames = c();
    organismnames = c();
    baselengths = c(NA);
    
    # load supported organisms
    scanSupportedOrganisms( refdir );    
    
    # read SDRF
    # Parse sdrf and create a list with a project per lane
    sdrffilename = dir( datadir, pattern="sdrf", full.names = TRUE )
    
    if( length(sdrffilename) > 0 ) {
        
        log.info( "Creating projects from SDRF file" );

        
        sdrf = readSDRF( sdrffilename[1] )
        
        index1 = grep("Array.Data.File", names(sdrf@data)) # direct submission, array-based format
        index2 = grep("ENA_RUN", names(sdrf@data))         # AE-ENA experiment format
        index3 = grep("Sample", names(sdrf@data))          # manual sdrf format
        
        # AE-ENA experiment sdrf format
        #
        if (length(index2) > 0 && length(projectnames) == 0) {
            # get names from ENA_RUN
            
            projectnames = sdrf@data[[index2]]
            projectnames = unique( projectnames )
            
            for( name in projectnames ) {
                organismnames = c(organismnames, sub( " ", "_", 
                sdrf@data$Characteristics.Organism.[grep( name, sdrf@data[[index2]] )[1]] ));
            }
        }
        
        # direct submission, array-based sdrf format
        #
        if (length(index1) > 0 && length(projectnames) == 0) {
            # try getting names from Array.Data.File
            #
            ind = regexpr( "_[12].f", sdrf@data$Array.Data.File ) # files have to has extension starting with 'f'
        
            if (length(ind) > 0) {
                ind = sapply( 1:length(ind), function(x) { 
                    if(ind[x] < 0) 
                    regexpr( ".f", sdrf@data$Array.Data.File[x], fixed = TRUE) 
                else ind[x] } )
                
                for( i in 1:length(ind) ) {
                    projectnames = c( projectnames, substr( sdrf@data$Array.Data.File[i], 1, ind[i]-1 ) )
                }
                
                projectnames = unique( projectnames )
                
                for( name in projectnames ) {
                    organismnames = c(organismnames, sub( " ", "_", 
                    sdrf@data$Characteristics.Organism.[grep( name, sdrf@data$Array.Data.File )[1]] ));
                }
            }
        }
        
        # manual sdrf format
        #
        # "Sample" "Organism"     "Base.Length"
        #  1286     Homo sapiens   150
        #  1287     Homo sapiens   150
        if (length(index3) > 0 && length(projectnames) == 0) {
            
            projectnames = sdrf@data[[index3]]
            #projectnames = unique( projectnames )
            
            i0 = grep("Organism", names(sdrf@data))
            
            if (length(i0 > 0)) {
                # read organisms
                #
                organismnames = sub( " ", "_", sdrf@data[[i0]] )
            } else {
                # duplicate master organism
                #
                organismnames = sub( " ", "_", rep(organism, length(projectnames)))
            }

            i1 = grep("Base.Length", names(sdrf@data))
            
            if (length(i1 > 0)) {
                # read real lengths
                #
                baselengths = sdrf@data[[i1]];
            }

        }
        
        if (length(projectnames) == 0) {
            
            # register step status
            #
            registerExpStepWarning(METADATA_SANITY, 
                    projname, "Cannot determine projects from SFRF");
            
        }
        
    }
     
    if (length(projectnames) == 0 || length(sdrffilename) == 0) {
    
        if (organism == "automatic") {
            # register step status
            #
            log.info( "'Automatic' organism can be used along with SDRF descriptors.");
            log.info( "Please specify 'organism' precisely. Supported organisms:" );
            log.info( "'", names(getSupportedOrganisms( refdir )), "'" );
            
            registerExpStepFailure(PARAMETER_SANITY, 
                    "No SDRF, cannot use 'automatic' organism. Please specify it.");
            
            stop();
        }
        
        log.info( "Creating projects from data files" );
        
        
        # list all .fastq and .fq files
        fastqfnames = c(dir( datadir, pattern="*.fq", full.names = FALSE ), 
        dir( datadir, pattern="*.fastq", full.names = FALSE ));
        
        if( length(fastqfnames) == 0 ) {
            
            # register step status
            #
            registerExpStepFailure(OBTAINING_RAWDATA, 
                    "No .fastq or .fq file found in ", datadir);
            
            stop();
        }
        
        if (!isOrganismSupported(organism, refdir)) {
            
            # register step status
            #
            log.info( "Organism '", organism,"' is not supported." );
            log.info( "Supported organisms are: " );
            log.info( "'", names(getSupportedOrganisms( refdir )), "'" );
            
            registerExpStepFailure(OBTAINING_METADATA, 
                    "Organism '", organism,"' is not supported.");
            
            stop();
        }
        
        # get organism version
        # version = getCurrentRefVersion(organism);
        
        # read data folder, get all fastq files
        #fnamematches = regexpr( "\\w*", fastqfnames );
        matches = regexpr( "_[12].f", fastqfnames );
        matches = sapply( 1:length(matches), function(x) { 
        if(matches[x] < 0) regexpr( ".f", fastqfnames[x], fixed = TRUE) else matches[x] } );
        
        # get the project names
        # projectnames = substr(fastqfnames, 1, attr(fnamematches, "match.length"));
        projectnames = unique(substr(fastqfnames, 1, matches-1));
        organismnames = rep(organism, length(projectnames));
    }

    log.info( "Found ", length(projectnames), " projects" );
    
    if (organism != "automatic") {
        indices = c();
        
        for(org in organism) {
            orgindices = grep(org, organismnames);
            
            if (length(orgindices) == 0) {
                
                # register step status
                #
                registerExpStepWarning(METADATA_SANITY, 
                        "Warning, organism ", org, " not found in SDRF");
                
            } else {
                indices = c(indices, grep(org, organismnames));
            }
        }
        
        # backup projects
        projectnamesBackup = projectnames;
        organismnamesBackup = organismnames;
        
        # correct projects
        projectnames = projectnames[indices];
        organismnames = organismnames[indices];
        
        if (length(projectnamesBackup) != length(projectnames)){
            log.info( "Selected ", length(projectnames), " projects using organism(s) ", 
                    paste(organism, collapse=" ") );
            
            # register step status
            #
            registerExpStepWarning(PARAMETER_SANITY, 
                    length(projectnamesBackup) - length(projectnames), 
                    " projects filtered out by organism filter and will not be computed");
            
        }
    }
    
    
    project0 = initDefaultProject();
    
    # init default project
    #
    project0$basedir = basedir;
    project0$psrdir = psrdir;
    project0$datadir = datadir;
    project0$refdir = refdir;
    
    project0 = setUserOptions( project0, options );
    
    for( i in 1:length(projectnames) ) {
        
        projname = projectnames[i];
        
        projects[[ projname ]] = createDefaultProject(projname, organismnames[i], project0);
        
        if (quality == "auto") {
            # detect quality
            #
            projects[[ projname ]]$qual_type = get_fastq_quality0( projects[[ projname ]]$fq_files[1] );
        } else {
            # set user selected quality
            projects[[ projname ]]$qual_type = quality;
        }
        
        # check if paired end info 
        # needs to be defined
        #
        if (!is.null(baselengths[i]) && !is.na(baselengths[i]) && projects[[ projname ]]$pairing$type == "PE") {
            # get fastq filename
            #
            fqfilename = projects[[ projname ]]$fq_files[1];
            
            # get read length
            #
            readlength = getReadLength0( fqfilename );
            
            log.info( fqfilename, " Determined Read Length: ", readlength);
            
            # if user has not defined 
            # the insert size using options
            #
            if(is.null(options$insize)) {
                # update the insert size
                #
                projects[[ projname ]]$pairing$insize = baselengths[i] - readlength * 2;
            }
            
            # if size deviation 
            # 
            if( is.null(options$insizedev) ) {
                projects[[projname]]$pairing$insizedev = as.integer( 25 );
            }
        }
    }
    
    if (quality == "auto") {
        
        scales0 = unique(unlist(sapply(projects, function(x){ x$qual_type; })));
    
        if ( length(scales0) > 1 && options$aligner == "tophat" 
                && !getPipelineOption("ignorequalityerrors") ) {
            
            
            # register step status
            #
            log.info( "\tQuality scales of several fastq files are different. " );
            log.info( "\tTry increasing the 'fastqreadmax' pipeline option, using " );
            log.info( "\tgetPipelineOptions('fastqreadmax' = 100000), or disable" ); 
            log.info( "\tquality errors using option 'ignorequalityerrors'." );
            
            registerExpStepFailure(OBTAINING_METADATA, 
                    "Quality canot be automatically determined.");
            
            stop();
        }
    }
    
    metadata = list();
    
    metadata$psrdir = psrdir;
    metadata$basedir = basedir;
    metadata$datadir = datadir;
    metadata$refdir = refdir;
    
    attr(projects, 'metadata') = metadata;
    
    # register step status
    #
    registerExpStepCompleted(OBTAINING_METADATA);
    registerExpStepCompleted(METADATA_SANITY);
    
    return( projects );
}

setUserOptions <- function( project, options ) {
    trace.enter("setUserOptions");
    on.exit({ trace.exit() })
    
    log.info( "\tProcessing Options:" );
    
    # Strand specific
    #
    #
    
    if( is.null(options$stranded) | !is.logical( options$stranded ) ) {
        project$stranded = FALSE;
    } else {
        project$stranded = options$stranded;
    }
    
    log.info( "\tStrand specific:\t",       project$stranded );
    
    # Reference Type
    #
    #
    
    if( is.null(options$reference) ) {
        project$reference$type = "genome"
    } else {
        project$reference$type = options$reference;
    }
    
    log.info( "\tReference type:\t",   project$reference$type);
    
    
    # Aligner
    #
    #
    
    if( is.null(options$aligner) ) {
        project$aligner$type = "tophat";
    } else {
        project$aligner$type = options$aligner;
    }
    
    log.info( "\tAligner:\t",               project$aligner$type);
    
    # Aligner options
    #
    #
    
    if( is.null(options$aligner_options) ) {
        project$aligner$options = getAlignerDefaultOptions( project$aligner$type, project$reference$type );
        project$aligner$override_options = FALSE;
    } else {
        project$aligner$options = options$aligner_options;
        project$aligner$override_options = TRUE;
    }
    
    log.info( "\tAligner options:\t");
    for(optname in names(project$aligner$options)) {
        log.info( "\t\t", optname, " = ", project$aligner$options[[optname]]);
    }
    
    # Insert Size
    #
    #
    
    # warn only for tophat because bwa and bowtie estimate these values automatically
    if(!is.null(options$insize)) {
        project$pairing$insize = options$insize;
        log.info( "\tInsert size:\t", project$pairing$insize);
    } else if((project$aligner$type != "bwa") & (project$aligner$type != "bowtie")) {
        log.info( "\tInsert size:\t", "Not defined. Assuming: ", "(default parameter for TopHat)");
    }
    
    
    # Insert Size Deviation
    #
    #
    
    # warn only for aligners other than bwa because bwa estimates this values automatically
    if( is.null(options$insizedev) & (project$aligner$type != "bwa") ) {
        log.info( "\tInsert stdev:\t", "Not defined. Assuming: ", "(default parameter for TopHat)");
    } else {
        project$pairing$insizedev = options$insizedev;
        log.info( "\tInsert stdev:\t", project$pairing$insizedev);
    }
    
    
    #if( is.null(options$reference_version) ) {
    #    project$reference$version = getCurrentRefVersion( project )
    #} else {
    #    # TODO: looks like a bug here, need to check
    #    project$reference$version = options$reference_version[project$organism]
    #}
    
    # Count Method
    #
    #
    
    if( is.null(options$count_method) ) {
        project$count$method = "cufflinks"
    } else {
        project$count$method = options$count_method
    }
    
    log.info( "\tCount method:\t",     project$count$method);

    
    # Count Feature & Count options
    #
    #
    
    if( is.null(options$count_feature) ) {
        project$count$feature = "transcript";
    } else {
        project$count$feature = options$count_feature;         # if transcript set also countOptions as options for cufflinks
        project$count$options = options$count_options;
    }

    log.info( "\tFeature to count:\t", project$count$feature);
    log.info( "\tCount options:\t",    project$count$options);

    # N11n (Normalization)
    #
    #
    
    if( is.null(options$normalisation) ) {
        project$count$normalisation = "rpkm"
    } else {
        project$count$normalisation = options$normalisation;
    }
    
    log.info( "\tNormalisation:\t",       project$count$normalisation);
    
    # S13n (Standartisation)
    # 
    #
    
    if( is.null(options$standardise) ) {
        project$count$standardise = FALSE;
    } else {
        project$count$standardise = options$standardise
    }
    
    log.info( "\tStandardisation:\t",     project$count$standardise);
    
    
    # Filtering
    #
    #
    
    if( is.null(options$filtering_options) ) {
        project$filtering_options = getDefaultFilteringOptions()
    } else {
        project$filtering_options = mergeOptions(getDefaultFilteringOptions(), options$filtering_options);
    }
    
    project$filter = options$filter;

    log.info( "\tFiltering :\t",          project$filter);
    log.info( "\tFiltering optons:");
    for(optname in names(project$filtering_options)) {
        log.info( "\t\t", optname, " = ", project$filtering_options[[optname]]);
    }
    
    # Validation
    #
    #
    
    if( project$count$method == "mmseq" && !(project$aligner$type %in% c("tophat", "bowtie")) ) {
        
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "MMSEQ can only be used with the output of Tophat or Bowtie");
        
        stop();
    }

    if( project$count$method == "cufflinks" && !(project$aligner$type == "tophat") ) {
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "Cufflinks can only be used with the output of Tophat");
        
        stop();
    }

    if( project$count$method == "count" && project$reference$type != "transcriptome" ) {
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "Count can only be used using a transcriptome as reference");
        
        stop();
    }

    if( project$reference$type == "transcriptome" && project$count$feature != "transcript" ) {
        # register step status
        #
        registerExpStepFailure(PARAMETER_SANITY, 
                "When using the transcriptome as the reference only transcript level counts are allowed" );
        
        stop();
    }
    
    return(project)
}

# used by createAEprojects and createProjectOptions
fixedOptions <- function( project ) {
    trace.enter("fixedOptions");
    on.exit({ trace.exit() })
    
    project$aligner$threads = 16;
    
    # stable version workflow
    # new version
    project$aligner$files = c( "accepted_hits.bam", "accepted_hits.filt.sorted.bam", "accepted_hits.sortednames.bam" )
    
    project$aligner$bamtags = c(  
    # SAM format defined tags, 
    "NM", #NM: number of nucleotide differences
    "MD", #MD: string for mismatching positions
    "XA", # BWA: alternative alignments
    "XS", # TopHat: strand, BWA: suboptimal alignments
    "NS" ) # ?
    
    project$annot$type = "gtf"            # gtf file or download from biomaRt
    
    project$count$files = c( exon="exon.expr", transcript="transcripts.expr", gene="genes.expr" )
    
    ###########################################################################
    # NO need to CHANGE beyond this point (in principle)
    ###########################################################################
    
    organismversion = paste(project$organism, ".", project$reference$version, sep='');
    
    project$aligner$indexes_dir = paste( project$refdir, '/', project$aligner$type, "_indexes/", organismversion, sep="" )
    
    if(!file.exists( project$aligner$indexes_dir ) && project$aligner$type != "custom" ) {
        # register step status
        #
        registerExpStepFailure(REFERENCE_SANITY, 
                'Could not find ', project$aligner$indexes_dir );
        
        stop();
    }
    
    project$reference$chrRepresentation = "chr"
    
    project$reference$file = paste( project$aligner$indexes_dir, "/", 
                                    organismversion, ".", reftype_to_filename( project$reference$type ), sep="" )
    
    project$annot$folder = paste( project$refdir, "/annotation/", organismversion, sep="" )
    
    fail = download_annotation( project, run = TRUE )
    
    if( fail ) {
        
        registerExpStepFailure(REFERENCE_SANITY, 
                'Could not download annotation');
        
        stop();
    }
    
    project$fq_files = c( dir( project$datadir, pattern=paste(project$name, "_[12]", ".f", sep=""), full.names = TRUE ),
                          dir( project$datadir, pattern=paste(project$name, ".f", sep=""), full.names = TRUE ) )
    
    project$fq_files = verify_fqfiles( project$fq_files )
    
    if(length(project$fq_files) > 1) {
        project$pairing$type = "PE"            # S(ingle)R(read) or P(aired)E(nd)
    } else {
        project$pairing$type = "SR"
    }
    
    # anything beyong here must go into the "ignore" list in compareProjects()
    project$aux_files = c( 
    fqAsShortReadQ     = paste( project$projectdir, "/reads.RData", sep="" ),
    shortReadQtab      = paste( project$projectdir, "/readstab.RData", sep="" ),
    fqAsDataFrame      = paste( project$projectdir, "/fqAsDataFrame.RData", sep="" ), 
    
    alignedRead        = paste( project$aligner$out_dir, "/alignedRead", project$pairing, ".RData", sep="" ),
    alignedGenRanges   = paste( project$aligner$out_dir, "/alignedGenRanges", project$pairing, ".RData", sep="" ), # unused
    hitCount           = paste( project$aligner$out_dir, "/hitCount", project$pairing$type, ".RData", sep="" ), 
    
    idsAsDataFrame     = paste( project$projectdir, "/idtab.RData", sep=""),
    annotAsIRange      = paste( project$annot$folder, "/annotIR.RData", sep="" ), 
    annot              = paste( project$annot$folder, "/annot.RData", sep="" ), 
    bamAsList          = paste( project$aligner$out_dir, "/bamAsList.RData", sep="" ), # unused 
    bamNamesAsList     = paste( project$aligner$out_dir, "/bamNamesAsList.RData", sep="" ),
    bamAsIR            = paste( project$aligner$out_dir, "/bamAsIR.RData", sep="" ),
    countsPerFeature   = paste( project$aligner$out_dir, "/counts.RData", sep="" ),
    uniqueReads        = paste( project$aligner$out_dir, "/unique.RData", sep="" ),
    indexes            = paste( project$aligner$out_dir, "/indexes.RData", sep=""), # unused ?
    
    seqOccurrenceTab   = paste( project$aligner$out_dir, "/seqOccurrenceTab.RData", sep=""),
    polyCount          = paste( project$aligner$out_dir, "/polyCount.RData", sep="" ), # unused
    pileup             = paste( project$aligner$out_dir, "/pileup.RData", sep=""),
    
    avgQuals           = paste( project$projectdir, "/avgQuals.RData", sep="" ),
    
    dusty              = paste( project$projectdir, "/dustyScore.RData", sep=""),
    efpkm              = paste( project$projectdir, "/efpkm.RData", sep=""),
    ecount             = paste( project$projectdir, "/ecount.RData", sep=""))
    
    return( project )
}

getDefaultProcessingOptions <- function() {
    trace.enter("getDefaultProcessingOptions");
    on.exit({ trace.exit() })
    
    options = list( 
            stranded              = FALSE,
            insize                = NULL,
            insizedev             = NULL,
            reference             = "genome",
            aligner               = "tophat",
            aligner_options       = NULL,
            count_feature         = "transcript",
            count_options         = "",
            count_method          = "cufflinks",
            filter                = TRUE,
            standardise           = FALSE,
            normalisation         = "rpkm")

    return( options )
}

mergeOptions <- function(baseoptions, newoptions) {
    trace.enter("mergeOptions");
    on.exit({ trace.exit() })
    
    for(n in names(newoptions)) {
        baseoptions[[n]] = newoptions[[n]];
    }
    
    return(baseoptions);
}


getDefaultFilteringOptions <- function() {
    trace.enter("getDefaultFilteringOptions");
    on.exit({ trace.exit() })
    
    opt = list( 
            mismatches = 2, 
            chr_ignore = c("MT"), 
            minqual    = 10, 
            minmapq    = 1, 
            maxN       = 2, 
            maxpol     = 0.75,         # either a proportion (<0) or an integer
            duplicates = "remove", 
            multihits  = "remove", 
            gapped     = "remove")     # either "keep" or "remove"
    return( opt )
}

getAlignerDefaultOptions <- function( aligner, reference ) {
    trace.enter("getAlignerDefaultOptions");
    on.exit({ trace.exit() })
    
    tophat_default_options = list() 
    bwa_default_options = list()
    bowtie_default_options = list()
    custom_default_options = list()
    
    tophat_default_options[["--segment-mismatches"]]=2   # mismatches per segment, so a read of length 'n' will have up to n/25*2 mismatches    
    tophat_default_options[["--segment-length"]]=25      # number of bases to chop the reads into...
    tophat_default_options[["--mate-inner-dist"]]=0      # required for PE, no default
    tophat_default_options[["--mate-std-dev"]]=20        # required for PE
    tophat_default_options[["--min-anchor"]]=8           # min 3
    tophat_default_options[["--splice-mismatches"]]=0        
    tophat_default_options[["--min-intron"]]=50        
    tophat_default_options[["--max-intron"]]=as.integer(5000000) # the larger this number the slower the algorithms run
    #tophat_default_options[[""]] = "--solexa1.3-quals"
    tophat_default_options[["--min-isoform-fraction"]]=0.15     # 0 disables the filter
    tophat_default_options[["--max-multihits"]]=40    
    
    bwa_default_options[["-n"]]=0.04      # maximum edit distance, either a percentage of the sequence or an integer, 
                                          # e.g. 0.02 means 2% of the sequence length, 2 means only 2 mismatches
    bwa_default_options[["-l"]]=NULL      # seed length, disabled by default
    bwa_default_options[["-k"]]=2         # seed mismatches
    bwa_default_options[["hits"]]=40      # in samse and sampe: maximum number of alignements to output
    
    # return only the best alignments with the least number of mismatches, 
    # --best --strata only works for single reads
    # but can be left as is for paired-end, only it won't work, further filtering is then needed
    # -a tells to output all valid alignments
    # -S tells bowtie to output in the SAM format
    bowtie_default_options[[""]]="--best --strata -a -S"
    bowtie_default_options[["-m"]]=40      # maximum number of alignements per read
    
    ##### mutually exclusive section
    bowtie_default_options[["-n"]]=2       # maximum number of mismatches in the seed
    bowtie_default_options[["-l"]]=50      # seed length
    bowtie_default_options[["-e"]]=70      # the sum of qualities at mismatches bases must be under this
    # OR
    #bowtie_default_options[["-v"]]=2      # maximum number of mismatches, ignore qualities
    #bowtie_default_options[["-I"]]=50     # minimum intron size
    #bowtie_default_options[["-X"]]=as.integer(500000)    # maximum intron size
    ##### 
    
    bwa_default_options[["-o"]]=NULL       # maximum number of gap opens
    bwa_default_options[["-e"]]=NULL       # maximum number of gap extensions
    if( reference == "transcriptome" ) {
        bwa_default_options[["-o"]]=0      # maximum number of gap opens
        bwa_default_options[["-e"]]=0      # maximum number of gap extensions
    } 
    
    return( get( paste(aligner, "_default_options", sep="") )  )
}

projectSummary <- function( project ) {
    trace.enter("projectSummary");
    on.exit({ trace.exit() })
    
    ignore = c( "aux_files", "report_dir" )
    
    for( n in names(project)[!(names(project) %in% ignore)] ) {
        if( length(grep("dir",n)) == 0 ) {
            log.info( n,":" )
            if( !is.null(names(project[[n]])) ) {
                log.info( "" )
                for( o in names(project[[n]]) )
                log.info( "\t", o, ":\t", project[[n]][[o]], "" )
            } else
            log.info( "\t", project[[n]], "" )
        }
    }
}

compareProjects <- function( p1, p2, verbose = FALSE ) {
    #trace.enter("compareProjects");
    #on.exit({ trace.exit() })
    
    ignore = c( "out_dir", "aux_files", "report", "sample", "count" )
    if( !is.list(p1) || !is.list(p2) ) {
        if( !is.list(p1) && !is.list(p2) ) {
            #if( verbose )
            #log.info( "\tComparing ", names(p1), ": ", p1, " with ", names(p2), ": ", p2, "" )
            # it will be accepted even if the old project is lacking some options 
            if( all(p1 == p2 || is.null(p2)) ) {
                return(TRUE)
            } else {
                if( verbose )
                    log.info( "\t", p1, " != ", p2 )
                return(FALSE)
            }
        } else {
            if( verbose )
                log.info( "\t", p1, " != ", p2 )
            return(FALSE)
        }
    } else {
        return( all( sapply( names(p1), function(sp) if( sp %in% ignore ) TRUE else compareProjects( p1[[sp]], p2[[sp]], verbose=verbose ) ) ) )
    }
}

scanSupportedOrganisms <- function( refdir ) {
    # beware: this function loads supported organisms and 
    # stores them in the package internal variable storage, 
    # whereas getSupportedOrganisms uses the stored variable
    # 
    # This is done so that calls to getSupportedOrganism 
    # don't result in that the folder is re-read multiple 
    # times
    # 

    trace.enter("scanSupportedOrganisms");
    on.exit({ trace.exit() })
    
    refdirfull = paste( refdir, "/reference_genomes/", sep="" )
    log.info( "\tSearching ", refdirfull, " for supported organisms" );
    files = dir( refdirfull )
    loc = regexpr("\\.", files)
    organisms = substr( files, 1, loc-1 )
    versions = substr( files, loc+1, nchar(files) )
    versions = versions[organisms != ""]
    organisms = organisms[organisms != ""]
    
    supportedVersions = versions
    names(supportedVersions) = organisms
    
    supportedOrganisms = sort(supportedVersions,decreasing = TRUE);
    
    setPackageVariable("supportedOrganisms" = supportedOrganisms);
    
    return(supportedOrganisms);
}

getSupportedOrganisms <- function( refdir ) {

    # uses the previously stored organisms
    # loads organisms if not yet loaded.
    #

    trace.enter("getSupportedOrganisms");
    on.exit({ trace.exit() })
    
    if (isPackageVariableDefined("supportedOrganisms")) {
        return(getPackageVariable("supportedOrganisms"));
    } else {
        log.info("getSupportedOrganisms called before scanSupportedOrganisms!");
        return(scanSupportedOrganisms(refdir));
    }
    
}

# put here the latest indexed genome/transcriptome version available
getCurrentRefVersion <- function( organism, refdir ) { # project
    trace.enter("getCurrentRefVersion");
    on.exit({ trace.exit() })
    
    supportedVersions = getSupportedOrganisms( refdir );
    version = supportedVersions[names(supportedVersions) %in% organism]
    
    if( length(version) == 0 ) {
        # register step status
        #
        log.info( "\tOrganism ", organism, " is not available on the filesystem." )
        log.info( "\tPlease use prepareReference and prepareAnnotation to get reference " )
        log.info( "\tand annotation for ", organism );
        
        registerExpStepFailure(REFERENCE_SANITY, 
                "No reference found for ", organism );
        
        stop();
    
    } else {
        version = version[1]
        log.info( "\tReference: ", version, "" )
    }
    
    #log.info( "\tUsing version ", version, " of the ", project$organism, 
    #        " ", project$reference$type, " reference" );
    return( version );
}

isOrganismSupported <- function( organism, refdir ) {
    trace.enter("isOrganismSupported");
    on.exit({ trace.exit() })
    
    supportedVersions = getSupportedOrganisms( refdir );
    supported = organism %in% names(supportedVersions);
    supported;
}
