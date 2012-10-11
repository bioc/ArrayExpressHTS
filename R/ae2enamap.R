# ENA ROUTINES

# get ftp fastQ URL base
#
getENAFastqURLBase = function(){
    trace.enter("getENAFastqURLBase");
    on.exit({ trace.exit() })

    return('ftp://ftp.sra.ebi.ac.uk/vol1');
}

# get ENA local storage path base
#
getENAFastqPathBase = function(){
    trace.enter("getENAFastqPathBase");
    on.exit({ trace.exit() })

    return('/nfs/era-pub/vol1');
}


#
# ENA EXPID -> ENA RUN
# ENA METADATA
#
# interface: HTTP
#

# create ENA xml view URL for EXPERIMENT ID
#
#
makeENAExpURL = function(enaexpid) {
    trace.enter("makeENAExpURL");
    on.exit({ trace.exit() })

    paste("http://www.ebi.ac.uk/ena/data/view/",enaexpid,"&display=xml",sep='');
}

# create url to AE experiment files
#
makeAEFilesURL = function(accession) {
    trace.enter("makeAEFilesURL");
    on.exit({ trace.exit() })

    paste("http://www.ebi.ac.uk/arrayexpress/files/",accession,"/",sep='');

    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.idf.txt
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.sdrf.txt
}

checkURLAvailable = function(url) {
    trace.enter("checkURLAvailable");
    on.exit({ trace.exit() })
    
    suppressWarnings({
        
        con.url = try(url(url, open='rb'), silent = TRUE);
        
        try.error = inherits(con.url, "try-error");
        
        try(close(con.url), silent = TRUE);
    
    })
    
    return(!try.error);
}

# create urls to AE SDRF
#
makeAESDRFURLs = function(accession) {
    trace.enter("makeAESDRFURLs");
    on.exit({ trace.exit() })
    
    return (c( paste(makeAEFilesURL(accession), accession, '.sdrf.txt', sep=''),
               paste(makeAEFilesURL(accession), accession, '.seq.sdrf.txt', sep='') ));
    
    
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.idf.txt
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.sdrf.txt
}


# check URL availability, 
# return first available SDRF URL
#
getAESDRFURL = function(accession) {
    trace.enter("getAESDRFURL");
    on.exit({ trace.exit() })
    
    # non trivial logic since there 
    # can be 2 types of SDRF files
    #
    
    sdrfurls = makeAESDRFURLs(accession);
    
    for (sdrfurl in sdrfurls) {
        if (checkURLAvailable(sdrfurl)) {
            return( sdrfurl );
        }
    }
    
    for (sdrfurl in sdrfurls) {
        log.info('SDRF ',sdrfurl, ' does not exist');
    }
    
    return (NULL);
}

# create url to AE IDF
#
makeAEIDFURL = function(accession) {
    trace.enter("makeAEIDFURL");
    on.exit({ trace.exit() })
    
    paste(makeAEFilesURL(accession), accession, ".idf.txt", sep='');
    
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.idf.txt
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.sdrf.txt
}


# check URL availability, 
# return first available SDRF URL
#
getAEIDFURL = function(accession) {
    trace.enter("getAEIDFURL");
    on.exit({ trace.exit() })
    
    # non trivial logic since there 
    # can be 2 types of SDRF files
    #
    
    idfurls = makeAEIDFURL(accession);
    
    if (checkURLAvailable(idfurls)) {
        return( idfurls );
    }
    
    log.info('IDF ',idfurls, ' does not exist');
    
    return (NULL);
}


makeENAFastqFileName = function(fname) {
    paste("fastq/", substr(fname, 1, 6), "/", substr(fname, 1, 9), "/", fname, sep="")
}


# read ENA experiment XML description
# get links to fastQ and library settings 
# for Paired-End experiments
#
readENADataForRunID = function(enarunid) {
    trace.enter("readENADataForRunID");
    on.exit({ trace.exit() })
    
    result = list();
    
    result$enarunid = enarunid;

    # read experiment info
    #
    #
    
    rundocurl = makeENAExpURL(enarunid);
    
    rundoc = xmlTreeParse(readLines(rundocurl), asText = TRUE, useInternalNodes = TRUE);
    
    expidxml = rundoc["//EXPERIMENT_REF"];
    
    if (length(expidxml) > 0)  {
        
        # FIX: incompatibility
        #
        #result$enaexpid = xmlToList(expidxml[[1]])[["accession"]];
        result$enaexpid = xmlToList( expidxml[[1]] )[[".attrs"]][["accession"]];
    
        
    } else {
        
        # unable to read library layout
        #
        # register step status
        #
        registerExpStepWarning(METADATA_SANITY, 
                enarunid, ' No Experiment ID Found in ', rundocurl);
        
        result$error$ENAExpIDNotFound = 1;
    }


    #
    # set errors beforehand
    #
    result$error$ENASampleIDNotFound = 1;
    result$error$ENAFastQFileNotFound = 1;

    # scan rundoc links
    #
    linksxml = rundoc["//RUN_LINK"];
    
    if (length(linksxml) > 0) {
        
        for(i in 1:length(linksxml)) {
        
            l = xmlToList(linksxml[[i]]);
            
            if (l$XREF_LINK$DB == "ENA-SAMPLE") {
                # store sampleid
                result$enasampleid = l$XREF_LINK$ID;
                
                # reset error
                #
                result$error$ENASampleIDNotFound = NULL;
            }
            
            if (l$XREF_LINK$DB == 'ENA-FASTQ-FILES') {
                
                descriptorurl = unlist(l$XREF_LINK$ID);
                
                cnt0025 = 3;
                while(cnt0025 >= 0) {
                    # Retries to overcome
                    # ENA interface bug
                    #
                    descriptor = try({ read.table(url(descriptorurl), header = TRUE, sep = "\t",
                                        row.names = NULL, stringsAsFactors = FALSE, fill = TRUE) });
                    
                    if (inherits(descriptor, "try-error")) {
                        sendReminder(descriptorurl, descriptor);
                        
                        if (cnt0025 == 0) {
                            #
                            #
                            registerExpStepFailure(OBTAINING_METADATA,
                                    "Failed ENA Interface: read.table(url(", descriptorurl,"))");
                            
                            stop();
                        }
                    } else {
                        #
                        # break from the loop
                        #
                        break;
                    }
                }
                
                
                indexesall = seq(descriptor[['File.Name']]);
                
                indexesin = grep(".fastq", descriptor[['File.Name']]);
                
                if (length(indexesin) != length(indexesall)) {
                    indexesout = indexesall[is.na(match(indexesall, indexesin))];
                    
                    for(i001 in indexesout) {
                        
                        registerExpStepWarning(METADATA_SANITY,
                                descriptor[['Run']][i001], ' Missing FastQ, ENA value: ', 
                                descriptor[['File.Name']][i001]);
                        
                    }
                }
                
                if (length(indexesin) > 0) {
                    # add only if there are fastq files
                    #
                    result$enafastq = makeENAFastqFileName(descriptor[['File.Name']][indexesin]);
                    
                    # reset the error
                    #
                    result$error$ENAFastQFileNotFound = NULL;
                }
                
            }
        }
    }
    
    # produce warning log
    #
    if (!is.null(result$error$ENASampleIDNotFound)) {
        
        registerExpStepWarning(METADATA_SANITY,
                enarunid, ' No Sample ID Found in ', rundocurl);
    }


    # produce warning log
    #
    if (!is.null(result$error$ENAFastQFileNotFound)) {
        registerExpStepWarning(METADATA_SANITY,
                enarunid, ' No FastQ File Found in ', rundocurl);
    }
    
    if (is.null(result$error$ENAExpIDNotFound)) {
        
        # get experiment descriptor
        #
        #
        expdocurl = makeENAExpURL(result$enaexpid);
        
        expdoc = xmlTreeParse(readLines(expdocurl), asText = TRUE, useInternalNodes = TRUE);
        
        # read layout info
        #
        #
        layoutxml = expdoc["//LIBRARY_LAYOUT"];
        
        if (length(layoutxml) > 0) {
            
            layout = xmlToList(layoutxml[[1]]);
            
            layoutname = names(layout);
            
            result$layout = list();
            
            result$layout$name = layoutname;
            
            if (layoutname == "PAIRED") {
                
                result$layout$info = layout$PAIRED;
                #result$layout$NOMINAL_LENGTH = layout$PAIRED['NOMINAL_LENGTH']
                #result$layout$ORIENTATION = layout$PAIRED['ORIENTATION']
                #result$layout$NOMINAL_SDEV = layout$PAIRED['NOMINAL_SDEV']
            }
        } else {
            registerExpStepWarning(METADATA_SANITY,
                    enarunid, " No Library Layout found in ", expdocurl);
            
            # set error
            #
            result$error$ENALibraryLayoutNotFound = 1;
        }
    }

    if (is.null(result$error$ENASampleIDNotFound)) {
    
        # read sample info
        #
        #
            
        sampledocurl = makeENAExpURL(result$enasampleid)
        
        sampledoc = xmlTreeParse(readLines(sampledocurl), asText = TRUE, useInternalNodes = TRUE);
        
        taxonxml = sampledoc["//TAXON_ID"];
        
        if (length(taxonxml) > 0) {
        
            result$sampletaxon = xmlToList(taxonxml[[1]]);
            
        } else {
            
            registerExpStepWarning(METADATA_SANITY,
                    enarunid, " No Taxon ID found in ", sampledocurl);
            
            # taxon not found
            # can't think of any case where it's really
            # used (so far) hence no error for this case
            # 
        }
        
        # get common sample name
        #
        commonnamexml = sampledoc["//COMMON_NAME"];
        
        if (length(commonnamexml) > 0) {
            # sample name found
            #
            result$samplename = xmlToList(commonnamexml[[1]]);
        } else {
            
            # try to get scientific sample name
            #
            commonnamexml = sampledoc["//SCIENTIFIC_NAME"];
            
            if (length(commonnamexml) > 0) {
                # found it
                #
                result$samplename = xmlToList(commonnamexml[[1]]);
                
            } else {
                # bad news
                #
                registerExpStepWarning(METADATA_SANITY,
                        enarunid, ' No Sample Name Found in ', sampledocurl);
                
                result$error$ENASampleNameNotFound = 1;
            }
        }
    }
    
    return(result);
    
}


# read ENA experiment XML description
# get links to fastQ and library settings 
# for Paired-End experiments
#
readENADataForExpID = function(enaexpid) {
    trace.enter("readENADataForExpID");
    on.exit({ trace.exit() })
    
    result = list();
    
    # read experiment info
    #
    #
    
    expdocurl = makeENAExpURL(enaexpid);
    
    expdoc = xmlTreeParse(readLines(expdocurl), asText = TRUE, useInternalNodes = TRUE);
    
    links = expdoc["//EXPERIMENT_LINK"];
    
    cnt = length(links);
    
    if (cnt > 0) {
        for(i in 1:cnt) {
            l = xmlToList(links[[i]]);
            
            if (l$XREF_LINK$DB == 'ENA-RUN') {
                result$enarunid = l$XREF_LINK$ID;
            }
            
            if (l$XREF_LINK$DB == 'ERA-FASTQ') { ## ENA backward compatible field name
                result$enafastq = unlist(strsplit(l$XREF_LINK$ID, ","));
            
            } else if (l$XREF_LINK$DB == 'ENA-FASTQ-FILES') { ## new interface
                
                descriptorurl = unlist(l$XREF_LINK$ID);
                
                cnt0025 = 3;
                while(cnt0025 >= 0) {
                    # Retries to overcome
                    # ENA interface bug
                    #
                    descriptor = try({ read.table(url(descriptorurl), header = TRUE, sep = "\t", 
                                        row.names = NULL, stringsAsFactors = FALSE, fill = TRUE) });
                    
                    if (inherits(descriptor, "try-error")) {
                        sendReminder(descriptorurl, descriptor);
                        
                        if (cnt0025 == 0) {
                            #
                            #
                            registerExpStepFailure(OBTAINING_METADATA,
                                    "Failed ENA Interface: read.table(url(", descriptorurl,"))");
                            
                            stop();
                        }
                    } else {
                        #
                        # break from the loop
                        #
                        break;
                    }
                }
                
                indexesall = seq(descriptor[['File.Name']]);
                
                indexesin = grep(".fastq", descriptor[['File.Name']]);
                
                if (length(indexesin) != length(indexesall)) {
                    
                    indexesout = indexesall[is.na(match(indexesall, indexesin))];
                    
                    for(i001 in indexesout) {
                        registerExpStepWarning(METADATA_SANITY,
                                descriptor[['Run']][i001], ' Missing FastQ, ENA value: ', 
                                descriptor[['File.Name']][i001]);
                    
                    }
                }
                
                if (length(indexesin) > 0) {
                    # add only if there are fastq files
                    #
                    result$enafastq = makeENAFastqFileName(descriptor[['File.Name']][indexesin]);
                }
                
            }
        }
    } else {
        registerExpStepWarning(METADATA_SANITY,
                enaexpid, ' No FastQ links found in ', expdocurl);
    }
    
    
    # read layout info
    #
    #
    
    # read layout info
    layoutxml = expdoc["//LIBRARY_LAYOUT"];
    
    if (length(layoutxml) > 0) {
        
        layout = xmlToList(layoutxml[[1]]);
        
        layoutname = names(layout);
        
        result$layout = list();
        
        result$layout$name = layoutname;
    
        if (layoutname == "PAIRED") {
            
            result$layout$info = layout$PAIRED;
            #result$layout$NOMINAL_LENGTH = layout$PAIRED['NOMINAL_LENGTH']
            #result$layout$ORIENTATION = layout$PAIRED['ORIENTATION']
            #result$layout$NOMINAL_SDEV = layout$PAIRED['NOMINAL_SDEV']
        }
    } else {
        registerExpStepWarning(METADATA_SANITY,
                enaexpid, ' No Library Layout found in ', expdocurl);
    }
    
    
    # read sample info
    #
    #
    samplexml = expdoc["//SAMPLE_DESCRIPTOR"];
    
    if (length(samplexml) > 0) {
        
        # FIX: incompatibility
        #
        #sampleid = xmlToList(samplexml[[1]])[['accession']];
        sampleid = xmlToList(samplexml[[1]])[[".attrs"]] [["accession"]]
        
        sampledocurl = makeENAExpURL(sampleid)
        
        sampledoc = xmlTreeParse(readLines(sampledocurl),asText = TRUE,useInternalNodes = TRUE);
    
        taxonxml = sampledoc['//TAXON_ID'];
        
        if (length(taxonxml) > 0) {
            result$sampletaxon = xmlToList(taxonxml[[1]]);
        } else {

            registerExpStepWarning(METADATA_SANITY,
                    enaexpid, ' No Taxon ID found in ', sampledocurl);
        }
        
        # get common sample name
        #
        commonnamexml = sampledoc['//COMMON_NAME'];
    
        if (length(commonnamexml) > 0) {
            result$samplename = xmlToList(commonnamexml[[1]]);
        
        } else {
        
            # get scientific sample name
            #
            commonnamexml = sampledoc['//SCIENTIFIC_NAME'];
            
            if (length(commonnamexml) > 0) {
                
                result$samplename = xmlToList(commonnamexml[[1]]);
            
            } else {
                registerExpStepWarning(METADATA_SANITY,
                        enaexpid, ' No Sample Name found in ', sampledocurl);
            }
        }
    
    
    } else {
        registerExpStepWarning(METADATA_SANITY,
                enaexpid, ' No Sample ID found in ', expdocurl);
    }
    
    return(result);
    
}

#
# AE -> ENA EXPID 
#
#  G E T   A C C E S S I O N
#
#

# get public AE path base 
#
#
getAEExperimentPathBase = function() {
    trace.enter("getAEExperimentPathBase");
    on.exit({ trace.exit() })
    
    #/ebi/ftp/pub/databases/microarray/data/experiment/MEXP/E-MEXP-511
    '/ebi/ftp/pub/databases/microarray/data/experiment';
}

# create SDRF file name from accession,
# check is file exists.
#
#
makeAESDRFFilenames = function(accession) {
    trace.enter("makeAESDRFFilenames");
    on.exit({ trace.exit() })
    
    # non trivial logic since there 
    # can be 2 types of SDRF files
    #
    paste(accession,'.sdrf.txt', sep='')
    
    filenamebase = paste(getAEExperimentPathBase(), substr(accession,3,6),accession, sep='/');
    
    return( c(paste(filenamebase,paste(accession,'.sdrf.txt', sep=''), sep='/'),
              paste(filenamebase,paste(accession,'.seq.sdrf.txt', sep=''), sep='/')) );
        
}


# create SDRF file name from accession,
# check is file exists.
#
#
getAESDRFFilename = function(accession) {
    trace.enter("getAESDRFFilename");
    on.exit({ trace.exit() })
    
    # there can be 2 types of SDRF files
    #
    
    fnames = makeAESDRFFilenames(accession);
    
    for(sdrffname in fnames) {
        if(file.exists(sdrffname)) {
            return (sdrffname);
        }
    }
    
    for(sdrffname in fnames) {
        log.info('SDRF ',sdrffname, ' does not exist');
    }
    
    return(NULL);
}


# get raw data filr names for AE experiment 
# if any raw data files exist
#
#
getAERawDataFilenames = function(accession) {
    trace.enter("getAERawDataFilenames");
    on.exit({ trace.exit() })
    
    folder = paste(getAEExperimentPathBase(), substr(accession,3,6), accession, sep='/');
    return(dir(folder, pattern=".raw", full.names = TRUE));
}

# create IDF file name from accession,
# check is file exists.
#
#
getAEIDFFilename = function(accession) {
    trace.enter("getAEIDFFilename");
    on.exit({ trace.exit() })
    
    filename = paste(getAEExperimentPathBase(),
    substr(accession,3,6),accession,paste(accession,'.idf.txt', sep=''), sep='/');
    
    if (file.exists(filename)) {
        return( filename );
    }
     
    log.info('IDF ',filename, ' does not exist');
    
    return(NULL);
}

# read SDRF
#
#
readSDRF = function(fname) {
    trace.enter("readSDRF");
    on.exit({ trace.exit() })
    
    sdrf = read.AnnotatedDataFrame(fname, path = "", 
    row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char="$");
}


# read SDRF field
#
getDataFromSDRF = function(accession, sdrf, record) {
    trace.enter("getDataFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        
        indexes0 = getSupportedStrategyRowIndexes(accession, sdrf, getPipelineOption("supportedLibrary"));
        #indexes0 = getSupportedSelectionRowIndexes(accession, sdrf, getPipelineOption("supportedSelection"));
        
        if (is.null(indexes0)) {
            return(record);
        }
        
        record$assayname = c();
        record$sdrfrunid = c();
        record$sdrfexpid = c();
        record$sdrforganism = c();
        
        assays   = getAssayNameFromSDRF(sdrf)[ indexes0 ];
        runid    = getRunIDFromSDRF(sdrf)[ indexes0 ];
        organism = getOrganismFromSDRF(sdrf)[ indexes0 ];
        
        index = grep("ENA_EXPERIMENT", names(sdrf@data));
        
        if (length(index) > 0) {
            for( i in 1:length(indexes0) ) {
                expidrow = sdrf@data[indexes0[[ i ]], index];
                expidrow = expidrow[ !is.na(expidrow) & expidrow != "" ]; 
                
                if (length(expidrow) > 0) {
                    for(expid in expidrow) {
                        record$assayname = c(record$assayname, assays[[i]]);
                        record$sdrfrunid = c(record$sdrfrunid, runid[[i]]);
                        record$sdrfexpid = c(record$sdrfexpid, expid);
                        record$sdrforganism = c(record$sdrforganism, organism[[i]]);
                    }
                }
            }
        }
    }
    
    return(record);
}


# read ENA accessions (experiments) 
# from certain SDRF field
#
#
getENAExpIDFromSDRF = function(sdrf) {
    trace.enter("getENAExpIDFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("ENA_EXPERIMENT", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }
}


# strip off quotas
#
#
stripOffQuotas = function(stringarray) {
    for(i in 1:length(stringarray)) {
        if (length(grep("\"", stringarray[i])) > 0) {
            stringarray[i] = substring(stringarray[i], 2, nchar(stringarray[i]) -1);
        }
    }
    return (stringarray);
}

# read organism from SDRF
#
#
getOrganismFromSDRF = function(sdrf) {
    trace.enter("getOrganismFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        #index = grep("Characteristics.Organism", names(sdrf@data));
        index = grep(".Organism.", names(sdrf@data));
        if (length(index) > 0) {
            
            return( stripOffQuotas( sdrf@data[[ sort(names(sdrf@data)[index])[1] ]] ) );
        }
    }
    
    return(NULL);
}

# read LIBRARY_STRATEGY
#
#
getLibraryStrategyFromSDRF = function(sdrf) {
    trace.enter("getLibraryStrategyFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("LIBRARY_STRATEGY", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }
    
    return(NULL);
}

# read LIBRARY_SELECTION
#
#
getLibrarySelectionFromSDRF = function(sdrf) {
    trace.enter("getLibrarySelectionFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("LIBRARY_SELECTION", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }
    
    return(NULL);
}


# read FactorValue
#
#
getFactorValueFromSDRF = function(sdrf) {
    trace.enter("getFactorValueFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("Factor", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }
    
    return(NULL);
}


# read Assay Name
#
#
getAssayNameFromSDRF = function(sdrf) {
    trace.enter("getAssayNameFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("Assay.Name", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }
    
    return(NULL);
}

# read Run ID
#
#
getRunIDFromSDRF = function(sdrf) {
    trace.enter("getRunIDFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("ENA_RUN", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ index[1] ]]);
        }
    }

    return(NULL);
}


# work out SDRF filename, read SDRF and get
# ENA experiment IDs from it
#
#
getENAExpIDForAEAcc = function(accession) {
    trace.enter("getENAExpIDForAEAcc");
    on.exit({ trace.exit() })
    
    sdrffname = getAESDRFFilename(accession);
    
    #log.info('reading SDRF..');
    
    sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
    if (inherits(sdrf0, 'try-error')) {
        log.info('Error reading SDRF ', sdrffname);
        
        return(NULL);
    } else {
        #log.info('getting EXP IDs..');
        
        record = getENAExpIDFromSDRF(sdrf0);
        
        return(record);
    }
}

# work out SDRF filename, 
# read SDRF and return it
#
#
getSDRFforAccession = function(accession) {
    trace.enter("getSDRFforAccession");
    on.exit({ trace.exit() })
    
    sdrffname = getAESDRFFilename(accession);
    
    #log.info('reading SDRF..');
    
    sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
    if (inherits(sdrf0, 'try-error')) {
        log.info('Error reading SDRF ', sdrffname);
        
        return(NULL);
    } else {
        #log.info('getting EXP IDs..');
        
        return(sdrf0);
    }
}



#
# get ENA experiment id and organism from SDRF 
# and for each experiment obtain library meta data
# and RUN IDs
# 
# interface : HTTP
#
mapAEtoENAviaHTTP = function(accession, sdrffname) {
    trace.enter("mapAEtoENAviaHTTP");
    on.exit({ trace.exit() })
    
    #sdrffname = getAESDRFFilename(accession);
    
    record = list(accession=accession, sdrffname=sdrffname);
    
    if (!file.exists(sdrffname)) {
        
        # register step status
        #
        registerExpStepFailure(OBTAINING_METADATA,
                accession, " Error: SDRF not found ", sdrffname);
        
        stop();
        
        record$error$SDRFNotFound = 1;
        
    } else {
        sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
        if (inherits(sdrf0, 'try-error')) {
        
            # register step status
            #
            registerExpStepFailure(OBTAINING_METADATA,
                    accession, " Failed to read SDRF ", sdrffname);
            
            stop();
            
            record$error$SDRFNotReadable = 1;
           
        } else {
            #log.info('getting EXP IDs..');
            
            record = getDataFromSDRF(accession, sdrf0, record);
            
            #record$assayname = getAssayNameFromSDRF(sdrf0)
            #record$sdrfrunid = getRunIDFromSDRF(sdrf0);
            #record$sdrfexpid = getENAExpIDFromSDRF(sdrf0);
            #record$sdrforganism = getOrganismFromSDRF(sdrf0);
            
            
            # try to read data from ENA for each RUNIDs
            #
            #
            if (!is.null(record$sdrfrunid)) {
                
                if (!all(is.na((record$sdrfrunid)))) {
                    # check if there are empty records
                    # among the RUNID values
                    #
                    if (any(record$sdrfrunid == "")) {
                        
                        # bad news
                        #
                        # register status
                        #
                        registerExpStepWarning(METADATA_SANITY,
                                accession," Missing RUNIDs in SDRF");
                        
                        record$error$SDRFSomeRunIDMissing = 1;
                    }
                    
                    for (runid in unique(record$sdrfrunid)) {
                        if (nchar(runid) > 0) {
                            #log.info('   x=-', x,'-');
                            if (is.null(record$enarunid[[ runid ]])) {
                                record$enarunid[[ runid ]] = readENADataForRunID( runid );
                            }
                        }
                    }
                } else {
                    # bad news
                    #
                    registerExpStepWarning(METADATA_SANITY,
                            accession," RUNID Section in SDRF is NA");
                    
                    record$error$SDRFRunIDNotFound = 1;
                }
                
            } else {
                # bad news
                #
                registerExpStepWarning(METADATA_SANITY,
                        accession," No RUNID Section in SDRF");
                
                record$error$SDRFRunIDNotFound = 1;
            }
            
            
            #log.info('getting RUN IDs..');
            if (!is.null(record$sdrfexpid)) {
            
                if (!all(is.na((record$sdrfexpid)))) {
                    # check if there are empty records
                    # among the EXPID values
                    #
                    if (any(record$sdrfexpid == "")) {
                        # bad news
                        #
                        registerExpStepWarning(METADATA_SANITY,
                                accession," Missing EXPIDs in SDRF");
                        
                        record$error$SDRFSomeExpIDMissing = 1;
                    }
                    
                    for (x in unique(record$sdrfexpid)) {
                        if (nchar(x) > 0) {
                            #log.info('   x=-', x,'-');
                            if (is.null(record$enaexpid[[x]])) {
                                record$enaexpid[[x]] = readENADataForExpID(x);
                                record$enaexpid[[x]]$enaexpid = x;
                            }
                        }
                    }
                
                } else {
                    # bad news
                    #
                    registerExpStepWarning(METADATA_SANITY,
                            accession," EXPID Section in SDRF is NA");
                    
                    record$error$SDRFExpIDNotFound = 1;
                }
            
            } else {
                # bad news
                #
                registerExpStepWarning(METADATA_SANITY,
                        accession," No EXPID Section in SDRF");

                record$error$SDRFExpIDNotFound = 1;
            }
        }
    }
    
    
    return(record);
}


#
# get SDRF rows that have "RNA-Seq" in the LIBRARY_STRATEGY
# 
getSupportedStrategyRowIndexes = function(accession, sdrf0, supportedLibraries) {
    trace.enter("getSupportedStrategyRowIndexes");
    on.exit({ trace.exit() })
    
    # get library strategies
    #
    strategies0 = getLibraryStrategyFromSDRF( sdrf0 );
    
    # check if the library is defined
    
    if (!is.null(strategies0)) {
        
        indexes0 = c(); 
        
        for(lib0 in supportedLibraries) { 
            indexes0 = c(indexes0, grep(tolower( lib0 ), tolower( strategies0 ))); 
        }
        
        #indexes0 = grep(tolower("RNA-Seq"), tolower( strategies0 ));
        
        if (length( indexes0 ) > 0) {
            return( indexes0 );
            
        } else {
            log.info(accession, "No supported LIBRARY_STRATEGY found");
            
            for(s0 in strategies0) {
                log.info("\t", s0);
            }
            
            registerExpStepFailure(METADATA_SANITY,
                    accession, " No supported LIBRARY_STRATEGY found");
            
            stop();
            # return(NULL);
        }
    }
    
    registerExpStepWarning(METADATA_SANITY,
            accession, " LIBRARY_STRATEGY is not defined in SDRF");
    
    
    stop();
    
    # return(NULL);
}

#
# get SDRF rows that have "cDNA" in the LIBRARY_SELECTION
# 
getSupportedSelectionRowIndexes = function(accession, sdrf0, supportedSelection) {
    trace.enter("getSupportedSelectionRowIndexes");
    on.exit({ trace.exit() })
    
    # get library strategies
    #
    selections0 = getLibrarySelectionFromSDRF( sdrf0 );
    
    # check if the library is defined
    
    if (!is.null(selections0)) {
        
        indexes0 = c(); 
        
        for(lib0 in supportedSelection) { 
            indexes0 = c(indexes0, grep(tolower( lib0 ), tolower( selections0 ))); 
        }
        
        #indexes0 = grep(tolower("RNA-Seq"), tolower( selections0 ));
        
        if (length( indexes0 ) > 0) {
            return( indexes0 );
            
        } else {
            log.info(accession, "No supported LIBRARY_SELECTION found");
            
            for(s0 in selections0) {
                log.info("\t", s0);
            }
            
            registerExpStepFailure(METADATA_SANITY,
                    accession, " No supported LIBRARY_SELECTION found");
            
            stop();
            # return(NULL);
        }
    }
    
    registerExpStepWarning(METADATA_SANITY,
            accession, " LIBRARY_SELECTION is not defined in SDRF");
    
    
    stop();
    
    # return(NULL);
}


# FUNCTIONS
#
#

#
# read table
# used to reasd ENA dada file
# 
readData <- function(filename){
    trace.enter("readData");
    on.exit({ trace.exit() })
    
    data = read.table( filename, header = TRUE, stringsAsFactors = FALSE, 
        fill = TRUE, row.names = NULL, blank.lines.skip = TRUE, sep="\t" );
    
    data;
}

# check whether the resource
# denoted in the URL is available
#
checkURLResourcesAvailable = function(urls){ 
    trace.enter("checkURLResourcesAvailable");
    on.exit({ trace.exit() })
    
    retcode = c();
    
    for (x in urls) {
        
        url0 = try(url(x));
        
        if (inherits(url0, 'try-error')) {
            
            retcode = c(retcode, F);
        } else {
            
            result = try( readLines(con = url0, n = 2), silent = TRUE )
            
            if ( inherits(result, 'try-error') || length(result) != 2 ) {
                #warning(paste(retcode, " does not exist or is empty."),sep="");
                
                retcode = c(retcode, F);
            } else {
                retcode = c(retcode, T);
            }
            
            close(url0);
        }
    }
    
    retcode;
}

makeENAFastqURL = function(run, paired){
    trace.enter("makeENAFastqURL");
    on.exit({ trace.exit() })
    
    urlbase = paste(getENAFastqURLBase(), 'fastq', substr(run, 1, 6), run, run, sep='/')
    
    if (paired) {
        
        url0 = paste(urlbase, '_1.fastq.gz', sep='');
        url1 = paste(urlbase, '_2.fastq.gz', sep=''); 
        c(url0, url1);
        
    } else {
        url0 = paste(urlbase, '.fastq.gz', sep='');
        c(url0);
    }
}

#
# ENA DATA FILE
#
#
#

# download and ungzip ENA data file
#
#

downloadENADataFile = function(){
    trace.enter("downloadENADataFile");
    on.exit({ trace.exit() })
    
    #download.
    enadatafilename = 'fastqFileReport.gz';
    enadatafileurl = 'ftp://ftp.sra.ebi.ac.uk/meta/list/fastqFileReport.gz';
    
    log.info('downloading.. ');
    
    result = try(download.file(enadatafileurl, enadatafilename, mode = "wb"))
    
    if (inherits(result, 'try-error')) {
        log.info('error downloading ', enadatafileurl, ' ', result)
        
    } else {
        log.info('unzipping..');
        
        system(paste('gunzip ', getwd(), '/', enadatafilename, sep=''))
        #zip.file.extract();
    }
}

# load ENA DATA
#
#

loadENAData = function() {
    trace.enter("loadENAData");
    on.exit({ trace.exit() })
    
    enadata = readData('fastqFileReport');
    enadata;
}

# filter out "not ILLUNIMA" records
# filter out "no fastQ" records
#
#
filterENAData = function(inputdata){ 
    trace.enter("filterENAData");
    on.exit({ trace.exit() })
    
    omit = FALSE;
    omit = omit | (inputdata$INSTRUMENT_PLATFORM != "ILLUMINA");
    filtered0 = inputdata[!omit,];
    
    log.info('filtered out ', sum(omit),' experiments of not ILLUMINA type');
    
    omit = FALSE;
    omit = omit | (filtered0$FASTQ_FILES == 0);
    filtered1 = filtered0[!omit,];
    
    log.info('filtered out ', sum(omit),' experiments not having fastq files');
    
    paired = FALSE;
    paired = paired | (filtered1$FASTQ_FILES == 2);
    log.info('paired: ', sum(paired),' experiments');
    log.info('single: ', sum(!paired),' experiments');
    
    filtered1;
    #ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456/
} 

# GET ENA DATA
#
# interface: ENA DATA FILE
#

# compose a record of fields of interest from an ENA data record, 
# create fastQ URLs, check if fastq files are available, check layout.
#
# input: a "record" from ENA data frame.
# 
getENAAccessionRecord = function(run) {
    trace.enter("getENAAccessionRecord");
    on.exit({ trace.exit() })
    
    runid = run$RUN_ID
    
    numfiles = run$FASTQ_FILES;
    
    layout = run$LIBRARY_LAYOUT;
    
    organism = run$COMMON_NAME;
    
    layoutvalid = TRUE;
    
    if (numfiles == 0) {
        log.info(runid,': NO FASTQ for ', runid);
    }
    
    if ((numfiles == 2 && layout == "SINGLE") || (numfiles == 1 && layout == "PAIRED")) {
        log.info(runid,': LAYOUT MISMATCH: numfiles = ', numfiles, ' layout = ', layout);
        layoutvalid = FALSE;
    }
    
    paired = (numfiles == 2);
    
    urls = makeENAFastqURL(runid, paired);
    
    #urlvalid = checkURLResourcesAvailable(urls);
    urlvalid = c(T,T);
    
    for(i in 1:length(urlvalid)) {
        if (!urlvalid[i]) {
            log.info(runid,': FASTQ URL is invalid, url = ', urls[i]);
        }
    }
    
    list(runid=runid, organism=organism, numfiles=numfiles, urls=urls, urlvalid=urlvalid, layoutvalid=layoutvalid);
    
    #data.frame(urls=urls, urlvalid=urlvalid, layoutvalid=layoutvalid);
    #list(runid=runid, numfiles=numfiles, urls=urls, valid=valid, organism=organism);
}

#
# do getENAAccessionRecord for all RUNS in ENA data frame
#
# getENAAccessionRecord checks consistency and returns validity
# info, which can be used to estimate number of valid runs.
#
scanAllENAData = function(data, max=400000) {
    trace.enter("scanAllENAData");
    on.exit({ trace.exit() })
    
    #result = NULL; #data.frame();
    
    cnt = nrow(data);
    
    result = vector("list", cnt);
    
    p = floor(cnt / 100);
    curp = p;
    
    for (x in 1:cnt) {
        
        exp = data[x,];
        
        record = getENAAccessionRecord(exp);
        
        #if (is.null(result)) {
            #    result = record;
        #} else {
            #    result = rbind(result, record);
        #}
        
        result[[x]] = record;
        
        if (x == curp) {
            log.info(curp / p, "%% ",x," out of ",cnt," complete.." );
            curp = curp + p;
        }
        
        if (x > max) {
            log.info('finishing...');
            return(result);
        }
    }
    
    result;
}

# do getENAAccessionRecord for all RUNS that 
# correspond to specific ENA EXPERIMENT ID
#

getENARunIDForExpID = function(enadata, expid) {
    trace.enter("getENARunIDForExpID");
    on.exit({ trace.exit() })
    
    cnt = nrow(enadata);
    result = list();
    count = 1;
    
    for(i in 1:cnt) { 
        if (enadata$EXPERIMENT_ID[[i]] == expid) {
            exp = enadata[i,];
            result[[count]] = getENAAccessionRecord(exp);
            count = count + 1;
        }
    };
    result;
}



#
#   G R E P   A C C E S S I O N S
#

# ENA ID prefixes for grepping SDRF fields 
#
#
getENAAccessionPrefixes = function() {
    trace.enter("getENAAccessionPrefixes");
    on.exit({ trace.exit() })
    
    c('SRX', 'ERX', 'DRX', 'SRP', 'ERP', 'DRP', 'SRR', 'ERR', 'DRR');
}

# grep all SDRF fields to find any ENA
# references, return all found
#
grepSDRFForENAAccessions = function(sdrf, regexp = getENAAccessionPrefixes()) {
    trace.enter("grepSDRFForENAAccessions");
    on.exit({ trace.exit() })
    
    result = list();
    colnames = names(sdrf@data);
    cnt = ncol(sdrf);
    
    
    for(i in 1:cnt) {
        record = sdrf[[i]];
        #record = paste(record, collapse=' ');
        
        for (exp in regexp) {
            match = any(grep(exp, record)); 
            if (!is.na(match) && match) {
                result[[colnames[i]]] = record;
            }
        }
    }
    
    result;
}

# similar to getENAExpIDForAEAcc, but greps all 
# possible references to ENA from AE SDRF
#
#
grepAllENAAccessions = function(accession) {
    trace.enter("grepAllENAAccessions");
    on.exit({ trace.exit() })
    
    sdrffname = getAESDRFFilename(accession);
    
    #log.info('reading SDRF..');
    
    sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
    if (inherits(sdrf0, 'try-error')) {
        log.info('Error reading SDRF ', sdrffname);
        
        return(NULL);
    } else {
        #log.info('getting EXP IDs..');
        
        record = grepSDRFForENAAccessions(sdrf0);
        
        if (is.null(record)) {
            return(list(accession=accession, sdrffname=sdrffname));
        } else {
            return(list(accession=accession, sdrffname=sdrffname, record=record));
        }
    }
}

#
# MAP AE -> ENA 
#
# get ENA experiment id and organism from SDRF 
# and for each experiment obtain library meta data
# and RUN IDs
# 
# interface: ENA DATA FILE
#
mapAEtoENAviaENADataFile = function(enadata, accession) {
    trace.enter("mapAEtoENAviaENADataFile");
    on.exit({ trace.exit() })
    
    sdrffname = getAESDRFFilename(accession);
    
    record = list(accession=accession, sdrffname=sdrffname);
    
    sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
    if (inherits(sdrf0, 'try-error')) {
        log.info('Error reading SDRF ', sdrffname);
    } else {
        log.info('getting EXP IDs..');
        
        enaexpid = getENAExpIDFromSDRF(sdrf0);
        organism = getOrganismFromSDRF(sdrf0);
        
        record$organism = organism;
        
        log.info('getting RUN IDs..');
        
        cnt = length(enaexpid);
        
        if (cnt > 0) {
            #for (i in 1:cnt) 
            record$enarunid = list();
            
            for (x in enaexpid) {
                log.info('x=',x);
                record$enaexpid[[x]] = getENARunIDForExpID(enadata, x);
            }
        }
    }
    
    return(record);
}

#
# GET ALL AE ACCESSIONS
#
#

#
# GET ALL HTS + RNA Seq ACCESSIONS
#
#

getAllHTSAccessionsFromAE = function(refresh = FALSE) {
    trace.enter("getAllHTSAccessionsFromAE");
    on.exit({ trace.exit() })
    
    queryfilename = paste("query", ".HTS.All", ".xml", sep = "")
    
    if (!file.exists(queryfilename) || refresh) {
        
        query = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?",
        "exptype[]=%22high+throughput+sequencing+assay%22",sep='');
        
        query = try(download.file(query, queryfilename, mode = "wb"))
    }
    
    x = xmlTreeParse(queryfilename);
    
    root = xmlRoot(x);
    cnt = length(root);
    per = floor(cnt / 10);
    cur = seq(per, per * 12, per);;
    
    log.info('getting accessions from XML .. takes a while' );
    
    accession = sapply(1:cnt, function(i) { 
        if (i %in% cur) {
            log.info( floor(i/per*10), '%% ', i, ' out of ', cnt, ' ..' );
        }
        
        unlist(xmlElementsByTagName(root[[i]], 
        "accession"))["accession.children.text.value"]; 
    });
    
    return(accession);
}


getHTSAccessionsFromAE = function(refresh = FALSE) {
    trace.enter("getHTSAccessionsFromAE");
    on.exit({ trace.exit() })
    
    queryfilename = paste("query", "HTS", ".xml", sep = "")
    
    if (!file.exists(queryfilename) || refresh) {
        
        #query = paste("http://www.ebi.ac.uk/microarray-as/ae/xml/v2/experiments",
        #"?keywords=&species=&array=&exptype[]=",
        #"&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        #"&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        query = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments",
        "?keywords=&species=&array=&exptype[]=",
        "&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        "&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        query = try(download.file(query, queryfilename, mode = "wb"))
    }
    
    x = xmlTreeParse(queryfilename);
    
    root = xmlRoot(x);
    cnt = length(root);
    per = floor(cnt / 10);
    cur = seq(per, per * 12, per);;
    
    log.info('getting accessions from XML .. takes a while' );
    
    accession = sapply(1:cnt, function(i) { 
        if (i %in% cur) {
            log.info( floor(i/per*10), '%% ', i, ' out of ', cnt, ' ..' );
        }
        
        unlist(xmlElementsByTagName(root[[i]], 
        "accession"))["accession.children.text.value"]; 
    });
    
    return(accession);
}

#
# GET ALL HTS + RNA Seq ACCESSIONS
#
#

getHTSRNASeqAccessionsFromAE = function(refresh = FALSE) {
    trace.enter("getHTSAccessionsFromAE");
    on.exit({ trace.exit() })
    
    queryfilename = paste("query", ".HTS.RNASeq", ".xml", sep = "")
    
    if (!file.exists(queryfilename) || refresh) {
        
        #query = paste("http://www.ebi.ac.uk/microarray-as/ae/xml/v2/experiments",
        #"?keywords=&species=&array=&exptype[]=",
        #"&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        #"&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        #query = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments",
        #"?keywords=&species=&array=&exptype[]=",
        #"&exptype[]=RNA+assay&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        #"&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        query = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?",
        "exptype[]=%22RNA+assay%22&exptype[]=%22high+throughput+sequencing+assay%22",sep='');
        
        #http://www.ebi.ac.uk/arrayexpress/browse.html?keywords=&species=Homo+sapiens&array=
        #&exptype[]=RNA+assay&exptype[]=high+throughput+sequencing+assay&pagesize=25&sortby=releasedate&sortorder=descending&expandefo=on
        
        query = try(download.file(query, queryfilename, mode = "wb"))
    }
    
    x = xmlTreeParse(queryfilename);
    
    root = xmlRoot(x);
    cnt = length(root);
    per = floor(cnt / 10);
    cur = seq(per, per * 12, per);;
    
    log.info('getting accessions from XML .. takes a while' );
    
    accession = sapply(1:cnt, function(i) { 
        if (i %in% cur) {
            log.info( floor(i/per*10), '%% ', i, ' out of ', cnt, ' ..' );
        }
        
        unlist(xmlElementsByTagName(root[[i]], 
        "accession"))["accession.children.text.value"]; 
    });
    
    return(accession);
}

getAEExpType = function(accession) {
    trace.enter("getAEExpType");
    on.exit({ trace.exit() })
    
    expurl = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments/", accession, sep="");
    
    suppressWarnings({
                expdoc = try(xmlTreeParse(readLines(expurl), asText = TRUE, 
                                useInternalNodes = TRUE), silent = TRUE);
                });
    
    if (inherits(expdoc, 'try-error')) {
        log.info(expdoc);
        return(NULL);
    } else {
        links = expdoc["//experimenttype"];
        return(sapply(links, function(x){ xmlToList(x)[[1]] } ));
    }
}


getAEExpDescription = function(accession) {
    trace.enter("getAEExpDescription");
    on.exit({ trace.exit() })
    
    expurl = paste("http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments/", accession, sep="");
    
    expdoc = try(xmlTreeParse(readLines(expurl), asText = TRUE, useInternalNodes = TRUE), silent = TRUE);
    
    if (inherits(expdoc, 'try-error')) {
        message("Error!");
        return(NULL);
    } else {
        desc = expdoc["//description"];
        #return(sapply(desc, function(x){ xmlToList(x)$text } ));
        return(sapply(desc, function(x){ xmlValue(x) } ));
    }
}

