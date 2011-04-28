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

# create url to AE SDRF
#
makeAESDRFURL = function(accession) {
    trace.enter("makeAESDRFURL");
    on.exit({ trace.exit() })
    
    paste(makeAEFilesURL(accession), accession, ".sdrf.txt", sep='');
    
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.idf.txt
    # http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-16190/E-GEOD-16190.sdrf.txt
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
            }
        }
    } else {
        log.info(enaexpid, ' Warning: No FastQ links found in ', expdocurl);
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
        log.info(enaexpid, ' Warning: No Library Layout found in ', expdocurl);
    }
    
    
    # read sample info
    #
    #
    samplexml = expdoc["//SAMPLE_DESCRIPTOR"];
    
    if (length(samplexml) > 0) {
        sampleid = xmlToList(samplexml[[1]])[['accession']];
        
        sampledocurl = makeENAExpURL(sampleid)
        
        sampledoc = xmlTreeParse(readLines(sampledocurl),asText = TRUE,useInternalNodes = TRUE);
    
        taxonxml = sampledoc['//TAXON_ID'];
        
        if (length(taxonxml) > 0) {
            result$sampletaxon = xmlToList(taxonxml[[1]]);
        } else {
            log.info(enaexpid, ' Warning: No Taxon ID found in ', sampledocurl);
        }
    
        commonnamexml = sampledoc['//COMMON_NAME'];
    
        if (length(commonnamexml) > 0) {
            result$samplename = xmlToList(commonnamexml[[1]]);
        } else {
            log.info(enaexpid, ' Warning: No Sample Name found in ', sampledocurl);
        }
    
    
    } else {
        log.info(enaexpid, ' Warning: No Sample ID found in ', expdocurl);
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
getAESDRFFilename = function(accession) {
    trace.enter("getAESDRFFilename");
    on.exit({ trace.exit() })
    
    # non trivial logic since there 
    # can be 2 types of SDRF files
    #
    
    filename1 = paste(getAEExperimentPathBase(),
        substr(accession,3,6),accession,paste(accession,'.sdrf.txt', sep=''), sep='/');
    
    if (!file.exists(filename1)) {
        
        # no usual SDRF, check if 
        # .seq.sdrf.txt SDRF is present 
        #
        
        filename2 = paste(getAEExperimentPathBase(),
            substr(accession,3,6),accession,paste(accession,'.seq.sdrf.txt', sep=''), sep='/');
        
        if (!file.exists(filename2)) {
            log.info('SDRF ',filename1, ' does not exist');
            log.info('SDRF ',filename2, ' does not exist');
            
            # return the filename of .sdrf.txt
            #
            return(filename1);
        } else {
            # return the filename of .seq.sdrf.txt 
            #
            return(filename2);
        }
    
    } else {
        # return the filename of .sdrf.txt 
        #
        return(filename1);
    }
    
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
    
    if (!file.exists(filename)) {
        log.info('IDF ',filename, ' does not exist');
    }
    
    return(filename);
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
    
    return(NULL);
}

# read organism from SDRF
#
#
getOrganismFromSDRF = function(sdrf) {
    trace.enter("getOrganismFromSDRF");
    on.exit({ trace.exit() })
    
    if (!is.null(sdrf)) {
        index = grep("Characteristics.Organism", names(sdrf@data));
        if (length(index) > 0) {
            return(sdrf@data[[ sort(names(sdrf@data)[index])[1] ]]);
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
        log.info(accession, " Error: SDRF not found ", sdrffname);
    } else {
        sdrf0 = try(readSDRF(sdrffname),silent = TRUE);
    
        if (inherits(sdrf0, 'try-error')) {
            log.info(accession, " Error: cannot read SDRF ", sdrffname);
        } else {
            #log.info('getting EXP IDs..');
            
            record$assayname = getAssayNameFromSDRF(sdrf0)
            
            record$sdrfrunid = getRunIDFromSDRF(sdrf0);
            
            record$sdrfexpid = getENAExpIDFromSDRF(sdrf0);
            
            record$sdrforganism = getOrganismFromSDRF(sdrf0);
            
            #log.info('getting RUN IDs..');
            if (!is.null(record$sdrfexpid)) {
            
                if(any(record$sdrfexpid == "")) {
                    log.info(accession," Warning: Missing EXPIDs in SDRF");
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
            
                if (!is.null(record$sdrfrunid)) {
                    #
                    #
                    if(any(record$sdrfrunid == "")) {
                        log.info(accession," Warning: Missing RUNIDs in SDRF");
                    }
                    
                } else {
                    log.info(accession," Warning: ENA RUNID section in SDRF not found");
                }
                 
                
                # check fastq files are there
                #
                fastq = unlist(sapply(record$enaexpid, function(x)
                { 
                    #[[1]]$enaexpid$SRX023137$enafastq
                    #[1] "fastq/SRR059/SRR059744/SRR059744_1.fastq.gz"
                    #[2] "fastq/SRR059/SRR059744/SRR059744_2.fastq.gz"
                    
                    if (is.null(x$enafastq)) {
                        log.info(accession,"-", x, "-  Error: No FASTQ reference found");
                    }
                    
                    fastqfnames = paste(getENAFastqPathBase(), x$enafastq, sep="/");
                    
                    sapply(fastqfnames, function(z){
                        if (!file.exists(z)) {
                            log.info(accession," Warning: FASTQ missing ", z);
                        }
                    })
                    
                    x$enafastq; 
                    
                }));
                
                record$allenarunid = sapply(fastq, function(x){ unlist(strsplit(x, '/'))[3] });
                 
                # TODO: check the number of runs 
                # corresponds to the number of fastq
                #
            
            } else {
                log.info(accession," Error: NO references to ENA ");
            }
        }
    }
    
    
    return(record);
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
            log.info('stopping...');
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

#accessionlist = getHTSAccessionsFromAE(T); ## HTS
#accessionlist = getHTSRNASeqAccessionsFromAE(T); ## HTS && RNASeq

#xx = sapply(1:length(accessionlist), function(i){ log.info('i=',i); grepAllENAAccessions(accessionlist[[i]]); })
#xx_NAMES = xx[sapply(xx, function(x){ length(names(x$record)) > 0 })]
#xx_NO_EXP = xx_NAMES[sapply(xx_NAMES, function(x) { !("X.Comment..ENA_EXPERIMENT.." %in% names(x$record)) } )]

#table(sapply(xx, function(x){ length(names(x$record)) }))
#xx_2 = xx[sapply(xx, function(x){ length(names(x$record)) == 2 })]
#xx_1 = xx[sapply(xx, function(x){ length(names(x$record)) == 1 })]
#xx_3 = xx[sapply(xx, function(x){ length(names(x$record)) == 3 })]
#xx_4 = xx[sapply(xx, function(x){ length(names(x$record)) == 4 })]
#xx_5 = xx[sapply(xx, function(x){ length(names(x$record)) == 5 })]


#zz = sapply(1:length(accessionlist), function(i){ log.info('i=',i); mapAEtoENAviaHTTP(accessionlist[[i]]); })
#zz = sapply(1:length(accessionlist), function(i){ log.info(accessionlist[i]," i=",i); mapAEtoENAviaHTTP(accessionlist[i]); })

# PURE Mus musculus
#
#zz_Mus = zz[sapply(zz, function(x){ !is.null(x$enaexpid) && all(sapply(x$enaexpid, function(y){ (!is.null(y$samplename) && y$samplename == "Mus musculus") })) })]
#any(sapply(zz_Mus, function(x){ any(sapply(x$enaexpid, function(y) { (y$layout$name == "PAIRED") } )) }))
#> length(zz_Mus)
#[1] 50
#table(unlist(sapply(zz_Mus, function(x){ sapply(x$enaexpid, function(y){ y$samplename }) })))
#Mus musculus 
#429 


# HTS + RNA Seq
#> length(zz_Mus_Hybrid)
#[1] 16
#> accessionlist
 
#> length(zz_Homo_Hybrid)
#[1] 31

# HYBRID Mus musculus + another organism
#
#zz_Mus_Hybrid = zz[sapply(zz, function(x){ any(sapply(x$enaexpid, function(y){ (!is.null(y$samplename) && y$samplename == "Mus musculus") }))})]
#length(zz_Homo_Hybrid)
#[1] 56

#> match(sapply(zz_Mus_Hybrid, function(x){ x$accession }) , sapply(zz_Mus, function(x){ x$accession }))
#[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 NA 24
#[26] 25 26 27 28 29 30 31 32 33 34 35 NA 36 37 38 39 40 41 42 43 NA 44 45 46 47
#[51] 48 49 50 51 52 53

#table(unlist(sapply(zz_Mus_Hybrid, function(x){ sapply(x$enaexpid, function(y){ y$samplename }) })))
#
#Apis mellifera      Arabidopsis thaliana          Canis familiaris 
#1                         4                        16 
#Chlamydomonas reinhardtii        Ciona intestinalis               Danio rerio 
#1                         1                         2 
#Gallus gallus              Homo sapiens     Monodelphis domestica 
#12                        21                        10 
#Mus musculus              Oryza sativa       Populus trichocarpa 
#468                         1                         1 



# PURE Homo sapiens
#
#zz_Homo = zz[sapply(zz, function(x){ !is.null(x$enaexpid) && all(sapply(x$enaexpid, function(y){ (!is.null(y$samplename) && y$samplename == ""Homo sapiens") })) })]
#> length(zz_Homo)
#[1] 56

# HYBRID Homo sapiens + another organism
#
#sum(sapply(zz, function(x){ any(sapply(x$enaexpid, function(y){ (!is.null(y$samplename) && y$samplename == "Homo sapiens") }))}))
#[1] 60
#zz_Homo_Hybrid = zz[sapply(zz, function(x){ !is.null(x$enaexpid) && any(sapply(x$enaexpid, function(y){ (!is.null(y$samplename) && y$samplename == "Homo sapiens") }))})]
#> length(zz_Homo_Hybrid)
#[1] 60

#> match(sapply(zz_Homo_Hybrid, function(x){ x$accession }) , sapply(zz_Homo, function(x){ x$accession }))
#[1]  1 NA  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
#[26] 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
#[51] 50 51 52 53 54 55 NA 56 57 58

#table(unlist(sapply(zz_Homo_Hybrid, function(x){ sapply(x$enaexpid, function(y){ y$samplename }) })))
#
#Homo sapiens Homo sapiens brain       Mus musculus 
#721                  1                 24 


#zz_Homo = zz[sapply(zz, function(x){ (!is.null(x$enaexpid[[1]]) && x$enaexpid[[1]]$samplename == "Homo sapiens"); })]
#any(sapply(zz_Homo, function(x){ any(sapply(x$enaexpid, function(y) { (y$layout$name == "PAIRED") } )) }))

#sapply(zz_ENA,function(x) { 
    #    log.info('accession = ', x$accession); 
    #    sapply(unlist(x$enaexpid), function(x) { log.info('runid = ',x, ' expid = ',enadatafull$EXPERIMENT_ID[which(enadatafull$RUN_ID == x)]) } )})


#total = length(unlist(sapply(zz, function(x){ x$enaexpid[[1]]$samplename })))
#sum(unlist(sapply(zz, function(x){ x$enaexpid[[1]]$samplename == "Homo sapiens" })))
#[1] 58


#rawlist = sapply(accessionlist, function(x){ a = list(); a$data=list(accession = x, rawdata=getAERawDataFilenames(x)); return(a); })
#zz_Incomplete = zz[sapply(zz, function(x){ (length(unique(x$sdrfrunid)) != length(unique(x$allenarunid))) })]
#sapply(zz, function(x){ sapply(x$enaexpid, function(y){ (y$samplename == "Homo sapiens"); })})

#sapply(accessionlist, function(x){ message(paste(getAEexptype(x), collapse=", ")); })

# SPREADSHEET
#
# LAYOUT
# sapply(zz1, function(x){   message(   unique(sapply(x$enaexpid, function(y) { y$layout$name; } ))     )    })
#
# SAMPLE NAME
# sapply(zz1, function(x){   message(   unique(sapply(x$enaexpid, function(y) { y$samplename; } ))     )    })
#
# SDRF SAMPLE NAME
# sapply(zz1, function(x){   message(   unique(x$sdrforganism)     )    })
#
# ENA RUNS
# sapply(zz1, function(x){   message( length(unique(unlist(sapply(x$enaexpid, function(y) {  sapply(y$enafastq, function(gg){ strsplit(gg, split="/")[[1]][3] })  } )))) )   })
# 
# SDRF RUNS
# sapply(zz1, function(x) { message(  length( unique(x$sdrfrunid[ which(x$sdrfrunid != "") ]) ) ) })
#
# EXP ID
# sapply(zz1, function(x) { message(length(unique(names(x$enaexpid)))) })
#
# SDRF ID
# sapply(zz1, function(x) { message(  length( unique(x$sdrfexpid[ which(x$sdrfexpid != "") ]) ) ) })
#
# CHECK EF
# sapply(zz1, function(x){ y = FALSE; if (!is.null(x$sdrffname) && file.exists(x$sdrffname)) { s0 = readSDRF(x$sdrffname); y = any(grep("Factor", names(s0@data))) }; message(y); y; })
#
#


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

getHTSAccessionsFromAE = function(refresh = FALSE) {
    trace.enter("getHTSAccessionsFromAE");
    on.exit({ trace.exit() })
    
    queryfilename = paste("query", "HTS", ".xml", sep = "")
    
    if (!file.exists(queryfilename) || refresh) {
        
        #query = paste("http://wwwdev.ebi.ac.uk/microarray-as/ae/xml/v2/experiments",
        #"?keywords=&species=&array=&exptype[]=",
        #"&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        #"&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        query = paste("http://wwwdev.ebi.ac.uk/arrayexpress/xml/v2/experiments",
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
        
        #query = paste("http://wwwdev.ebi.ac.uk/microarray-as/ae/xml/v2/experiments",
        #"?keywords=&species=&array=&exptype[]=",
        #"&exptype[]=high+throughput+sequencing+assay&pagesize=25",
        #"&sortby=releasedate&sortorder=descending&expandefo=on",sep='');
        
        #query = paste("http://wwwdev.ebi.ac.uk/arrayexpress/xml/v2/experiments",
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
    
    expurl = paste("http://wwwdev.ebi.ac.uk/arrayexpress/xml/v2/experiments/", accession, sep="");
    
    expdoc = try(xmlTreeParse(readLines(expurl), asText = TRUE, useInternalNodes = TRUE), silent = TRUE);
    
    if (inherits(expdoc, 'try-error')) {
        message("Error!");
        return(NULL);
    } else {
        links = expdoc["//experimenttype"];
        return(sapply(links, function(x){ xmlToList(x)[[1]] } ));
    }
}

getAEExpDescription = function(accession) {
    
    expurl = paste("http://wwwdev.ebi.ac.uk/arrayexpress/xml/v2/experiments/", accession, sep="");
    
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

