flagToBitStr = function(flag) {
    
    bits01 = intToBits(flag);
    str01 = "";
    
    for (i in 1:11) {
        if (bits01[i] == 1) {
            str01 = paste(str01, "1", sep = "");
        } else {
            str01 = paste(str01, "0", sep = "");
        }
    }
    
    return (str01);
}


parseFlag = function(flag) {
    
    trace.enter("parseFlag");
    on.exit({ trace.exit() })
    
    record01 = list();
    
    # Flag     Chr   Description
    # 0x0001    p    the read is paired in sequencing
    # 0x0002    P    the read is mapped in a proper pair
    # 0x0004    u    the query sequence itself is unmapped
    # 0x0008    U    the mate is unmapped
    # 0x0010    r    strand of the query (1 for reverse)
    # 0x0020    R    strand of the mate
    # 0x0040    1    the read is the first read in a pair
    # 0x0080    2    the read is the second read in a pair
    # 0x0100    s    the alignment is not primary
    # 0x0200    f    the read fails platform/vendor quality checks
    # 0x0400    d    the read is either a PCR or an optical duplicate    
    
    record01$readIsPaired = bitAnd(flag, 0x0001) == 0x0001;
    record01$mappedInPair = bitAnd(flag, 0x0002) == 0x0002;
    record01$readUnmapped = bitAnd(flag, 0x0004) == 0x0001;
    record01$mateUnmapped = bitAnd(flag, 0x0008) == 0x0008;
    record01$readStrand = bitAnd(flag, 0x0010) == 0x0010;
    record01$mateStrand = bitAnd(flag, 0x0020) == 0x0020;
    record01$readIsFirstInPair = bitAnd(flag, 0x0040) == 0x0040;
    record01$readIsSecondInPair = bitAnd(flag, 0x0080) == 0x0080;
    record01$alignmentNotPrimary = bitAnd(flag, 0x0100) == 0x0100;
    record01$readQualityFailure = bitAnd(flag, 0x0200) == 0x0200;
    record01$readIsPCRDuplicate = bitAnd(flag, 0x0400) == 0x0400;
    
    return(record01);
}

printFlag = function(flag) {
    record01 = parseFlag(flag);
    
    message("LEGEND:    pPuUrR12sfd");
    message("FLAGS:     ", flagToBitStr(flag));
    
    if(record01$readIsPaired) {
        message("p = 1, the read is paired in sequencing");
    } else {
        message("p = 0, the read is NOT paired in sequencing");
    }
    
    if(record01$mappedInPair) {
        message("P = 1, the read is mapped in a proper pair");
    } else {
        message("P = 0, the read is NOT mapped in a proper pair");
    }
    
    if(record01$readUnmapped) {
        message("u = 1, the query sequence itself is unmapped");
    } else {
        message("u = 0, the query sequence itself is NOT unmapped (mapped)");
    }
    
    if(record01$mateUnmapped){
        message("U = 1, the mate is unmapped");
    } else {
        message("U = 0, the mate is NOT unmapped (mapped)");
    }
    
    if (record01$readStrand) {
        message("r = 1, strand of the query is reverse");
    } else {
        message("r = 0, strand of the query NOT reverse");
    }
    
    if (record01$mateStrand) {
        message("R = 1, strand of the mate is reverse");
    } else {
        message("R = 0, strand of the mate NOT reverse");
    }
    
    if (record01$readIsFirstInPair) {
        message("1 = 1, the read is the first read in a pair");
    } else {
        message("1 = 0, the read is NOT the first read in a pair");
    }
    
    if (record01$readIsSecondInPair) {
        message("2 = 1, the read is the second read in a pair");
    } else {
        message("2 = 0, the read is NOT the second read in a pair");
    }
    
    if (record01$alignmentNotPrimary) {
        message("s = 1, the alignment is not primary");
    } else {
        message("s = 0, the alignment is primary");
    }
    
    if (record01$readQualityFailure) {
        message("f = 1, the read fails platform/vendor quality checks");
    } else {
        message("f = 0, the read passes platform/vendor quality checks");
    }
    
    if (record01$readIsPCRDuplicate) {
        message("d = 1, the read is either a PCR or an optical duplicate");
    } else {
        message("d = 0, the read is NOT a PCR or an optical duplicate");
    }
}

getSAMLineInfo = function(samline) {
    
    trace.enter("getSAMLineInfo");
    on.exit({ trace.exit() })
    
    fields = unlist(strsplit(samline, "\t"));
    
    info = list();
    info$flags   = parseFlag(fields[2]);
    info$refname = fields[3];
    info$refpos  = fields[4];
    
    return(info);
}

findUsefulRange = function(samtools, bamfile, range) {
    
    trace.enter("findUsefulRange");
    on.exit({ trace.exit() })
    
    lines01 = system(paste(samtools, " view ",bamfile," | head -n 100", sep=""), intern = TRUE)
    
    if (length(lines01) > 0) {
        info = list();
        
        for (line in lines01) {
            info = getSAMLineInfo(line);
            
            # first mapped read
            if (info$flags$readUnmapped == FALSE) {
                return(paste(info$refname, ":", info$refpos, "-", as.integer(info$refpos) + range ,sep=""));
            }
        }
        
        # return info
        return(paste(info$refname, ":", info$refpos, "-", as.integer(info$refpos) + range ,sep=""));
    } else {
        return(NULL);
    }
}

writePropertyFile <- function(organism, version, region, filename) {
    
    trace.enter("writePropertyFile");
    on.exit({ trace.exit() })
    
    # "Homo sapiens"
    # "GRCh37.64
    #"1:69284467-69290714"
    
    write.table(data.frame(list("Organism" = organism, "Version" = version, "Region" = region)), 
            row.names = FALSE, col.names = TRUE, file = filename, sep="\t", quote = FALSE);
    
}

recordAlignedProperties <- function( ) {
    
    trace.enter("recordAlignedProperties");
    on.exit({ trace.exit() })
    
    sortedbam      = paste(.project$aligner$out_dir, "/",  "accepted_hits.sorted.bam", sep="");
    sortedbamprops = paste(sortedbam,  ".prop", sep="");
    
    if(file.exists(sortedbam)) {
        
        samtoolscmd = paste(getPipelineOption("ArrayExpressHTS.samtools"), "/samtools", sep="");
        
        usefulrange = findUsefulRange(samtoolscmd, sortedbam, 10000);
        
        if (is.null(usefulrange)) {
            log.warning("BAM is EMPTY ", sortedbam)
            
            usefulrange = "1:1000-10000"; # hardcoded
        }
        
        writePropertyFile(sub("_", " ", .project$organism), 
                .project$reference$version, usefulrange, sortedbamprops); # "1:1-10000"
        
    }
}

