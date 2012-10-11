#include <stdio.h> 
#include <string.h>  
#include <stddef.h> 
#include <R.h>
#include <Rdefines.h>


/*
    Andrew Tikhonov <andrew@ebi.ac.uk>
*/

#define LINESIZE 1000
#define NEWLINE '\n'


void getReadLength(const char** infname , int* length) {  
    char line[LINESIZE]; 
    
    FILE *infile;
    
    /* init to "undefined" */
    length[0] = 0;
    
    if((infile = fopen(infname[0], "r")) == NULL) { 
        Rprintf("Error opening file %s\n", infname[0]); 
        length[0] = -1;
        
    } else {
        fgets(line, sizeof(line), infile); // read first line and discard it
        fgets(line, sizeof(line), infile); // read second line and check its length
        
        fclose(infile); 
        
        for (int i = 0; i < LINESIZE;i++) {
            if (line[i] == NEWLINE || line[i] == 0) {
                length[0] = i;
                break;
            }
        }
    }
}

#define UNDEF   -1
#define PHRED33 33
#define PHRED64 64
#define SOLEXA  59

void checkQuality(const char** infname , int* readmax, int* qualityscore) {  
    char line[LINESIZE]; 
    
    FILE *infile;
    
    /* init to "undefined" */
    qualityscore[0] = 0;
    
    int stop = 0; 
    
    if((infile = fopen(infname[0], "r")) == NULL) { 
        Rprintf("Error opening file %s\n", infname[0]); 
        qualityscore[0] = UNDEF;
        
    } else {
        qualityscore[0] = PHRED64;
    
        for(int i = 0;(i < readmax[0] && !stop); i++) {
            /* Rprintf("i=%d\n", i); */
            
            fgets(line, sizeof(line), infile); // read first line and discard it
            fgets(line, sizeof(line), infile); // read second line
            fgets(line, sizeof(line), infile); // read third line
            fgets(line, sizeof(line), infile); // and forth line
            
            for (int j = 0; j < LINESIZE; j++) {
                
                if (line[j] == NEWLINE || line[j] == 0) {
                    break;
                }
                
                if (line[j] < SOLEXA) {
                    qualityscore[0] = PHRED33;
                    
                    /* Rprintf("detected %d at read %d\n", line[j], i); */
                    
                    stop = 1;
                    break;
                }
            }
        }
        
        fclose(infile); 
    }
    
}



//void count_polyL( char **letter, int *len, char **seq, int *res ) {
SEXP count_polyL(SEXP letter, SEXP len, SEXP seq) {
    SEXP result;
    int *p_result; 
    
    //Rprintf("Printed from countL2:\n");
    
    int length = INTEGER_VALUE(len);
    //Rprintf("length= %d\n", length);
    
    char char0 = CHAR(STRING_ELT(letter, 0))[0];
    //Rprintf("char0= %c\n", char0);
    
    PROTECT(result = NEW_INTEGER(length));
    p_result = INTEGER_POINTER(result);

    //for (int i0=0; i0 < length; i0++) {
    //    Rprintf("i0=%d\n", i0);
    //    char *p_seq = CHAR(STRING_ELT(seq, i0));
    //    int p_seq_strlen = strlen(p_seq);
    //    Rprintf("p_seq_strlen=%d\n", p_seq_strlen);
    //    Rprintf("p_seq[%d]=%s\n", i0, p_seq);
    //}
    //Rprintf("p_result[0]= %d\n", p_result[0]);
    
    int j, count;
 
    for(int i = 0; i < length; i++) { 
        //Rprintf("i= %d\n", i);
    
        count = 0;
        j = 0;
        
        char *p_seq = (char *) CHAR(STRING_ELT(seq, i));
        
        int strlength = strlen(p_seq);
        
        //Rprintf("strlength= %d\n", strlength);
        //Rprintf("p_seq[j]= %c\n", p_seq[j]);
        
        while( j < strlength && p_seq[j] == char0 ) {
            count++;
            j++;
        
            //Rprintf("count= %d\n", count);
            //Rprintf("j= %d\n", j);
        }
        
        p_result[i] = count;
        
        //Rprintf("p_result[i]= %d\n", p_result[i]);
        
        count = 0;
        j = strlength -1;
        
        //Rprintf("j= %d\n", j);
        //Rprintf("p_seq[j]= %c\n", p_seq[j]);

        while( j >= 0 && p_seq[j] == char0 ) {
            count++;
            j--;
        
            //Rprintf("count= %d\n", count);
            //Rprintf("j= %d\n", j);
        }
        
        if( p_result[i] < count ) {
            p_result[i] = count; 
        }
    
        //Rprintf("p_result[i]= %d\n", p_result[i]);
    }
    
    UNPROTECT(1);
    return result;
}


//void phred_to_average_qual( int *len, char **quals, double *res ) {
SEXP phred_to_average_qual(SEXP len, SEXP quals) {
    
    SEXP result;
    double *p_result; 
    
    //Rprintf("Printed from phred_to_average_qual2:\n");
    
    int length = INTEGER_VALUE(len);
    //Rprintf("length= %d\n", length);
    
    PROTECT(result = NEW_NUMERIC(length));
    p_result = NUMERIC_POINTER(result);

    int i, j;
    double r;
    
    for( i = 0; i < length; i++ ) { 
        
        char *p_qual = (char *) CHAR(STRING_ELT(quals, i));
        int strlength = strlen(p_qual);
        
        r = 0;
    
        for( j = 0; j < strlength; j++ ) {
            r = r + p_qual[j] - 33;
        }
        
        p_result[i] = r / strlength; 
    }
    
    UNPROTECT(1);
    return result;
    
}

// int *len, int *flag
SEXP is_firstmate( SEXP len, SEXP flag ) {
    SEXP result;
    int *p_result;

    int length = INTEGER_VALUE(len);

    PROTECT(result = NEW_INTEGER(length));
    p_result = INTEGER_POINTER(result);

    int i;

    for( i = 0; i < length; i++ ) {
        p_result[i] = (INTEGER(flag)[i] & 0x41) == 0x41;
    }

    UNPROTECT(1);
    return result;
}

SEXP is_secondmate( SEXP len, SEXP flag ) {
    SEXP result;
    int *p_result;

    int length = INTEGER_VALUE(len);

    PROTECT(result = NEW_INTEGER(length));
    p_result = INTEGER_POINTER(result);

    int i;

    for( i = 0; i < length; i++ ) {
        p_result[i] = (INTEGER(flag)[i] & 0x81) == 0x81;
    }

    UNPROTECT(1);
    return result;
}

//void is_polyX( int *len, int *limit, char *base, char **seq, int *res ) {
SEXP is_polyX( SEXP len, SEXP lim, SEXP base, SEXP seq) {

    SEXP result;
    int *p_result; 
    
    int length = INTEGER_VALUE(len);
    //Rprintf("length= %d\n", length);
    
    int limit = INTEGER_VALUE(lim);
    //Rprintf("limit= %d\n", limit);
    
    PROTECT(result = NEW_INTEGER(length));
    p_result = INTEGER_POINTER(result);
    
    int i, j, count;
    
    for( i = 0; i < length; i++ ) { 
        count = 0;
        j = 0;
        
        char *p_seq = (char *) CHAR(STRING_ELT(seq, i));
        int strlength = strlen(p_seq);
        
        while( count < limit && j < strlength && p_seq[j] == 'A' ) {
            count++;
            j++;
        }
        if( count >= limit ) {
            p_result[i] = 1;
        } else {
            count = 0;
            j = strlength -1;
            while( count < limit && j >= 0 && p_seq[j] == 'A' ) {
                count++;
                j--;
            }
            if( count >= limit ) {
                p_result[i] = 1;
            }
            else {
                p_result[i] = 0;
            }
        }
    }
    
    UNPROTECT(1);
    return result;

}

