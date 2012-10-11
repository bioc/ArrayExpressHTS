/*
 Author: Andrew Tikhonov <andrew@ebi.ac.uk>
 
*/

#include <stdio.h> 
#include <string.h>  
#include <stddef.h>
#include <stdlib.h>
//#include <R.h>
//#include <Rdefines.h>

void pastefield(char *dest, char **last, char *text1, char *text2){
    int length1 = strlen(text1);
    int length2 = strlen(text2);
    
    strncpy(dest, text1, length1);
    strncpy(dest + length1, text2, length2);
    
    *last = dest + length1 + length2; 
}


#define LINESIZE 1000
#define LINESEP "\t"

void addXS( const char *infname, const char *outfname, int *res ) {
    char line[LINESIZE]; 
    char result[LINESIZE]; 
    
    FILE *infile;
    FILE *outfile; 
    
    char *pattern = "XS:"; 
    
    char *ptr; 
    char *last;
    char *last2; 
    
    int index; 
    
    Rprintf("What?\n");
    if((infile = fopen(infname, "r")) == NULL) { 
        Rprintf("Error Opening File %s\n", infname); 
        res[0]=1; 
    } else if((outfile = fopen(outfname, "w")) == NULL) { 
        Rprintf("Error Creating File %s\n", outfname); 
           res[0]=1;  
    } else {
        while( fgets(line, sizeof(line), infile) != NULL ) { 
            //Rprintf("l=%s", line);
            
            pastefield(result, &last2, strtok_r (line, LINESEP, &last), LINESEP); // first field
            
            ptr = strtok_r (NULL, LINESEP, &last); // second field
            pastefield(last2, &last2, ptr, LINESEP);
            
            index = atoi( ptr );
            
            int i;
            for (i = 2; i < 11; i++) { // skip 9 fields
                pastefield(last2, &last2, strtok_r (NULL, "\t\n", &last), LINESEP);
            }
            
            while( (ptr = strtok_r (NULL, "\t\n", &last)) != NULL ) {
                if (ptr[0] == pattern[0] && ptr[1] == pattern[1] && ptr[2] == pattern[2]) {
                    // skip it
                } else {
                    pastefield(last2, &last2, ptr, LINESEP);
                }
            }
            
            pastefield(last2, &last2, ((index & 0x10) == 0 ? "XS:A:+": "XS:A:-"), "\n");
            
            // terminate
            last2[0] = 0;
            
            //Rprintf("r=%s", result);
            //fprintf(outfile, "%s", result);
            fputs(result, outfile);
        } 
        
        fclose(infile); 
        fclose(outfile); 
        res[0] = 1;
    }
}

