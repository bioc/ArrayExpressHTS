#include <stdio.h>  
#include "samtools/sam.h"
#include "hashtable.h"
#include "hashtable_itr.h"
#include "hashtable_utility.h"

typedef struct {  
    int beg, end;  
    samfile_t *in;  
} tmpstruct_t;  

struct key {
    unsigned char str[50];
    int length;
};

struct value {
    int mm;
    unsigned int index[80]; // maximum 40 alignments 
    int length;
};

static unsigned int hashfromkey(void *ky) {
    struct key *k = (struct key *)ky;
    unsigned int hash = 5381;
    int c, i;
    
    for( i = 0; i < k->length; i++ ) {
        c = k->str[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    
    return hash;
}

static int equalkeys(void *ky1, void *ky2) {
    struct key *k1 = (struct key *)ky1;
    struct key *k2 = (struct key *)ky2;
    //Rprintf( "Comparing |%s| with |%s| results in %d\n", k1->str, k2->str, strcmp(k1->str, k2->str) );
    return(strcmp(k1->str, k2->str) == 0);
}

void print_key(void *ky) {
    struct key *k = (struct key *)ky;
    Rprintf("ID: %s, %d\n", k->str, k->length);
}

void print_value( void *va ) {
    struct value *v = (struct value *)va;
    int i;
    Rprintf( "\tIndexes for %d mismatches of length %d: ", v->mm, v->length );
    for( i = 0; i < v->length; i++ ) {
        Rprintf( "%d ", v->index[i] );
    }
    Rprintf( "\n" );
}

int int_cmp(const void *a, const void *b) 
{ 
    const int *ia = (const int *)a; // casting pointer types 
    const int *ib = (const int *)b;
    return *ia  - *ib; 
    /* integer comparison: returns negative if b > a 
     and positive if a > b */ 
} 

void fltbam(int* argc, char *argv[], int* result) {  
    /* the real processor */

    Rprintf( "Welcome, creating structures...\n" );

    Rprintf( "argc = %d\n", argc[0] );
    int dd0 = 0;
    for (dd0 = 0;dd0 < argc[0];dd0++) {
        Rprintf( "argv[%d] = %s\n", dd0, argv[dd0] );
    }    

    tmpstruct_t tmp;  
    tmpstruct_t tmp_out;
    Rprintf( "Checking args...\n" );
    if (argc[0] != 3) {  
        Rprintf("Usage: fltbam <in.bam> <out.bam>\n");  
        result[0] = 2;
        return;  
    }  
    
    Rprintf( "Args OK, reading input file...\n" );
    tmp.beg = 0; tmp.end = 0x7fffffff;  
        tmp.in = samopen(argv[1], "rb", 0);  
        if (tmp.in == 0) {  
        Rf_error("Fail to open BAM file %s\n", argv[1]);  
        result[0] = 3;
        return;
    }
    
    struct key *k, *kk;
    struct value *found, *v;
    struct hashtable *h;
    //struct hashtable_itr *itr;
    
    Rprintf( "Creating hashtable...\n" );
    h = create_hashtable(1000, hashfromkey, equalkeys);
    
    bam1_t *line = bam_init1();
    int count = 0, i = 0;
    int tindex = 0;
    int mate = 0;

    Rprintf( "Populating hashtable...\n" );
    while( 0 < samread(tmp.in, line) ) {
        //Rprintf( "Type %d. File pointer %x\n", tmp.in->type&1, *(tmp.in->x.bam) );
        //Rprintf("flag %d\n", line->core.flag);
        if( i % 1000000 == 0 ) 
        Rprintf( "read %d\n", i );
        
        if( line->core.flag & BAM_FPROPER_PAIR ) {
            if( mate ) {    
                k = (struct key *)malloc(sizeof(struct key));
                if( NULL == k ) {
                    Rprintf("Ran out of memory allocating a key\n");
                    result[0] = 4;
                    return;
                }
                
                strcpy(k->str, bam1_qname(line));
                k->length = strlen(k->str);
                
                //print_key( k ); 
                //Rprintf("\tAuxdata is: %d\n", bam_aux2i(bam_aux_get(line,"NM")) );    
                found = hashtable_search(h,k);
                count = count + bam_aux2i(bam_aux_get(line,"NM"));
                if( NULL == found ) {
                    //Rprintf("\tKey not found, inserting.\n");
                    found = (struct value *)malloc(sizeof(struct value));
                    found->mm = count;
                    found->index[0] = i-1;
                    found->length = 1;
                    tindex = tindex + 1;
                    if( !hashtable_insert(h,k,found) ) { result[0] = 5; return; }
                    //print_value(found);
                    //Rprintf( "\tIndexes: %d\n", tindex );
                } else {
                    free(k);
                    if( count < found->mm ) {
                        //Rprintf("\tFound key less than count, replacing.\n");
                        found->mm = count;
                        found->index[0] = i-1;
                        tindex = tindex-(found->length) + 1;
                        found->length = 1;
                        //hashtable_change(h, k, found); 
                        //print_value(found);
                        //Rprintf( "\tIndexes: %d\n", tindex );
                    } else if( count == found->mm ) {
                        //Rprintf("\tFound key same as count, appending.\n");
                        found->index[found->length] = i-1;
                        tindex = tindex + 1;
                        found->length = found->length+1;
                        //print_value(found);
                        //Rprintf( "\tIndexes: %d\n", tindex );
                        //} else {
                        //    Rprintf("\tBigger count, doing nothing.\n");
                        //    Rprintf( "\tIndexes: %d\n", tindex );
                    }
                }
                
            } else {
                count = bam_aux2i(bam_aux_get(line,"NM"));
            }
        } else {
            count = 0;
        }
        i++;
        mate = !mate;
    }
    Rprintf("Hashtable contains %u keys, and %d indexes.\n", hashtable_count(h), tindex);
    Rprintf( "Closing input file...\n" );
    samclose(tmp.in); 
    
    Rprintf( "Creating new index array...\n" );    
    //int allind[tindex];    
    int *allind = malloc(tindex * sizeof (int));
    Rprintf( "Creating iterator...\n" );
    struct hashtable_itr *itr = hashtable_iterator(h);
    int j;
    i = 0;
    Rprintf( "Copying indexes...\n" );
    if (hashtable_count(h) > 0) {
        do {
            kk = hashtable_iterator_key(itr);
            v = hashtable_iterator_value(itr);
            for( j = 0; j < v->length; j++ )  
                allind[i+j] = v->index[j];
            i = i + v->length;
        } while (hashtable_iterator_advance(itr));
    }
    Rprintf( "%d indexes copied...\n" , tindex );
    //for( i = 0; i < tindex; i++ )
    //    Rprintf( "%d ", allind[i] );
    //Rprintf( "\n" );
    Rprintf( "Sorting...\n" );
    qsort (allind, tindex, sizeof(int), int_cmp);
    //Rprintf( "%d indexes copied:" , sizeof(allind)/sizeof(int) );
    //for( i = 0; i < tindex; i++ )
    //    Rprintf( "%d ", allind[i] );
    //Rprintf( "\n" );
    Rprintf( "Destroying hashtable...\n" );
    hashtable_destroy(h, 1);
    Rprintf( "Destroying iterator...\n" );
    free(itr);
    
    Rprintf( "Opening input file again...\n" );
    tmp.beg = 0; tmp.end = 0x7fffffff;
    tmp.in = samopen(argv[1], "rb", 0);
    if (tmp.in == 0) {
        Rf_error("Fail to open BAM file %s\n", argv[1]);
        result[0] = 6;
        return;
    }
    i = 0;
    j = 0;
    Rprintf( "Opening output file %s...\n", argv[2] );
    tmp_out.beg = 0; tmp_out.end = 0x7fffffff;
    tmp_out.in = samopen(argv[2], "wb", tmp.in->header );
    if (tmp_out.in == 0) {
        Rf_error("Fail to open BAM file %s\n", argv[2]);
        result[0] = 7;
        return;
    }
    Rprintf( "Outputing indexes...\n" );
    while( 0 < samread(tmp.in, line) ) {
        if( j >= tindex ) {
            Rprintf( "Found all indexes, bye!\n" );
            break;
        }
        if( allind[j] == i ) {
            //Rprintf( "\tFound an index %d, writting\n", allind[j] );
            samwrite(tmp_out.in,line);
            samread(tmp.in, line);
            samwrite(tmp_out.in,line);
            i++;
            j++;
        }    
        i++;    
    }
    bam_destroy1(line);
    samclose(tmp.in);
    samclose(tmp_out.in);
    free(allind);
    //free(k);
    //free(found);
    result[0] = 0;
    return;
}

int main(int argc, char *argv[]) {
    /* wrapper */
    int result[1] = { 0 };
    fltbam(&argc, argv, result);
    return (result[0]);
}


