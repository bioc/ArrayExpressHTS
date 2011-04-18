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
    //printf( "Comparing |%s| with |%s| results in %d\n", k1->str, k2->str, strcmp(k1->str, k2->str) );
    return(strcmp(k1->str, k2->str) == 0);
}

void print_key(void *ky) {
    struct key *k = (struct key *)ky;
    printf("ID: %s, %d\n", k->str, k->length);
}

void print_value( void *va ) {
    struct value *v = (struct value *)va;
    int i;
    printf( "\tIndexes for %d mismatches of length %d: ", v->mm, v->length );
    for( i = 0; i < v->length; i++ ) {
        printf( "%d ", v->index[i] );
    }
    printf( "\n" );
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

    fprintf( stdout, "Welcome, creating structures...\n" );

    fprintf( stdout, "argc = %d\n", argc[0] );
    int dd0 = 0;
    for (dd0 = 0;dd0 < argc[0];dd0++) {
        fprintf( stdout, "argv[%d] = %s\n", dd0, argv[dd0] );
    }    

    tmpstruct_t tmp;  
    tmpstruct_t tmp_out;
    fprintf( stdout, "Checking args...\n" );
    if (argc[0] != 3) {  
        fprintf(stderr, "Usage: fltbam <in.bam> <out.bam>\n");  
        result[0] = 2;
        return;  
    }  
    
    fprintf( stdout, "Args OK, reading input file...\n" );
    tmp.beg = 0; tmp.end = 0x7fffffff;  
        tmp.in = samopen(argv[1], "rb", 0);  
        if (tmp.in == 0) {  
        fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);  
        result[0] = 3;
        return;
    }
    
    struct key *k, *kk;
    struct value *found, *v;
    struct hashtable *h;
    //struct hashtable_itr *itr;
    
    fprintf( stdout, "Creating hashtable...\n" );
    h = create_hashtable(1000, hashfromkey, equalkeys);
    
    bam1_t *line = bam_init1();
    int count = 0, i = 0;
    int tindex = 0;
    int mate = 0;

    printf( "Populating hashtable...\n" );
    while( 0 < samread(tmp.in, line) ) {
        //printf( "Type %d. File pointer %x\n", tmp.in->type&1, *(tmp.in->x.bam) );
        //printf("flag %d\n", line->core.flag);
        if( i % 1000000 == 0 ) 
            printf( "read %d\n", i );
        
        if( line->core.flag & BAM_FPROPER_PAIR ) {
            if( mate ) {    
                k = (struct key *)malloc(sizeof(struct key));
                if( NULL == k ) {
                    printf("Ran out of memory allocating a key\n");
                    result[0] = 4;
                    return;
                }
                
                strcpy(k->str, bam1_qname(line));
                k->length = strlen(k->str);
                
                //print_key( k ); 
                //printf("\tAuxdata is: %d\n", bam_aux2i(bam_aux_get(line,"NM")) );    
                found = hashtable_search(h,k);
                count = count + bam_aux2i(bam_aux_get(line,"NM"));
                if( NULL == found ) {
                    //printf("\tKey not found, inserting.\n");
                    found = (struct value *)malloc(sizeof(struct value));
                    found->mm = count;
                    found->index[0] = i-1;
                    found->length = 1;
                    tindex = tindex + 1;
                    if( !hashtable_insert(h,k,found) ) { result[0] = 5; return; }
                    //print_value(found);
                    //printf( "\tIndexes: %d\n", tindex );
                } else {
                    free(k);
                    if( count < found->mm ) {
                        //printf("\tFound key less than count, replacing.\n");
                        found->mm = count;
                        found->index[0] = i-1;
                        tindex = tindex-(found->length) + 1;
                        found->length = 1;
                        //hashtable_change(h, k, found); 
                        //print_value(found);
                        //printf( "\tIndexes: %d\n", tindex );
                    } else if( count == found->mm ) {
                        //printf("\tFound key same as count, appending.\n");
                        found->index[found->length] = i-1;
                        tindex = tindex + 1;
                        found->length = found->length+1;
                        //print_value(found);
                        //printf( "\tIndexes: %d\n", tindex );
                        //} else {
                        //    printf("\tBigger count, doing nothing.\n");
                        //    printf( "\tIndexes: %d\n", tindex );
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
    printf("Hashtable contains %u keys, and %d indexes.\n", hashtable_count(h), tindex);
    printf( "Closing input file...\n" );
    samclose(tmp.in); 
    
    printf( "Creating new index array...\n" );    
    //int allind[tindex];    
    int *allind = malloc(tindex * sizeof (int));
    printf( "Creating iterator...\n" );
    struct hashtable_itr *itr = hashtable_iterator(h);
    int j;
    i = 0;
    printf( "Copying indexes...\n" );
    if (hashtable_count(h) > 0) {
        do {
            kk = hashtable_iterator_key(itr);
            v = hashtable_iterator_value(itr);
            for( j = 0; j < v->length; j++ )  
                allind[i+j] = v->index[j];
            i = i + v->length;
        } while (hashtable_iterator_advance(itr));
    }
    printf( "%d indexes copied...\n" , tindex );
    //for( i = 0; i < tindex; i++ )
    //    printf( "%d ", allind[i] );
    //printf( "\n" );
    printf( "Sorting...\n" );
    qsort (allind, tindex, sizeof(int), int_cmp);
    //printf( "%d indexes copied:" , sizeof(allind)/sizeof(int) );
    //for( i = 0; i < tindex; i++ )
    //    printf( "%d ", allind[i] );
    //printf( "\n" );
    printf( "Destroying hashtable...\n" );
    hashtable_destroy(h, 1);
    printf( "Destroying iterator...\n" );
    free(itr);
    
    printf( "Opening input file again...\n" );
    tmp.beg = 0; tmp.end = 0x7fffffff;
    tmp.in = samopen(argv[1], "rb", 0);
    if (tmp.in == 0) {
        fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
        result[0] = 6;
        return;
    }
    i = 0;
    j = 0;
    printf( "Opening output file %s...\n", argv[2] );
    tmp_out.beg = 0; tmp_out.end = 0x7fffffff;
    tmp_out.in = samopen(argv[2], "wb", tmp.in->header );
    if (tmp_out.in == 0) {
        fprintf(stderr, "Fail to open BAM file %s\n", argv[2]);
        result[0] = 7;
        return;
    }
    printf( "Outputing indexes...\n" );
    while( 0 < samread(tmp.in, line) ) {
        if( j >= tindex ) {
            printf( "Found all indexes, bye!\n" );
            break;
        }
        if( allind[j] == i ) {
            //printf( "\tFound an index %d, writting\n", allind[j] );
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


