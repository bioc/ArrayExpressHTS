#
# .C and .Call wrappers
#

count_polyL <- function(letter, num, seq){
    trace.enter("count_polyL");
    on.exit({ trace.exit() })
    
    .Call("count_polyL", 
        as.character(letter), 
        as.integer(num), 
        as.character(seq), 
        PACKAGE="ArrayExpressHTS")
} 


phred_to_average_qual <- function(num, quals){
    trace.enter("phred_to_average_qual");
    on.exit({ trace.exit() })
    
    .Call("phred_to_average_qual", 
        as.integer(num), 
        as.character(quals), 
        PACKAGE="ArrayExpressHTS")
} 


is_polyX <- function(num, limit, base, seq){
    trace.enter("is_polyX");
    on.exit({ trace.exit() })
    
    .Call("is_polyX", 
        as.integer(num), 
        as.integer(limit), 
        as.character(base), 
        as.character(seq), 
        PACKAGE="ArrayExpressHTS")
} 

