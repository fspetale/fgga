library(fgga)
## create data for testing
test_check_graph <- function(){
    if(requireNamespace("graph", quietly=TRUE)){
        TRUE
    }else{
        FALSE
    }
}

test_tpg <- function(){
    x <- tableTPG(3)
    checkEquals(length(x), 20);
}

