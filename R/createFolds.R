createFolds <- function(target, k_fold = 10) {
    target <- as.factor(target)
    if (k_fold < length(target)) {
        target <- factor(as.character(target))
        numInClass <- table(target)
        foldVector <- vector(mode = "integer", length(target))
        for (i in seq_len(length(numInClass))) {
            min_reps <- numInClass[i] %/% k_fold
            if (min_reps > 0) {
                spares <- numInClass[i] %% k_fold
                seqVector <- rep(seq_len(k_fold), min_reps)
                if (spares > 0) {
                    seqVector <- c(seqVector, sample(
                        seq_len(k_fold),
                        spares
                    ))
                }
                foldVector[which(target == names(numInClass)[i])] <-
                    sample(seqVector)
            }
            else {
                foldVector[which(target == names(numInClass)[i])] <-
                    sample(seq_len(k_fold), size = numInClass[i])
            }
        }
    }
    outfoldVector <- split(seq(along = target), foldVector)
    names(outfoldVector) <- paste("Fold", gsub(" ", "0", format(
        seq(along = outfoldVector)
    )), sep = "")
    outfoldVector
}
