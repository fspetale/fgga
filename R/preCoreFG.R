.onLoad <- function(libname = NULL, pkgname="fgga"){
    btc <- NULL
    assign('bfc',BiocFileCache() , envir = .GlobalEnv )
}

.getCache <- function(url, nameGO){
    verbose <- has_internet()
    rid <- bfcquery(bfc, nameGO, "rname")$rid
    if (!length(rid)) {
        if(verbose) {
            message("Downloading ", nameGO)
            rid <- names(bfcadd(bfc, nameGO, url))}}

    if (!isFALSE(bfcneedsupdate(bfc, rid))) bfcdownload(bfc, rid)

    return (bfcrpath(bfc, rids = rid))
}

.goPlusList <- function(relationship, goPlus) {
    goRelation <- lapply(goPlus, FUN = function(x) {
        (x[intersect(intersect(grep( relationship, x$pred), grep("GO_", x$sub)),
                    grep("GO_", x$obj)), ])})
    indexList <- which(lapply(goRelation,
                        FUN = function(x) (length(x$obj))) > 0)
    goRelation <- do.call("rbind", goRelation[indexList])
    goRelation <- apply(goRelation, MARGIN = 2, FUN = function(x) {
        (gsub("http://purl.obolibrary.org/obo/GO_", "GO:", x))})
    return(goRelation)
}

.goPlusInfo <- function() {
    message("Read  GO-plus")
    url <- .getCache(
        "http://purl.obolibrary.org/obo/go/extensions/go-plus.json", "goPlus")
    if (length(url)==0) stop("Neither internet connection nor cached GO data")
    goTermInfo <- fromJSON(url, simplifyVector = TRUE)
    codeRelations <- c(
        "is_a", "BFO_0000050", "BFO_0000066", "RO_0002211",
        "RO_0002212", "RO_0002213", "RO_0002215", "RO_0002216")
    nameRelations <- c(
        "is_a", "part_of", "occurs_in", "regulates",
        "negatively_regulates", "positively_regulates",
        "capable_of", "capable_of_part_of")
    #   goPlus <- sapply(codeRelations,
    goPlus <- lapply(codeRelations,
                    FUN = .goPlusList, goPlus = goTermInfo$graphs$edges)
    indexGOTerms <- seq(2, length(nameRelations))
    goPlusEdges <- lapply(seq_len(length(nameRelations)),
                            FUN = function(x, y, z) {x[[z]][, 2] <- y[z]
                            return(x[[z]])}, x = goPlus, y = nameRelations)
    goPlusMatrixEdges <- do.call("rbind", goPlusEdges)
    goPlusMatrixEdges <- cbind(goPlusMatrixEdges,
                                matrix(0, dim(goPlusMatrixEdges)[1], 3))

    goTerms <- vapply(goTermInfo$graphs$nodes, FUN = function(x) {
        (gsub("http://purl.obolibrary.org/obo/GO_", "GO:", x$id))},
        FUN.VALUE = character(dim(goTermInfo$graphs$nodes[[1]])[1]))
    goDomains <- vapply(goTermInfo$graphs$nodes[[1]]$meta$basicPropertyValues,
            FUN = function(x) {y <- grep(
            "biological_process|molecular_function|cellular_component", x$val)
            if (length(y) > 0) {return(
                switch(x$val[y], "biological_process" = "BP",
                    "molecular_function" = "MF", "cellular_component" = "CC"))
                    } else {return("")}}, FUN.VALUE = character(1))

    goTerms <- cbind(goTerms, goDomains)
    goTerms <- goTerms[-which(goDomains == ""), ]
    rownames(goTerms) <- goTerms[, 1]
    goTerms <- goTerms[intersect(union(
        goPlusMatrixEdges[, 1], goPlusMatrixEdges[, 3]), rownames(goTerms)), ]
    indexGOTerms <- 0

    for (i in seq_len(length(codeRelations))) {
        goRelationLength <- length(goPlusEdges[[i]]) / 3
        indexGOTerms <- seq.int(max(indexGOTerms) + 1,
                                goRelationLength + max(indexGOTerms))
        goPlusMatrixEdges[indexGOTerms, 4] <- vapply(goPlusEdges[[i]][
            seq_len(goRelationLength)], FUN = function(x, y) {y[x, 2]},
            y = goTerms, FUN.VALUE = character(1))
        if (i == 1) {
            goPlusMatrixEdges[indexGOTerms, 5] <- goPlusMatrixEdges[
                indexGOTerms, 4]
        } else {
            goPlusMatrixEdges[indexGOTerms, 5] <- vapply(
                goPlusEdges[[i]][
                    seq.int(2 * goRelationLength + 1, 3 * goRelationLength)],
                FUN = function(x, y) {y[x, 2]}, y = goTerms,
                FUN.VALUE = character(1))
        }
        goPlusMatrixEdges[indexGOTerms, 6] <- paste(
            goPlusMatrixEdges[indexGOTerms, 4],
            goPlusMatrixEdges[indexGOTerms, 5], sep = "-")
    }

    colnames(goPlusMatrixEdges) <- c(
    "Son", "Relationship", "Father", "Domain_Son", "Domain_Father", "Domains")
    return(goPlusMatrixEdges)
}

.readGO <- function(file) {
    minimal <- TRUE
    goTermRegexp <- "^\\[(Term|Typedef|Instance)\\]"
    tagRegexp <- "^(relationship: )?([^ \t]*[^:]):?\\s+(.+)"
    propagate_relationships <- c("is_a", "part_of", "regulates",
                "negatively_regulates", "positively_regulates", "occurs_in")
    rawLines <- readLines(file)
    m <- regexpr(text = rawLines, pattern = "^([^!{]+[^!{ \t])")
    lines <- regmatches(x = rawLines, m = m)
    termLines <- grep(pattern = goTermRegexp, x = lines)
    if (length(termLines) == 0) {
        stop("No GO-terms detected in Gene Ontology")}
    taggedLines <- grep(pattern = tagRegexp, x = lines)
    tagMatches <- regmatches(
        x = lines[taggedLines], regexec(text = lines[taggedLines],
                                        pattern = tagRegexp))
    tags <- vapply(tagMatches, "[", 3, FUN.VALUE = character(1))
    values <- vapply(tagMatches, "[", 4, FUN.VALUE = character(1))
    allTagTypes <- unique(tags)
    useTags <- if (minimal) {
        intersect(c("id", "name", "is_obsolete"), allTagTypes)
    } else {allTagTypes}
    propagateLines <- which(tags %in% propagate_relationships)
    parents <- unname(lapply(FUN = unique, split(
        values[propagateLines], cut(taggedLines[propagateLines],
            breaks = c(termLines, Inf), labels = seq(length(termLines))))))
    tagLines <- which(tags %in% useTags)
    properties <- mapply(SIMPLIFY = FALSE, FUN = function(vals, lns) {
        unname(split(vals, cut(lns, breaks = c(termLines, Inf), labels =
                seq(length(termLines)))))}, split(values[tagLines],
                tags[tagLines]), split(taggedLines[tagLines], tags[tagLines]))
    simplify <- intersect(names(properties), c("id", "name",
                "def", "comment", "is_obsolete", "created_by", "creation_date"))
    properties[simplify] <- lapply(properties[simplify], function(lst) {
        vapply(lst, "[", 1, FUN.VALUE = character(1))})
    names(properties) <- gsub(x = names(properties), pattern =
                "^((parents)|(children)|(ancestors))$", replacement = "\\1_OBO")
    do.call(what = .ontologyIndex, c(list(version = substr(
        lines[seq(termLines[1] - 1)], 1, 1000), parents = parents,
        id = properties[["id"]], name = properties[["name"]],
        obsolete = if ("is_obsolete" %in% names(properties)) {
            (!is.na(properties[["is_obsolete"]])) &
                properties[["is_obsolete"]] == "true"
        } else {rep(FALSE, length(properties[["id"]]))}), properties[
            -which(names(properties) %in% c("id", "name", "is_obsolete"))]))
}

.strAncsFromPars <- function(id, pars, chld) {
    stopifnot(all(vapply(list(pars, chld), function(x) {
        is.null(names(x)) | identical(names(x), id)}, FUN.VALUE = logical(1))))
    int.pars <- c(split(as.integer(
        factor(unlist(use.names = FALSE, pars), levels = id)), unlist(
            use.names = FALSE, mapply(SIMPLIFY = FALSE, FUN = rep, id, vapply(
            pars, length, FUN.VALUE = 0L)))), setNames(nm = setdiff(id, unlist(
                    use.names = FALSE, pars)), rep(list(integer(0)), length(
                        setdiff(id, unlist(use.names = FALSE, pars))))))[id]
    int.chld <- c(split(as.integer(
        factor(unlist(use.names = FALSE, chld), levels = id)), unlist(
            use.names = FALSE, mapply(SIMPLIFY = FALSE, FUN = rep, id,
            vapply(chld, length, FUN.VALUE = 0L)))), setNames(nm = setdiff(id,
            unlist(use.names = FALSE, chld)), rep(list(integer(0)),
            length(setdiff(id, unlist(use.names = FALSE, chld))))))[id]
            setNames(nm = id, lapply(.ancsFromPars(int.pars, int.chld),
                            function(x) id[x]))
}

.ancsFromPars <- function(pars, chld) {
    ancs <- as.list(seq(length(pars)))
    done <- vapply(pars, function(x) length(x) == 0, FUN.VALUE = logical(1))
    cands <- which(done)
    new.done <- seq_len(length(cands))
    while (!all(done)) {
        cands <- unique(unlist(use.names = FALSE, chld[cands[new.done]]))
        v <- vapply(pars[cands], function(x) all(done[x]),
                    FUN.VALUE = logical(1))
        if (!is.logical(v)) {
            stop("Can't get ancestors for items ",
                paste0(collapse = ", ", which(!done)))}
        new.done <- which(v)
        done[cands[new.done]] <- TRUE
        ancs[cands[new.done]] <- mapply(
        SIMPLIFY = FALSE, FUN = c, lapply(cands[new.done], function(x) {
        unique(unlist(use.names = FALSE, ancs[pars[[x]]]))}), cands[new.done])
    }
    ancs
}

.ontologyIndex <- function(parents, id = names(parents), name = id,
                        obsolete = setNames(nm = id, rep(FALSE, length(id))),
                        version = NULL) {
    if (is.null(id)) {
        stop("Must give non-NULL term IDs: either as 'id' argument or as the
            names of the 'parents' argument")
    }
    if (!is.character(id)) {
        stop("'id' argument must be of class 'character'")
    }
    if (!((is.null(names(parents)) & length(parents) == length(id)) |
            identical(names(parents), id))) {
        stop("`parents` argument must have names attribute identical to `id`
        argument or be the same length")
    }
    missing_terms <- setdiff(unlist(use.names = FALSE, parents), id)
    if (length(missing_terms) > 0) {
        warning(paste0("Some parent terms not found: ",paste0(collapse = ", ",
            missing_terms[seq(min(length(missing_terms), 3))]),
            if (length(missing_terms) > 3) {
                paste0(" (", length(missing_terms) - 3, " more)")
            } else {""}))
        parents <- lapply(parents, intersect, id)
    }
    children <- c(lapply(FUN = as.character, X = split(unlist(
        use.names = FALSE, rep(id,
                    times = vapply(parents, FUN = length, FUN.VALUE = 0L))),
        unlist(use.names = FALSE, parents))), setNames(nm = setdiff(id,
                unlist(use.names = FALSE, parents)), rep(list(character(0)),                                                                                                             length(setdiff(id, unlist(use.names = FALSE, parents))))))[id]
    structure(lapply(FUN = setNames, nm = id, X = list(id = id, name = name,
            parents = parents, children = children,
            ancestors = .strAncsFromPars(id, unname(parents), unname(children)),
            obsolete = obsolete)), version = version)
}

preCoreFG <- function(myGOnames, domains = "goPlus") {
    goPlusRef <- .goPlusInfo()
    nameRelations <- c("is_a", "part_of", "occurs_in", "regulates",
                        "negatively_regulates", "positively_regulates",
                        "capable_of", "capable_of_part_of")
    tableInference <- data.frame(unlist(lapply(nameRelations,
                        FUN = function(x, y) {rep(x, y)},
                        y = length(nameRelations))),
                        rep(nameRelations, length(nameRelations)), 0)
    colnames(tableInference) <- c("fatherRelation", "grandfatherRelation", "h")
    url <- .getCache("http://purl.obolibrary.org/obo/go.obo", "goBasic")
    if (length(url)==0) stop("Neither internet connection nor cached GO data")
    ontologyGO <- .readGO(url)
    tableInference[c(seq_len(11), 17, 18, seq.int(24, 26),33, 34, 41, 42, 49,
                    50, 57, seq.int(60, 62)), 3] <- 1
    if (domains == "BP-MF") crossOntology <- c("BP-MF", "MF-BP")
    if (domains == "BP-CC") crossOntology <- c("BP-CC", "CC-BP")
    if (domains == "MF-CC") crossOntology <- c("MF-CC", "CC-MF")
    if (domains == "goPlus") {
        indexGOs <- unique(unlist(
            lapply(myGOnames, FUN = function(x, y) which(y[, 1] == x),
                    y = goPlusRef)))
    } else {
        gominPlusRef <- rbind(
            goPlusRef[which(goPlusRef[, 6] == crossOntology[1]), ],
            goPlusRef[which(goPlusRef[, 6] == crossOntology[2]), ])
        indexGOs <- unique(unlist(lapply(myGOnames,
                FUN = function(x, y) which(y[, 1] == x), y = gominPlusRef[, ])))
    }
    myGOnames <- sort(c(myGOnames, "GO:FES"))
    matrixGO <- matrix(0, length(myGOnames), length(myGOnames))
    colnames(matrixGO) <- rownames(matrixGO) <- myGOnames
    for (i in seq_len(length(myGOnames))) {
        nodesGO <- ontologyGO$parents[[myGOnames[i]]]
        indexGOs <- which(goPlusRef[, 1] == myGOnames[i])
        if (length(nodesGO) > 0) {
            for (j in seq_len(length(nodesGO))) {
                if (nodesGO[j] %in% myGOnames) {
                    matrixGO[nodesGO[j], i] <- 1
                }
                if (length(indexGOs) > 0) {
                    for (j in seq_len(length(indexGOs))) {
                        if (goPlusRef[indexGOs[j], 4] !=
                            goPlusRef[indexGOs[j], 5]) {
                            if (length(which(myGOnames ==
                                            goPlusRef[indexGOs[j], 3])) == 1) {
                                matrixGO[goPlusRef[indexGOs[j], 3],
                                            goPlusRef[indexGOs[j], 1]] <- 1
                                relationGOs <- goPlusRef[
                                    which(goPlusRef[indexGOs[j], 3] ==
                                            goPlusRef[, 1]), ]
                                if (is.matrix(relationGOs)) {
                                    for (k in seq_len(dim(relationGOs)[1])) {
                                        hValue <- tableInference[intersect(
                                            which(tableInference[, 2] ==
                                                        relationGOs[k, 2]),
                                            which(tableInference[, 1] ==
                                                goPlusRef[indexGOs[j], 2])), 3]
                                        if (hValue == 0) {
                                            matrixGO[relationGOs[k, 1],
                                                goPlusRef[indexGOs[j], 3]] <- 0
                                        }
                                    }
                                } else {
                                    hValue <- tableInference[intersect(
                                        which(tableInference[, 2] ==
                                                    relationGOs[2]),
                                        which(tableInference[, 1] ==
                                                goPlusRef[indexGOs[j], 2])), 3]
                                    if (hValue == 0) {
                                        matrixGO[relationGOs[1],
                                                goPlusRef[indexGOs[j], 3]] <- 0
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (domains == "BP-MF")
        matrixGO["GO:FES", c("GO:0008150", "GO:0003674")] <- 1
    if (domains == "BP-CC")
        matrixGO["GO:FES", c("GO:0008150", "GO:0005575")] <- 1
    if (domains == "MF-CC")
        matrixGO["GO:FES", c("GO:0003674", "GO:0005575")] <- 1
    if (domains == "goPlus") {
        matrixGO["GO:FES", c("GO:0008150", "GO:0003674", "GO:0005575")] <- 1
    }

    return(as(matrixGO, "graphNEL"))
}
