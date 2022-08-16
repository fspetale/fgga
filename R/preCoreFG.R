.onLoad <- function(libname = NULL, pkgname="fgga"){
    btc <- NULL
    assign('bfc',BiocFileCache() , envir = .GlobalEnv)
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
    codeRelations <- c("is_a", "BFO_0000050", "BFO_0000066", "RO_0002211",
        "RO_0002212", "RO_0002213", "RO_0002215", "RO_0002216")
    nameRelations <- c("is_a", "part_of", "occurs_in", "regulates",
                        "negatively_regulates", "positively_regulates",
                        "capable_of", "capable_of_part_of")
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
                switch(x$val[y], "biological_process" = "GOBP",
                "molecular_function" = "GOMF", "cellular_component" = "GOCC"))
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
            goPlusMatrixEdges[indexGOTerms, 5] <-
                goPlusMatrixEdges[indexGOTerms, 4]
        } else {
        goPlusMatrixEdges[indexGOTerms, 5] <- vapply( goPlusEdges[[i]][
            seq.int(2 * goRelationLength + 1, 3 * goRelationLength)],
            FUN = function(x, y) {y[x, 2]}, y = goTerms,
            FUN.VALUE = character(1))
        }
        goPlusMatrixEdges[indexGOTerms, 6] <- paste(
            goPlusMatrixEdges[indexGOTerms, 4],
            goPlusMatrixEdges[indexGOTerms, 5], sep = "-")
    }

    colnames(goPlusMatrixEdges) <- c("Son", "Relationship", "Father",
                                    "Domain_Son", "Domain_Father", "Domains")
    return(goPlusMatrixEdges)
}

.readOntology <- function(file, nameOnto, relationshipsOntology) {
    minimal <- TRUE
    goTermRegexp <- "^\\[(Term|Typedef|Instance)\\]"
    tagRegexp <- "^(relationship: )?([^ \t]*[^:]):?\\s+(.+)"
    rawLines <- readLines(file)
    m <- regexpr(text = rawLines, pattern = "^([^!{]+[^!{ \t])")
    lines <- regmatches(x = rawLines, m = m)
    termLines <- grep(pattern = goTermRegexp, x = lines)
    if (length(termLines) == 0) {
        stop("No terms detected in the Ontology")}
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
    propagateLines <- which(tags %in% relationshipsOntology)
    parents <- unname(lapply(FUN = unique, split(
        values[propagateLines], cut(taggedLines[propagateLines],
                breaks = c(termLines, Inf), labels = seq(length(termLines))))))
    tagLines <- which(tags %in% useTags)
    properties <- mapply(SIMPLIFY = FALSE, FUN = function(vals, lns) {
        unname(split(vals, cut(lns, breaks = c(termLines, Inf), labels =
                seq(length(termLines)))))}, split(values[tagLines],
                tags[tagLines]), split(taggedLines[tagLines], tags[tagLines]))
    simplify <- intersect(names(properties), c("id", "name", "def", "comment",
                        "is_obsolete", "created_by", "creation_date"))
    properties[simplify] <- lapply(properties[simplify], function(lst) {
        vapply(lst, "[", 1, FUN.VALUE = character(1))})
    names(properties) <- gsub(x = names(properties), pattern =
                "^((parents)|(children)|(ancestors))$", replacement = "\\1_OBO")
    do.call(what = .ontologyIndex, c(list(ontology = nameOnto ,version = substr(
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
    int.chld <- c(split(as.integer( factor(unlist(use.names = FALSE, chld),
                    levels = id)), unlist(use.names = FALSE, mapply(
                    SIMPLIFY = FALSE, FUN = rep, id, vapply(chld, length,
                    FUN.VALUE = 0L)))), setNames(nm = setdiff(id,
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
                        version = NULL, ontology = "GO") {
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
    if (ontology == "GO"){
        warning(paste0("Some parent terms not found: ",paste0(collapse = ", ",
                missing_terms[seq(min(length(missing_terms), 3))]),
                if (length(missing_terms) > 3) {
                    paste0(" (", length(missing_terms) - 3, " more)")
                    } else {""}))}
    parents <- lapply(parents, intersect, id)
    }
    children <- c(lapply(FUN = as.character, X = split(unlist(
        use.names = FALSE, rep(id, times = vapply(parents, FUN = length,
        FUN.VALUE = 0L))), unlist(use.names = FALSE, parents))),
        setNames(nm = setdiff(id, unlist(use.names = FALSE, parents)),
                rep(list(character(0)), length(setdiff(id,
                unlist(use.names = FALSE, parents))))))[id]
        structure(lapply(FUN = setNames, nm = id, X = list(id = id, name = name,
            parents = parents, children = children,
            ancestors = .strAncsFromPars(id, unname(parents), unname(children)),
            obsolete = obsolete)), version = version)
}

.createTableRef <- function(domains){
    if (domains == "PO"){
        ontoPlusRef <- matrix(c("","part_of","", "GOCC", "PO", "GOCC-PO"), 16,
                                6, byrow = TRUE)
    ontoPlusRef[, 1] <-
        c("GO:0090406", "GO:0043667", "GO:0090404", "GO:0043673", "GO:0048226",
        "GO:0043670", "GO:0043671", "GO:0043669", "GO:0043672", "GO:0043676",
        "GO:0043668", "GO:0043078", "GO:0043680", "GO:0043678", "GO:0035619",
        "GO:0031143")
    ontoPlusRef[, 3] <-
        c("PO:0025195", "PO:0025281", "PO:0025195", "PO:0025195", "PO:0025281",
        "PO:0025281", "PO:0025281", "PO:0025281", "PO:0025281", "PO:0025281",
        "PO:0025281", "PO:0009089", "PO:0000078", "PO:0025281", "PO:0000256",
        "PO:0030055")
    ontoPlusRef[c(12, 13, 15, 16), 2] <- c("located_in", "develops_from",
                                            "develops_from", "participates_in")
    }
    if (domains == "ZFA"){
        ontoPlusRef <- matrix(c("","is_a","", "GOCC", "ZFA", "GOCC-ZFA"), 4, 6,
                                byrow = TRUE)
        ontoPlusRef[, 1] <-
            c("GO:0060417", "GO:0044301", "GO:1990032", "GO:0031633")
        ontoPlusRef[, 3] <-
            c("ZFA:0000084", "ZFA:0001711", "ZFA:0001713", "ZFA:0009198")
    }
    if (domains == "HPO"){
        ontoPlusRef <- matrix(c("", "characteristic_of", "", "GOBP", "HPO",
                            "GOBP-HPO"), 457, 6, byrow = TRUE)
        ontoPlusRef[, 1] <-
        c("GO:0006954", "GO:0006954", "GO:0060073", "GO:0006954", "GO:0042703",
        "GO:0042703", "GO:0006954", "GO:0006954", "GO:0006954", "GO:0007605",
        "GO:0007605", "GO:0002526", "GO:0006954", "GO:0002544", "GO:0006954",
        "GO:0055127", "GO:0007608", "GO:0006954", "GO:0007601", "GO:0007601",
        "GO:0006954", "GO:0006954", "GO:0007601", "GO:0006954", "GO:0002118",
        "GO:0019226", "GO:0072350", "GO:0045136", "GO:0030252", "GO:0045136",
        "GO:0061696", "GO:0030252", "GO:0042703", "GO:0003008", "GO:0070977",
        "GO:0043473", "GO:0043473", "GO:0043473", "GO:0006954", "GO:0043473",
        "GO:0006582", "GO:0001503", "GO:0050890", "GO:0060004", "GO:0060004",
        "GO:0006954", "GO:0048513", "GO:0060004", "GO:0001503", "GO:0005739",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0009790", "GO:0007565",
        "GO:0007567", "GO:0060047", "GO:0007565", "GO:0007565", "GO:0006954",
        "GO:0006954", "GO:0007565", "GO:0007565", "GO:0008152", "GO:0046951",
        "GO:0046323", "GO:0030097", "GO:0006954", "GO:0030421", "GO:0050892",
        "GO:0006954", "GO:0032098", "GO:0006954", "GO:0060464", "GO:0007585",
        "GO:0042065", "GO:0022010", "GO:0001964", "GO:0001764", "GO:0050881",
        "GO:0030431", "GO:0050881", "GO:0016049", "GO:0042465", "GO:0042747",
        "GO:0060004", "GO:0006954", "GO:0006954", "GO:0001503", "GO:0001503",
        "GO:0042742", "GO:0071746", "GO:0006955", "GO:0045728", "GO:0070227",
        "GO:0070977", "GO:0006954", "GO:0007585", "GO:0007585", "GO:0045453",
        "GO:0030282", "GO:0006954", "GO:0071754", "GO:0007585", "GO:0008483",
        "GO:0007059", "GO:0008203", "GO:0042592", "GO:0006629", "GO:0022011",
        "GO:0019226", "GO:0014732", "GO:0008203", "GO:0004494", "GO:0071743",
        "GO:0051319", "GO:0002185", "GO:0071736", "GO:0006281", "GO:0019172",
        "GO:0071746", "GO:0070288", "GO:0005739", "GO:0008152", "GO:0001958",
        "GO:0034435", "GO:0004658", "GO:0030424", "GO:0031594", "GO:0005604",
        "GO:0030261", "GO:0071746", "GO:0042552", "GO:0030424", "GO:0071754",
        "GO:0006555", "GO:0004333", "GO:0070527", "GO:0005622", "GO:0008482",
        "GO:0016137", "GO:0014734", "GO:0005739", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0003010", "GO:0042571",
        "GO:0071736", "GO:0001503", "GO:0001503", "GO:0006520", "GO:0009072",
        "GO:0000096", "GO:0009235", "GO:0008152", "GO:0006687", "GO:0001573",
        "GO:0006144", "GO:0006206", "GO:0019752", "GO:0006551", "GO:0006801",
        "GO:0006631", "GO:0055074", "GO:0006568", "GO:0006096", "GO:0001659",
        "GO:0030203", "GO:0006954", "GO:0007608", "GO:0007608", "GO:0030282",
        "GO:0001503", "GO:0006012", "GO:0010960", "GO:0060343", "GO:0002544",
        "GO:0006954", "GO:0001503", "GO:0043473", "GO:0030183", "GO:0006955",
        "GO:0042113", "GO:0042110", "GO:0001503", "GO:0001503", "GO:0071743",
        "GO:0035265", "GO:0043473", "GO:0070977", "GO:0001503", "GO:0035282",
        "GO:0070977", "GO:0001503", "GO:0005929", "GO:0006094", "GO:0001503",
        "GO:0001503", "GO:0002544", "GO:0070166", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0005739", "GO:0043209", "GO:0022410",
        "GO:0032289", "GO:0043209", "GO:0043473", "GO:0043473", "GO:0043473",
        "GO:0043473", "GO:0043473", "GO:0043473", "GO:0043473", "GO:0043473",
        "GO:0043473", "GO:0008283", "GO:0043473", "GO:0042552", "GO:0043473",
        "GO:0060004", "GO:0043473", "GO:0043473", "GO:0043473", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503", "GO:0004565",
        "GO:0046544", "GO:0045136", "GO:0055069", "GO:0005739", "GO:0032286",
        "GO:0005739", "GO:0005739", "GO:0001503", "GO:0001503", "GO:0006776",
        "GO:0035878", "GO:0001503", "GO:0055127", "GO:0046543", "GO:0007283",
        "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0005746", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0043473", "GO:0007605", "GO:0006954", "GO:0042695",
        "GO:0046543", "GO:0006783", "GO:0046903", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0042571", "GO:0042571", "GO:0019230", "GO:0009081",
        "GO:0006558", "GO:0009069", "GO:0006544", "GO:1901052", "GO:0009066",
        "GO:0006566", "GO:0006555", "GO:0009064", "GO:0006541", "GO:0006547",
        "GO:0006560", "GO:0006553", "GO:0006525", "GO:0006549", "GO:0006573",
        "GO:0006570", "GO:0006534", "GO:0050667", "GO:0072507", "GO:0055080",
        "GO:0055067", "GO:0009112", "GO:0046110", "GO:0001676", "GO:0008152",
        "GO:0008152", "GO:0009437", "GO:1903509", "GO:0006664", "GO:0030218",
        "GO:0007598", "GO:0007597", "GO:0072377", "GO:0043648", "GO:0032787",
        "GO:0005976", "GO:0005975", "GO:0042593", "GO:0007049", "GO:0030261",
        "GO:1903510", "GO:0033559", "GO:0006693", "GO:0008015", "GO:0050878",
        "GO:0006000", "GO:0003014", "GO:0060073", "GO:0070293", "GO:0043209",
        "GO:0006954", "GO:0001816", "GO:0032602", "GO:0006954", "GO:0022010",
        "GO:0022011", "GO:0055064", "GO:0016477", "GO:0043473", "GO:0043473",
        "GO:0060047", "GO:0001503", "GO:0001503", "GO:0001503", "GO:0008152",
        "GO:0005753", "GO:0000052", "GO:0061116", "GO:0006591", "GO:0060004",
        "GO:0005739", "GO:0005739", "GO:0005739", "GO:0006600", "GO:0006954",
        "GO:0031514", "GO:0006563", "GO:0043473", "GO:0001503", "GO:0043473",
        "GO:0043473", "GO:0070288", "GO:0070085", "GO:0006486", "GO:0006487",
        "GO:0006493", "GO:0036066", "GO:0004095", "GO:0042552", "GO:0042552",
        "GO:0006954", "GO:0006954", "GO:0007268", "GO:0042696", "GO:0008283",
        "GO:0060073", "GO:0043473", "GO:0006954", "GO:0006954", "GO:0006954",
        "GO:0001503", "GO:0001503", "GO:0006954", "GO:0042713", "GO:0006954",
        "GO:0042756", "GO:0006954", "GO:0022011", "GO:0001503", "GO:0061696",
        "GO:0061696", "GO:0043473", "GO:0006954", "GO:0070227", "GO:0070227",
        "GO:0002250", "GO:0019221", "GO:0019221", "GO:0043473", "GO:0006953",
        "GO:0043473", "GO:0031052", "GO:0005739", "GO:0005739", "GO:0005746",
        "GO:0043473", "GO:0043473", "GO:0046903", "GO:0002185", "GO:0070288",
        "GO:0006954", "GO:0042730", "GO:0046323", "GO:0001503", "GO:0001503",
        "GO:0001503", "GO:0001503", "GO:0001843", "GO:0006954", "GO:0016234",
        "GO:0006955", "GO:0006766", "GO:0019852", "GO:0042359", "GO:0042360",
        "GO:0060073", "GO:0055062", "GO:0006954", "GO:0006954", "GO:0050890",
        "GO:0040007", "GO:0001503", "GO:0006954", "GO:0006954", "GO:0006954",
        "GO:0006954", "GO:0006954", "GO:0006954", "GO:0043473", "GO:0042755",
        "GO:0046541", "GO:0043473", "GO:0042373", "GO:0001503", "GO:0043084",
        "GO:0000819", "GO:0043473", "GO:0043473", "GO:0042110", "GO:0001503",
        "GO:0006954", "GO:0044307")
        ontoPlusRef[, 3] <-
        c("HP:0000024", "HP:0000031", "HP:0000103", "HP:0000123", "HP:0000140",
        "HP:0000141", "HP:0000230", "HP:0000246", "HP:0000265", "HP:0000364",
        "HP:0000365", "HP:0000371", "HP:0000388", "HP:0000389", "HP:0000403",
        "HP:0000405", "HP:0000458", "HP:0000498", "HP:0000504", "HP:0000505",
        "HP:0000509", "HP:0000554", "HP:0000618", "HP:0000704", "HP:0000718",
        "HP:0000762", "HP:0000816", "HP:0000823", "HP:0000824", "HP:0000826",
        "HP:0000837", "HP:0000845", "HP:0000858", "HP:0000871", "HP:0000927",
        "HP:0000953", "HP:0001000", "HP:0001010", "HP:0001101", "HP:0001106",
        "HP:0001107", "HP:0001216", "HP:0001249", "HP:0001265", "HP:0001284",
        "HP:0001287", "HP:0001338", "HP:0001347", "HP:0001362", "HP:0001427",
        "HP:0001507", "HP:0001508", "HP:0001510", "HP:0001511", "HP:0001561",
        "HP:0001622", "HP:0001635", "HP:0001660", "HP:0001663", "HP:0001701",
        "HP:0001733", "HP:0001735", "HP:0001928", "HP:0001939", "HP:0001946",
        "HP:0001952", "HP:0001978", "HP:0001997", "HP:0002019", "HP:0002024",
        "HP:0002037", "HP:0002039", "HP:0002090", "HP:0002101", "HP:0002104",
        "HP:0002171", "HP:0002188", "HP:0002267", "HP:0002269", "HP:0002304",
        "HP:0002360", "HP:0002375", "HP:0002446", "HP:0002487", "HP:0002494",
        "HP:0002522", "HP:0002583", "HP:0002586", "HP:0002663", "HP:0002703",
        "HP:0002718", "HP:0002720", "HP:0002721", "HP:0002723", "HP:0002731",
        "HP:0002750", "HP:0002754", "HP:0002789", "HP:0002791", "HP:0002797",
        "HP:0002832", "HP:0002840", "HP:0002850", "HP:0002883", "HP:0002910",
        "HP:0002916", "HP:0003107", "HP:0003110", "HP:0003119", "HP:0003130",
        "HP:0003134", "HP:0003202", "HP:0003209", "HP:0003210", "HP:0003212",
        "HP:0003214", "HP:0003236", "HP:0003237", "HP:0003254", "HP:0003258",
        "HP:0003261", "HP:0003281", "HP:0003287", "HP:0003287", "HP:0003336",
        "HP:0003349", "HP:0003353", "HP:0003390", "HP:0003398", "HP:0003400",
        "HP:0003451", "HP:0003460", "HP:0003469", "HP:0003477", "HP:0003496",
        "HP:0003524", "HP:0003536", "HP:0003540", "HP:0003575", "HP:0003643",
        "HP:0003649", "HP:0003712", "HP:0003737", "HP:0003892", "HP:0003893",
        "HP:0003894", "HP:0003897", "HP:0003914", "HP:0004020", "HP:0004051",
        "HP:0004052", "HP:0004053", "HP:0004233", "HP:0004246", "HP:0004254",
        "HP:0004257", "HP:0004274", "HP:0004280", "HP:0004305", "HP:0004313",
        "HP:0004315", "HP:0004330", "HP:0004331", "HP:0004337", "HP:0004338",
        "HP:0004339", "HP:0004341", "HP:0004342", "HP:0004343", "HP:0004345",
        "HP:0004352", "HP:0004353", "HP:0004354", "HP:0004357", "HP:0004358",
        "HP:0004359", "HP:0004363", "HP:0004365", "HP:0004366", "HP:0004370",
        "HP:0004371", "HP:0004386", "HP:0004408", "HP:0004409", "HP:0004605",
        "HP:0004606", "HP:0004915", "HP:0004921", "HP:0005089", "HP:0005231",
        "HP:0005263", "HP:0005275", "HP:0005336", "HP:0005357", "HP:0005374",
        "HP:0005384", "HP:0005419", "HP:0005451", "HP:0005474", "HP:0005479",
        "HP:0005484", "HP:0005599", "HP:0005616", "HP:0005623", "HP:0005776",
        "HP:0005832", "HP:0005885", "HP:0005938", "HP:0005959", "HP:0006016",
        "HP:0006257", "HP:0006280", "HP:0006285", "HP:0006454", "HP:0006598",
        "HP:0006607", "HP:0006628", "HP:0006789", "HP:0006808", "HP:0006979",
        "HP:0007266", "HP:0007305", "HP:0007406", "HP:0007440", "HP:0007513",
        "HP:0007542", "HP:0007661", "HP:0007680", "HP:0007703", "HP:0007730",
        "HP:0007832", "HP:0007850", "HP:0007894", "HP:0007922", "HP:0007988",
        "HP:0008000", "HP:0008001", "HP:0008002", "HP:0008034", "HP:0008087",
        "HP:0008103", "HP:0008108", "HP:0008134", "HP:0008142", "HP:0008166",
        "HP:0008185", "HP:0008197", "HP:0008277", "HP:0008306", "HP:0008311",
        "HP:0008316", "HP:0008322", "HP:0008369", "HP:0008371", "HP:0008372",
        "HP:0008383", "HP:0008477", "HP:0008513", "HP:0008647", "HP:0008669",
        "HP:0008747", "HP:0008785", "HP:0008788", "HP:0008797", "HP:0008820",
        "HP:0008828", "HP:0008829", "HP:0008972", "HP:0009105", "HP:0009106",
        "HP:0009107", "HP:0009887", "HP:0009900", "HP:0010280", "HP:0010314",
        "HP:0010465", "HP:0010472", "HP:0010514", "HP:0010656", "HP:0010660",
        "HP:0010675", "HP:0010701", "HP:0010702", "HP:0010831", "HP:0010892",
        "HP:0010893", "HP:0010894", "HP:0010895", "HP:0010898", "HP:0010899",
        "HP:0010900", "HP:0010901", "HP:0010902", "HP:0010903", "HP:0010904",
        "HP:0010907", "HP:0010908", "HP:0010909", "HP:0010912", "HP:0010914",
        "HP:0010917", "HP:0010918", "HP:0010919", "HP:0010927", "HP:0010929",
        "HP:0010930", "HP:0010932", "HP:0010933", "HP:0010964", "HP:0010965",
        "HP:0010966", "HP:0010967", "HP:0010968", "HP:0010969", "HP:0010972",
        "HP:0010988", "HP:0010989", "HP:0010990", "HP:0010995", "HP:0010996",
        "HP:0011012", "HP:0011013", "HP:0011014", "HP:0011018", "HP:0011019",
        "HP:0011020", "HP:0011022", "HP:0011023", "HP:0011028", "HP:0011032",
        "HP:0011033", "HP:0011036", "HP:0011037", "HP:0011038", "HP:0011096",
        "HP:0011110", "HP:0011113", "HP:0011115", "HP:0011123", "HP:0011400",
        "HP:0011401", "HP:0011422", "HP:0011489", "HP:0011509", "HP:0011512",
        "HP:0011675", "HP:0011836", "HP:0011849", "HP:0011863", "HP:0011922",
        "HP:0011925", "HP:0011965", "HP:0012021", "HP:0012025", "HP:0012046",
        "HP:0012087", "HP:0012102", "HP:0012103", "HP:0012113", "HP:0012115",
        "HP:0012261", "HP:0012278", "HP:0012293", "HP:0012306", "HP:0012319",
        "HP:0012320", "HP:0012343", "HP:0012345", "HP:0012346", "HP:0012347",
        "HP:0012358", "HP:0012359", "HP:0012380", "HP:0012447", "HP:0012448",
        "HP:0012486", "HP:0012490", "HP:0012535", "HP:0012569", "HP:0012574",
        "HP:0012590", "HP:0012643", "HP:0012647", "HP:0012648", "HP:0012649",
        "HP:0012791", "HP:0012792", "HP:0012819", "HP:0012875", "HP:0025439",
        "HP:0030082", "HP:0030151", "HP:0030172", "HP:0030290", "HP:0030338",
        "HP:0030339", "HP:0030493", "HP:0030683", "HP:0030886", "HP:0030887",
        "HP:0031404", "HP:0031406", "HP:0031407", "HP:0031605", "HP:0033331",
        "HP:0040007", "HP:0040012", "HP:0040013", "HP:0040014", "HP:0040015",
        "HP:0040030", "HP:0040031", "HP:0040075", "HP:0040081", "HP:0040133",
        "HP:0040165", "HP:0040224", "HP:0040270", "HP:0045001", "HP:0045002",
        "HP:0045003", "HP:0045004", "HP:0045005", "HP:0045073", "HP:0100299",
        "HP:0100326", "HP:0100508", "HP:0100509", "HP:0100511", "HP:0100514",
        "HP:0100519", "HP:0100529", "HP:0100533", "HP:0100537", "HP:0100543",
        "HP:0100555", "HP:0100569", "HP:0100577", "HP:0100584", "HP:0100614",
        "HP:0100646", "HP:0100653", "HP:0100662", "HP:0100669", "HP:0100738",
        "HP:0100755", "HP:0100816", "HP:0100831", "HP:0100856", "HP:0200023",
        "HP:0200024", "HP:0200064", "HP:0200098", "HP:0410035", "HP:0430013",
        "HP:0500006", "HP:0500032")
    ontoPlusRef[c(10, 37, 63, 64, 81, 83, 85, 106, 107, 141, 179, 181, 220,
                    358, 384, 385, 417, 435), 2] <- "characteristic_of_part_of"
    ontoPlusRef[c(15, 33, 34, 49, 66, 91, 258, 361, 418), 2] <- "towards"
    ontoPlusRef[c(25, 58, 112, 142, 405, 407), 2] <- "equivalentTo"
    ontoPlusRef[55, 2] <- "existence_overlaps"
    ontoPlusRef[c(105, 120, 127, 140, 373), 2] <- "capable_of"
    ontoPlusRef[c(139, 249), 2] <- "part_of"
    ontoPlusRef[c(158, 105, 113, 114, 120, 127, 137, 140, 245, 373), 4] <-
        "GOMF"
    ontoPlusRef[c(158, 105, 113, 114, 120, 127, 137, 140, 245, 373), 6] <-
        "GOMF-HPO"
    ontoPlusRef[c(31, 50, 92, 103, 115, 117, 118, 121, 122, 123, 128, 129, 130,
                132, 134, 135, 139, 143, 160, 161, 200, 208, 218, 219, 222,
                249, 251, 252, 268, 282, 283, 335, 351, 356, 357, 358, 361,
                367, 395, 396, 408, 409, 410, 414, 415, 425, 457), 4] <- "GOCC"
    ontoPlusRef[c(31, 50, 92, 103, 115, 117, 118, 121, 122, 123, 128, 129, 130,
                132, 134, 135, 139, 143, 160, 161, 200, 208, 218, 219, 222,
                249, 251, 252, 268, 282, 283, 335, 351, 356, 357, 358, 361,
                367, 395, 396, 408, 409, 410, 414, 415, 425, 457), 6] <-
                "GOCC-HPO"
    }
    colnames(ontoPlusRef) <- c("Son", "Relationship", "Father", "Domain_Son",
                                "Domain_Father", "Domains")
    return(ontoPlusRef)
}

.createTableInference <- function(nameRelations, domains = "GO"){
    tableInference <- data.frame(unlist(lapply(nameRelations,
                        FUN = function(x, y) {rep(x, y)},
                        y = length(nameRelations))),
                        rep(nameRelations, length(nameRelations)), 0)
    colnames(tableInference) <- c("fatherRelation", "grandfatherRelation", "h")
    if (domains == "GO"){
        tableInference[c(seq_len(11), 17, 18, seq.int(24, 26),33, 34, 41, 42,
                        49, 50, 57, seq.int(60, 70), seq.int(73, 80)), 3] <- 1
    }
    if (domains == "PO"){
        tableInference[c(seq_len(12), 16, 17, 21, 22),3] <- 1
        tableInference <- tableInference[-c(1,2,6,7),]
    }
    if (domains == "ZFA"){
        tableInference[c(seq_len(12), 16, 17, 21, 22),3] <- 1
        tableInference <- tableInference[-c(1,2,6,7),]
    }
    if (domains == "HPO"){
        tableInference <- tableInference[-c(1, 2, 3, 7, 9, 10, 11, 15, 17, 18,
                                            19, 23, 49, 50, 51, 55),]
        tableInference[c(seq_len(4), 10, 11, 12, 13, 21, 23, 29, 31, 37, 38, 39,
                                            40, 47, 43), 3] <- 1
    }
    return(tableInference)
}

preCoreFG <- function(ontoTerms, domains = "GO") {
    domains <- switch(domains, "GOCC"="GO", "GOBP"="GO", "GOMF"="GO","GO"="GO",
                    "GOCC-PO"="PO","GO-PO"="PO", "GOCC-ZFA"="ZFA",
                    "GO-ZFA"="ZFA", "GOCC-HPO"="HPO", "GOMF-HPO"="HPO",
                    "GOBP-HPO"="HPO", "GO-HPO"="HPO",
                    stop("INCORRECT DOMAIN; \"", domains, "\", ",
                        "SELECTED! --> INVALID OPTION"))

    ontologyPlusRef <- .goPlusInfo()
    nameRelations <- c("is_a", "part_of", "occurs_in", "regulates",
                    "negatively_regulates", "positively_regulates",
                    "capable_of", "capable_of_part_of", "equivalentTo")
    url <- .getCache("http://purl.obolibrary.org/obo/go.obo", "goBasic")
    if (length(url)==0) stop("Neither internet connection nor cached GO data")
    infoOntologies <- .readOntology(url, "GO", nameRelations)
    tableInference <- .createTableInference(nameRelations, "GO")

    if (domains=="PO"){
        newNameRelations <- c("is_a", "part_of", "participates_in",
                            "develops_from", "located_in")
        url <- .getCache("http://purl.obolibrary.org/obo/po.obo", "poBasic")
        if (length(url)==0)
            stop("Neither internet connection nor cached PO data")
        newOntology <- .readOntology(url, domains, newNameRelations)
    }

    if (domains=="ZFA"){
        newNameRelations <- c("is_a", "part_of", "overlaps", "develops_from",
                            "continuous_with")
        url <- .getCache("http://purl.obolibrary.org/obo/zfa.obo", "zfaBasic")
        if (length(url)==0)
            stop("Neither internet connection nor cached ZFA data")
        newOntology <- .readOntology(url, domains, newNameRelations)
    }

    if (domains=="HPO"){
        newNameRelations <- c("is_a")
        url <- .getCache("http://purl.obolibrary.org/obo/hp.obo", "hpoBasic")
        if (length(url)==0)
            stop("Neither internet connection nor cached HPO data")
        newOntology <- .readOntology(url, domains, newNameRelations)
        newNameRelations <- c("is_a", "capable_of", "part_of", "towards",
                        "characteristic_of", "characteristic_of_part_of",
                        "equivalentTo", "existence_overlaps")}

    if (domains!="GO"){
        addTableInference <- .createTableInference(newNameRelations, domains)
        ontologyPlusRef <- rbind(ontologyPlusRef, .createTableRef(domains))
        for (i in seq_len(length(infoOntologies))) infoOntologies[[i]] <-
            c(infoOntologies[[i]], newOntology[[i]])
        tableInference <- rbind(tableInference, addTableInference)
        rm(newOntology, addTableInference)}

    allTerms <- c()
    for (i in seq_len(length(ontoTerms))){
        if (!is.na(infoOntologies$id[ontoTerms[i]]) &&
            infoOntologies$obsolete[ontoTerms[i]][[1]] == FALSE) allTerms <-
            unique(union(allTerms, infoOntologies$ancestors[ontoTerms[i]][[1]]))
    }

    indexOnto <- c()
    for (i in seq_len(length(allTerms))){
        indexOnto <- union(indexOnto, ontologyPlusRef[
            which(ontologyPlusRef[,1] == allTerms[i]),3])
    }
    for (i in seq_len(length(indexOnto))){
        if (!is.na(infoOntologies$id[indexOnto[i]]) &&
            infoOntologies$obsolete[indexOnto[i]][[1]] == FALSE) allTerms <-
            unique(union(allTerms, infoOntologies$ancestors[indexOnto[i]][[1]]))
    }

    myTerms <- sort(c(allTerms, "Ontology:FES"))
    matrixOntology <- matrix(0, length(myTerms), length(myTerms))
    colnames(matrixOntology) <- rownames(matrixOntology) <- myTerms
    for (i in seq_len(length(myTerms))) {
        ontologyNodes <- infoOntologies$parents[[myTerms[i]]]
        indexTerms <- which(ontologyPlusRef[, 1] == myTerms[i])
        if (length(ontologyNodes) > 0) {
            for (j in seq_len(length(ontologyNodes))) {
                if (ontologyNodes[j] %in% myTerms) {
                    matrixOntology[ontologyNodes[j], i] <- 1
                }
            }
        if (length(indexTerms) > 0) {
            for (j in seq_len(length(indexTerms))) {
                if (ontologyPlusRef[indexTerms[j], 4] !=
                    ontologyPlusRef[indexTerms[j], 5]) {
                        if (length(which(myTerms ==
                            ontologyPlusRef[indexTerms[j], 3])) == 1) {
                            matrixOntology[ontologyPlusRef[indexTerms[j], 3],
                            ontologyPlusRef[indexTerms[j], 1]] <- 1
                            relationships <- ontologyPlusRef[
                            which(ontologyPlusRef[indexTerms[j], 3] ==
                            ontologyPlusRef[, 1]), ]
                        if (is.matrix(relationships)) {
                            for (k in seq_len(dim(relationships)[1])) {
                                hValue <- tableInference[intersect(
                                    which(tableInference[, 2] ==
                                    relationships[k, 2]),
                                    which(tableInference[, 1] ==
                                    ontologyPlusRef[indexTerms[j], 2])), 3]
                                if (hValue == 0) {
                                    matrixOntology[relationships[k, 1],
                                    ontologyPlusRef[indexTerms[j], 3]] <- 0
                                }
                            }
                    } else {
                    hValue <- tableInference[intersect(
                        which(tableInference[, 2] == relationships[2]),
                        which(tableInference[, 1] ==
                            ontologyPlusRef[indexTerms[j], 2])), 3]
                    if (hValue == 0) {
                        matrixOntology[relationships[1],
                                ontologyPlusRef[indexTerms[j], 3]] <- 0
                    }}}}
            }}
        }
    }

    rootGO <- c("GO:0008150", "GO:0003674", "GO:0005575")
    rootGO <- intersect(myTerms, rootGO)

    if (domains == "PO") {
        newRoots <- c()
        if (length(which(allTerms=="PO:0025131")) > 0) newRoots <- "PO:0025131"
        if (length(which(allTerms=="PO:0009012")) > 0)
            newRoots <- c(newRoots, "PO:0009012")
        matrixOntology["Ontology:FES",  c(rootGO, newRoots)] <- 1
    }
    if (domains=="ZFA") {
        newRoots <- c()
        if (length(which(allTerms == "ZFA:0100000")) > 0)
            newRoots <- "ZFA:0100000"
        matrixOntology["Ontology:FES",  c(rootGO, newRoots)] <- 1
    }
    if (domains=="HPO") {
        newRoots <- c()
        if (length(which(allTerms == "HP:0000001")) > 0)
            newRoots <- "HP:0000001"
        matrixOntology["Ontology:FES",  c(rootGO, newRoots)] <- 1
    }
    if (domains == "GO") {
        matrixOntology["Ontology:FES", rootGO] <- 1
    }

    return(as(matrixOntology, "graphNEL"))
}
