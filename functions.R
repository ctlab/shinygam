removeNAColumns <- function(d) {
    keep <- !sapply(d, all %o% is.na)
    d[, keep, with=F]
}

gmwcs.solver <- function (gmwcs, nthreads = 1, timeLimit = -1) {
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        system2(gmwcs, c("-n", nodes.file, "-e", edges.file, 
                         "-m", nthreads, "-t", timeLimit
                         ,             "-b"
        ))
        solution.file <- paste0(nodes.file, ".out")
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- GAM:::readGraph(node.file = solution.file,
                               edge.file = paste0(edges.file, ".out"),
                               network = network)
        return(res)
    }
}

heinz21.solver <- function(heinz2, nthreads = 1, timeLimit = -1) 
{
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (score.edges) {
            network <- GAM:::MWCSize(network.orig)
        }
        BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
        BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
        solution.file <- file.path(graph.dir, "sol.txt")
        system2(paste0(heinz2), c("-n", nodes.file, "-e", edges.file, 
                                  "-o", solution.file, "-m", nthreads, "-v", 
                                  0, "-t", timeLimit))
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- BioNet::readHeinzGraph(node.file = solution.file, network = network, 
                                      format = "igraph")
        if (score.edges) {
            res <- GAM:::deMWCSize(res, network.orig)
        }
        return(res)
    }
}

normalizeName <- function(x) {
    gsub("[^a-z0-9]", "", tolower(x)) 
}


read.table.smart <- function(path, ...) {
    fields <- list(...)    
    conn <- file(path)
    header <- readLines(conn, n=1)
    close(conn)
    
    seps <- c("\t", " ", ",", ";")
    sep <- seps[which.max(table(unlist(strsplit(header, "")))[seps])]
    
    res <- read.table(path, sep=sep, header=T, stringsAsFactors=F, check.names=F, quote='"')
    res <- as.data.table(res, keep.rownames=is.character(attr(res, "row.names")))
    
    oldnames <- character(0)
    newnames <- character(0)
    
    for (field in names(fields)) {        
        if (field %in% colnames(res)) {
            next
        }
        
        z <- na.omit(
            match(
                normalizeName(c(field, fields[[field]])),
                normalizeName(colnames(res))))
        if (length(z) == 0) {
            next
        }
        
        oldnames <- c(oldnames, colnames(res)[z[1]])
        newnames <- c(newnames, field)
    }
    
    logdebug("smart renaming")
    logdebug("from: %s", paste0(oldnames, collapse=" "))
    logdebug("to: %s", paste0(newnames, collapse=" "))
    setnames(res, oldnames, newnames)
    res
}

read.table.smart.de <- function(path, ID=ID) {
    read.table.smart(path, ID=ID, pval=c("pvalue"), log2FC=c("log2foldchange", "logfc"), name.orig=c("name"))
}


read.table.smart.de.gene <- function(path, idsList) {
    res <- read.table.smart.de(path, ID="ID")
    idColumn <- findIdColumn(res, idsList)
    if (idColumn$matchRatio < 0.1) {
        z <- na.omit(
            match(
                normalizeName(c("ID", "gene id", "gene", "entrez", "", "rn", "symbol")),
                normalizeName(colnames(res))))
        if (length(z) == 0) {
            setnames(res, colnames(res)[1], "ID")    
        } else {
            setnames(res, colnames(res)[z[1]], "ID")    
        }
        
    } else {
        if (idColumn$column != "ID" && "ID" %in% colnames(res)) {
            setnames(res, "ID", "ID.old")
        }
        setnames(res, idColumn$column, "ID")
    }
    res
}

read.table.smart.de.met <- function(path) {
    res <- read.table.smart.de(path, ID=c("metabolite", "kegg", "hmdb", "", "rn"))
    if (!"ID" %in% colnames(res)) {
        setnames(res, colnames(res)[1], "ID")
    }
    res
}

renderGraph <- function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)
    
    function() {
        val <- func()
        if (is.null(val)) {
            return(list(nodes=list(), links=list()));
        }
        for (a in list.vertex.attributes(val)) {
            if (!is.numeric(get.vertex.attribute(val, a))) {
                next
            }
            print(a)
            vs <- get.vertex.attribute(val, a)
            vs[which(vs == Inf)] <- 1e100
            vs[which(vs == -Inf)] <- -1e100
            val <- set.vertex.attribute(val, a, index=V(val), value=vs)
        }
        for (a in list.edge.attributes(val)) {
            if (!is.numeric(get.edge.attribute(val, a))) {
                next
            }
            print(a)
            vs <- get.edge.attribute(val, a)
            vs[which(vs == Inf)] <- 1e100
            vs[which(vs == -Inf)] <- -1e100
            val <- set.edge.attribute(val, a, index=E(val), value=vs)
        }
        module2list(val)
    }
}

necessary.de.fields <- c("ID", "pval")

vector2html <- function(v) {
    paste0("<ul>\n",
           paste("<li>", names(v), ": ", v, "</li>\n", collapse=""),
           "</ul>\n")
}

renderJs <- function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)
    
    function() {
        val <- func()
        paste0(val, ";", 
               paste(sample(1:20, 10, replace=T), collapse=""))
    }
}

toJsLiteral <- function(x) {
    if (is(x, "numeric")) {
        if (x == Inf) {
            return("Infinity")
        } else {
            return("-Infinity")
        }
        return(as.character(x))
    } else if (is(x, "character")) {
        return(shQuote(x));
    } else if (is(x, "logical")) {
        return(if (x) "true" else "false")
    } else {
        stop(paste0("can't convert ", x, " to JS literal"))
    }
}

makeJsAssignments  <- function(...) {
    args <- list(...)
    values <- sapply(args, toJsLiteral)
    paste0(names(values), " = ", values, ";\n", collapse="")
}

# adapted from shiny
simpleSelectInput <- function (inputId, choices, selected = NULL) 
{
    selectTag <- tags$select(id = inputId)
    optionTags <- mapply(choices, names(choices), SIMPLIFY = FALSE, 
                         USE.NAMES = FALSE, FUN = function(choice, name) {
                             optionTag <- tags$option(value = choice, name)
                             if (choice %in% selected) 
                                 optionTag$attribs$selected = "selected"
                             optionTag
                         })
    selectTag <- tagSetChildren(selectTag, list = optionTags)
    selectTag
}

generateFDRs <- function(es) {
    res <- ""
    num.positive = 150
    if (is.null(es$fb.met) != is.null(es$fb.rxn)) {
        num.positive <- num.positive / 2
    }
    
    if (!is.null(es$fb.met)) {
        fb <- es$fb.met
        pvals <- with(es$met.de.ext, { x <- pval; names(x) <- ID; na.omit(x) })            
        recMetFDR <- GAM:::recommendedFDR(fb, pvals, num.positive=num.positive)
        recAbsentMetScore <- min(GAM:::scoreValue(fb, 1, recMetFDR), -0.1)
        if (!is.null(es$fb.rxn)) {
            recAbsentMetScore <- recAbsentMetScore * 5
        }
        res <- paste0(res, sprintf('$("#metLogFDR").val(%.1f).trigger("change");', log10(recMetFDR)))
        res <- paste0(res, sprintf('$("#absentMetScore").val(%.1f).trigger("change");', recAbsentMetScore))
    }
    
    if (!is.null(es$fb.rxn)) {
        fb <- es$fb.rxn
        pvals <- with(es$rxn.de.ext, { x <- pval; names(x) <- ID; na.omit(x) })            
        recRxnFDR <- GAM:::recommendedFDR(fb, pvals, num.positive=num.positive)
        res <- paste0(res, sprintf('$("#geneLogFDR").val(%.1f).trigger("change");', log10(recRxnFDR)))
    }
    
    res
}


.intersectionSize <- function(...) { length(intersect(...))}

findIdColumn <- function(de, idsList,
                         sample.size=1000,
                         match.threshold=0.6,
                         remove.ensembl.revisions=TRUE) {
    # first looking for column with base IDs
    de.sample <- if (nrow(de) < sample.size) {
        copy(de)
    } else {
        de[sample(seq_len(nrow(de)), sample.size), ]
    }
    columnSamples <- lapply(de.sample, as.character)
    
    
    if (remove.ensembl.revisions) {
        columnSamples <- lapply(columnSamples, gsub,
                                pattern="(ENS\\w*\\d*)\\.\\d*",
                                replacement="\\1")
    }
    
    ss <- sapply(columnSamples,
                 .intersectionSize, idsList[[1]])
    
    if (max(ss) / nrow(de.sample) >= match.threshold) {
        # we found a good column with base IDs
        return(list(column=colnames(de)[which.max(ss)],
                    type=names(idsList)[1],
                    matchRatio=max(ss) / nrow(de.sample)))
    }
    
    z <- .pairwiseCompare(.intersectionSize,
                          columnSamples,
                          idsList)
    
    bestMatch <- which(z == max(z), arr.ind = TRUE)[1,]
    return(list(column=colnames(de)[bestMatch["row"]],
                type=names(idsList)[bestMatch["col"]],
                matchRatio=max(z) / nrow(de.sample)))
}

.pairwiseCompare <- function (FUN, list1, list2 = list1, ...)
{
    additionalArguments <- list(...)
    f1 <- function(...) {
        mapply(FUN = function(x, y) {
            do.call(FUN, c(list(list1[[x]], list2[[y]]), additionalArguments))
        }, ...)
    }
    z <- outer(seq_along(list1), seq_along(list2), FUN = f1)
    rownames(z) <- names(list1)
    colnames(z) <- names(list2)
    z
}