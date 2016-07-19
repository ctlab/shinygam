# :ToDo: import it from r-utils
write.tsv <- function(table, dir, file=NULL, gzip=FALSE, row.names=NA, col.names=NA, ...) {
    name <- deparse(substitute(table))
    table <- as.data.frame(table) 
    
    if (is.null(file)) {
        file <- file.path(dir, paste0(name, ".tsv", if (gzip) ".gz"))        
    }

    if (is.na(row.names)) {
        row.names <- is.character(attr(table, "row.names"))
    }

    if (!row.names && is.na(col.names)) {
        col.names=T
    }
    
    for (c in colnames(table)) {
        if (is.character(table[[c]])) {
            table[[c]] <- sub("#", "", table[[c]])            
        }
    }
    
    if (gzip) {
        file <- gzfile(file, "w")
    }
    write.table(table, file, quote=F,
                row.names=row.names, col.names=col.names, sep="\t")
    if (gzip) {
        close(file)
    }
}

removeNAColumns <- function(d) {
    keep <- !sapply(d, all %o% is.na)
    d[, keep, with=F]
}

writeGmwcsInstance <- function(graph.dir, network,
                   nodes.group.by=NULL, 
                   edges.group.by=NULL,
                   group.only.positive=F) {
    
    
    
    dir.create(graph.dir, showWarnings = FALSE)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    synonyms.file <- file.path(graph.dir, "synonyms.txt")
    
    synonyms <- c()
    
    nt <- get.vertex.attributes(network)
    if (!is.null(nodes.group.by)) {
        ntx <- if (group.only.positive) { 
            synonyms <- c(synonyms, with(nt[nt$score <= 0,], name))
            nt[nt$score > 0,] 
        } else { 
            nt 
        }
        if (nrow(ntx) > 0) {
            synonyms <- c(synonyms, aggregate(as.formula(sprintf("name ~ %s", nodes.group.by)),
                                              data=ntx, paste0, collapse=" ")$name)
        }
    } else {
        synonyms <- c(synonyms, nt$name)
    }
    nt <- rename(nt[, c("name", "score")], c("name"="#name"))
    
    et <- get.edge.attributes(network, include.ends = T)
    if (!is.null(edges.group.by)) {
        etx <- if (group.only.positive) { 
            synonyms <- c(synonyms, with(et[et$score <= 0,], sprintf("%s -- %s", from, to)))
            et[et$score > 0,] 
        } else { 
            et 
        }
        if (nrow(etx) > 0) {
            synonyms <- c(synonyms, 
                          aggregate(name ~ edges.group.by,
                                    data=list(
                                        name=sprintf("%s -- %s", etx$from, etx$to),
                                        edges.group.by=etx[[edges.group.by]]),
                                    paste0, collapse=" ")$name
            )    
        }
        
        
    }        
    et <- rename(et[, c("from", "to", "score")], c("from"="#from"))
    
    write.tsv(nt, file=nodes.file)
    write.tsv(et, file=edges.file)
    writeLines(sprintf("%s", synonyms), con=synonyms.file)
    
    if (length(synonyms) == 0) {
        synonyms.file <- NULL
    }
    
    list(nodes.file=nodes.file, 
         edges.file=edges.file, 
         synonyms.file=synonyms.file)
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


sgmwcs.solver <- function(gmwcs, nthreads = 1, timeLimit = -1, nodes.group.by=NULL, edges.group.by=NULL, c.size=NULL, group.only.positive=F, quite=T) {
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (!score.edges) {
            E(network)$score <- 0
        }
        
        graph.dir <- tempfile("graph")
        
        instance <- writeGmwcsInstance(graph.dir = graph.dir,
                                       network=network,
                                       nodes.group.by = nodes.group.by,
                                       edges.group.by = edges.group.by,
                                       group.only.positive = group.only.positive
        )
        
        system2(gmwcs, c("--nodes", instance$nodes.file,
                                "--edges", instance$edges.file,
                                if (!is.null(instance$synonyms.file)) c("--synonyms", instance$synonyms.file, "-B", 1) else NULL,
                                "--threads", nthreads, 
                                "--timelimit", timeLimit,
                                if (!is.null(c.size)) c("-c", c.size) else NULL,
                                "--break"
        ))
        solution.file <- paste0(instance$nodes.file, ".out")
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- GAM:::readGraph(node.file = solution.file,
                               edge.file = paste0(instance$edges.file, ".out"),
                               network = network)
        #attr(res, "optimal") <- any(grepl("SOLVED TO OPTIMALITY", out))
        return(res)
    }
}

sgmwcs.solver.eps <- function(..., eps=-1e-2) {
    solve <- sgmwcs.solver(..., c.size=50)
    solve2 <- sgmwcs.solver(..., c.size=10)
    function(network) {
        m <- solve(network)
        V(m)[score == 0]$score <- eps
        E(m)[score == 0]$score <- eps
        res <- solve2(m)
        return(res)
    }
}


read.table.smart <- function(path, ...) {
    fields <- list(...)    
    conn <- file(path)
    header <- readLines(conn, n=1)
    close(conn)
    
    sep <- "\t"
    for (s in c("\t", " ", ",")) {
        if (grepl(s, header)) {
            sep <- s
            break
        } 
    }

    res <- read.table(path, sep=sep, header=T, stringsAsFactors=F, check.names=F, quote='"')
    res <- as.data.table(res, keep.rownames=is.character(attr(res, "row.names")))
    
    oldnames <- character(0)
    newnames <- character(0)

    normalizeName <- function(x) {
        gsub("[^a-z0-9]", "", tolower(x)) 
    }
    
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
        
    setnames(res, oldnames, newnames)
    res
}

read.table.smart.de <- function(path, ID=ID) {
    read.table.smart(path, ID=ID, pval=c("pvalue"), log2FC=c("log2foldchange", "logfc"), name.orig=c("name"))
}


read.table.smart.de.gene <- function(path) {
    read.table.smart.de(path, ID=c("gene", "entrez", "", "rn", "symbol"))
}

read.table.smart.de.met <- function(path) {
    read.table.smart.de(path, ID=c("metabolite", "kegg", "hmdb", "", "rn"))
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
