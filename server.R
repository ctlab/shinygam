library(shiny)
library(data.table)
library(igraph)
library(GAM)


options(shiny.error=traceback)
load("./data/kegg.human.network.rda")
load("./data/kegg.mouse.network.rda")

networks <- list(
    "Mouse musculus"=kegg.mouse.network,
    "Homo sapiens"=kegg.human.network)

heinz.py <- "/usr/local/lib/heinz/heinz.py"
mwcs.path <- "/usr/local/bin/mwcs"

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

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {
    
    longProcessStart <- function() {
        session$sendCustomMessage(type='showWaitMessage', list(value=T))
    }

    longProcessStop <- function() {
        session$sendCustomMessage(type='showWaitMessage', list(value=F))
    }

    geneDEInput <- reactive({
        if (is.null(input$geneDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        res <- data.table(read.table(input$geneDE$datapath, sep="\t", header=T, stringsAsFactors=F))
        if (!all(necessary.de.fields %in% names(res))) {
            stop(paste0("Genomic differential expression data should contain at least these fields: ", 
                        paste(necessary.de.fields, collapse=", ")))
        }
        attr(res, "name") <- input$geneDE$name
        res
    })
    
    geneIdsType <- reactive({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        network <- networks[[isolate(input$network)]]
        gene.id.map <- network$gene.id.map
        res <- getIdType(data$ID, gene.id.map)
        if (length(res) != 1) {
            stop("Can't determine type of IDs for genes")
        }
        res
    })
    
    output$geneDESummary <- renderUI({
        gene.de <- geneDEInput()
        ids.type <- geneIdsType()
        if (is.null(gene.de)) {
            return("There is no genomic data")
        }
        
        div(
            HTML(
                vector2html(c(
                    "name" = attr(gene.de, "name"),
                    "length" = nrow(gene.de),
                    "ID type" = ids.type
                ))),
            p("Top DE genes:"))
    })
    
    output$geneDETable <- renderTable({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(as.data.frame(head(data[order(pval)])), digits=3)
    })
    
    
    metDEInput <- reactive({
        if (is.null(input$metDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        res <- data.table(read.table(input$metDE$datapath, sep="\t", header=T, stringsAsFactors=F))
        if (!all(necessary.de.fields %in% names(res))) {
            stop(paste0("Metabolic differential expression data should contain at least these fields: ", 
                        paste(necessary.de.fields, collapse=", ")))
        }
        attr(res, "name") <- input$metDE$name
        res
    })
    
    metIdsType <- reactive({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        GAM:::lazyData("met.id.map")
        res <- getIdType(data$ID, met.id.map)
        if (length(res) != 1) {
            stop("Can't determine type of IDs for metabolites")
        }
        res
    })
    
    output$metDESummary <- renderUI({
        met.de <- metDEInput()
        ids.type <- metIdsType()
        if (is.null(met.de)) {
            return("There is no metabolic data")
        }
        
        div(
            HTML(
                vector2html(c(
                    "name" = attr(met.de, "name"),
                    "length" = nrow(met.de),
                    "ID type" = ids.type
                ))),
            p("Top DE metabolites:"))
    })
    
    output$metDETable <- renderTable({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(as.data.frame(head(data[order(pval)])), digits=3)
    })
    
    esInput <- reactive({
        input$preprocess
        network <- networks[[isolate(input$network)]]
        gene.de <- isolate(geneDEInput())
        gene.ids <- isolate(geneIdsType())

        met.de <- isolate(metDEInput())
        met.ids <- isolate(metIdsType())

        if (is.null(gene.de) && is.null(met.de)) {
            return(NULL)
        }

        longProcessStart()

        if (!is.null(gene.de)) {
            gene.de <- gene.de[which(gene.de$pval < 1),]
        }

        if (!is.null(met.de)) {
            met.de <- met.de[which(met.de$pval < 1),]
        }
        
        reactions.as.edges = isolate(input$reactionsAs) == "edges"
        collapse.reactions = isolate(input$collapseReactions)
        use.rpairs = isolate(input$useRpairs)
        
        es <- makeExperimentSet(
            network=network,
            met.de=met.de, gene.de=gene.de,
            met.ids=met.ids, gene.ids=gene.ids,
            reactions.as.edges=reactions.as.edges,
            collapse.reactions=collapse.reactions,
            use.rpairs=use.rpairs,
            plot=F)
        longProcessStop()
        es
    })
    
    output$networkSummary <- reactive({
        es <- esInput()
        net <- es$subnet
        if (is.null(net)) {
            return("There is no built network")
        }
        
        vector2html(c(
            "number of nodes" = length(V(net)),
            "number of edges" = length(E(net))
            ))
    })
    
    output$networkParameters <- reactive({
        es <- NULL
        tryCatch({
            es <- esInput()
        }, error=function(e) {})
        
        paste0(
            makeJsAssignments(
                network.available = !is.null(es),
                network.hasReactionsAsNodes = !is.null(es) && !es$reactions.as.edges,
                network.hasReactionsAsEdges = !is.null(es) && es$reactions.as.edges,
                network.hasGenes = !is.null(es$fb.rxn),
                network.usesRpairs = !is.null(es) && es$use.rpairs
            ),
            "showFastHeinzAndMWCS(network.hasReactionsAsNodes);"
        )
    })
    
    output$showModulePanel <- renderJs({
        if (!is.null(esInput())) { return("mp = $('#module-panel'); mp[0].scrollIntoView();")
        }
        # return("mp = $('#module-panel'); mp.hide();")
        return("")
    })
    
    solver <- reactive({
        solverName <- input$solver
        if (solverName == "mwcs") {
            solver <- mwcs.solver(mwcs.path, timeLimit=min(input$mwcsTimeLimit, 120))
        } else if (solverName == "heinz") {
            solver <- heinz.solver(heinz.py, timeLimit=min(input$heinzTimeLimit, 240))
        } else if (solverName == "fastHeinz") {
            solver <- fastHeinz.solver
        } else {
            stop(paste("There is no solver called", solverName))
        }
        solver
    })
    
    rawModuleInput <- reactive({
        input$find
        met.fdr <- isolate(input$metFDR)
        gene.fdr <- isolate(input$geneFDR)
        absent.met.score=isolate(input$absentMetScore)
        absent.rxn.score=isolate(input$absentRxnScore)
        
        es <- isolate(esInput())
        
        if (is.null(es)) {
            return(NULL)
        }

        longProcessStart()
        res <- findModule(es,
                    met.fdr=met.fdr,
                    gene.fdr=gene.fdr,
                    absent.met.score=absent.met.score,
                    absent.rxn.score=absent.rxn.score,
                    solver=isolate(solver()))
        
        if (is.null(res) || length(V(res)) == 0) {
            stop("No module found")
        }
        res$description.string <- paste0(".mp", # min p-value
                                         if (es$reactions.as.edges) ".re" else ".rn",
                                         ".mf=", format(met.fdr, scientific=T),
                                         ".rf=", format(gene.fdr, scientific=T),
                                         ".ams=", absent.met.score,
                                         ".ars=", absent.rxn.score)
        longProcessStop()
        res
    })
    
    moduleInput <- reactive({
        module <- rawModuleInput()
        if (is.null(module)) {
            return(NULL)
        }
        
        es <- isolate(esInput())
        
        if (es$reactions.as.edges) {
            if (isolate(input$useRpairs)) {
                if (input$addTransPairs) {
                    module <- addTransEdges(module, es)
                }
            }
        } else {
            if (input$addMetabolitesForReactions) {
                module <- addMetabolitesForReactions(module, es)
            }
            if (input$addInterconnections) {
                module <- addInterconnections(module, es)
            }
            
            if ("log2FC" %in% list.vertex.attributes(module))
            module <- addNormLogFC(module)
            
            if (input$removeHangingNodes) {
                module <- removeHangingNodes(module)
            }
            
            if (input$simplifyReactionNodes) {
                module <- simplifyReactionNodes(module, es)
            }
            module <- expandReactionNodeAttributesToEdges(module)
        }
            
        module$description.string <- rawModuleInput()$description.string

        layout <- layout.norm(layout.kamada.kawai(module), 0, 1, 0, 1)
        V(module)$x <- layout[,1]
        V(module)$y <- layout[,2]

        module
    })
    
    output$moduleSummary <- reactive({
        module <- moduleInput()
        if (is.null(module)) {
            return("There is no module yet")
        }
        
        vector2html(c(
            "number of nodes" = length(V(module)),
            "number of edges" = length(E(module))
            ))
    })
    
    output$moduleParameters <- reactive({
        m <- NULL
        tryCatch({
            m <- moduleInput()
        }, error=function(e) {})
        makeJsAssignments(
            module.available = !is.null(m)
            )        
    })
    
     output$module <- renderGraph({
         moduleInput()
     })
    
    output$downloadNetwork <- downloadHandler(
        filename = reactive({ paste0("network.", tolower(esInput()$network$organism), ".xgmml") }),
        content = function(file) {
            saveModuleToXgmml(esInput()$subnet, "network", file)
        })
    
    output$downloadModule<- downloadHandler(
        filename = reactive({ paste0("module", moduleInput()$description.string, ".xgmml") }),
        content = function(file) {
            saveModuleToXgmml(moduleInput(), "module", file)
        })
    
    output$GAMVersion <- renderUI({
        p(paste("GAM version:", sessionInfo()$otherPkgs$GAM$Version))
    })
})
