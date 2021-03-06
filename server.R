library(shiny)
library(data.table)
library(igraph)
library(GAM)
library(GAM.db)
library(GAM.networks)
library(RCurl)
library(parallel)
#library(xlsx)
library(pryr)
library(logging)
addHandler(writeToFile, file="/dev/stderr", level='DEBUG')

"%o%" <- pryr::compose

options(shiny.error=traceback)
#options(shiny.trace=TRUE)
options(shiny.fullstacktrace=TRUE)

source("./functions.R")

#data("met.id.map")
#data("kegg.human.network")
#data("kegg.mouse.network")
#data("kegg.arabidopsis.network")
#data("kegg.yeast.network")

uploaderPath <- "/mnt/data/upload"

networks <- list(
    "mmu"="kegg.mouse.network",
    "hsa"="kegg.human.network",
    "ath"="kegg.arabidopsis.network",
    "sce"="kegg.yeast.network"
    )



heinz2 <- "/usr/local/lib/heinz2/heinz21"
h.solver <- heinz.solver("/usr/local/lib/heinz/heinz.py", timeLimit=4*60)
attr(h.solver, "description") <- "Heinz (time limit = 4m)"

h2.solver <- heinz21.solver(heinz2, timeLimit=30, nthreads=detectCores())
attr(h2.solver, "description") <- "Heinz2 (time limit = 30s)"
g.solver <- gmwcs.solver("gmwcs", timeLimit=30, nthreads=detectCores())
attr(g.solver, "description") <- "gmwcs (time limit = 30s)"



example.gene.de.path <- "http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"
example.met.de.path <- "http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv"


# not sure if works
example.gene.de <- force(as.data.table(read.table(text=getURL(example.gene.de.path), stringsAsFactors=FALSE, header=1)))
attr(example.gene.de, "name") <- basename(example.gene.de.path)
example.met.de <- force(as.data.table(read.table(text=getURL(example.met.de.path), stringsAsFactors=FALSE, header=1)))
attr(example.met.de, "name") <- basename(example.met.de.path)


shinyServer(function(input, output, session) {
    longProcessStart <- function() {
        session$sendCustomMessage(type='showWaitMessage', list(value=T))
    }

    longProcessStop <- function() {
        session$sendCustomMessage(type='showWaitMessage', list(value=F))
    }

    output$initialized <- renderJs('$("#initializing").hide()')

    queryValues <- reactiveValues()

    observe({
        query <- parseQueryString(session$clientData$url_search)
        if ("organism" %in% names(query)) {
            updateSelectInput(session, "network", selected=query$organism)
            #values$queryOrganism <- query$organism
        } 
        
        if ("geneDE_key" %in% names(query)) {
            geneDE_key <- gsub("[^a-z0-9]", "", query$geneDE_key)
            loginfo("found key: %s", geneDE_key)
            if (geneDE_key == "") {
                geneDE_key <- NULL
            }
            queryValues$geneDE_key <- geneDE_key
        }
     })

    loadExample <- reactive({
        input$loadExampleGeneDE || input$loadExampleMetDE
    })

    getNetwork <- reactive({
        net <- if (loadExample()) {
            "kegg.mouse.network"
        } else {
            networks[[input$network]]
        }

        GAM:::lazyData(net)
        res <- get(net)
        # :ToDo: it's a hack
        res$rxn2name$name <- ""
        res
    })

    geneDEInput <- reactive({
        if (loadExample()) {
            if (input$loadExampleGeneDE) {
                return(example.gene.de)
            }
            return(NULL)
        }
        
        if (!is.null(input$geneDE) && !is(input$geneDE, "data.frame")) {
            # Value was reset
            return(NULL)
        }
        
        if (!is.null(input$geneDE)) {
            loginfo("GeneDE file:")
            loginfo(capture.output(str(input$geneDE)))
            loginfo("reading gene.de: %s", input$geneDE$name)
            loginfo("      from file: %s", input$geneDE$datapath)
            
            path <- input$geneDE$datapath
            deName <- input$geneDE$name
        } else if (!is.null(queryValues$geneDE_key)) {
            # User has not uploaded a file yet but provided value by key
            path <- file.path(uploaderPath, queryValues$geneDE_key)
            loginfo("GeneDE from key, path: %s", path)
            if (!file.exists(path)) {
                return(NULL)
            }
            deName <- queryValues$geneDE_key
        } else {
            # User has not uploaded a file yet and we don't have a key
            return(NULL)
        }
        
        network <- getNetwork()
        gene.id.map <- network$gene.id.map
        
        res <- read.table.smart.de.gene(path, idsList=gene.id.map)
        logdebug(capture.output(str(res)))
        if (!all(necessary.de.fields %in% names(res))) {
            loginfo("not all fields in DE file: %s", input$geneDE$datapath)
            if (grepl("xlsx?$", input$geneDE$name)) {
                stop("We do not support excel files yet, please, use tab-separated files instead")
            } else{
                stop(paste0("Genomic differential expression data should contain at least these fields: ", 
                            paste(necessary.de.fields, collapse=", ")))
            }
        }
        attr(res, "name") <- deName
        res
    })
    
    geneIdsType <- reactive({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        network <- getNetwork()
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
    
    notMappedGenes <- reactive({
        network <- getNetwork()
        gene.id.map <- network$gene.id.map
        geneIT <- geneIdsType()

        if (is.null(geneIT)) {
            return(NULL)
        }

        if (geneIT == network$gene.ids) {
            return(NULL)
        }
        notMapped <- setdiff(geneDEInput()$ID, gene.id.map[[geneIT]])
        notMapped
    })

    output$geneDENotMapped <- renderUI({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        notMapped <- notMappedGenes()
        network <- getNetwork()

        div(
            p(sprintf("Not mapped to %s: %s", network$gene.ids, length(notMapped))),
            if (length(notMapped) > 0) { 
                p("Top unmapped genes:",
                  a("show",
                    id="geneDENotMappedTableShowBtn",
                    href="javascript:void()",
                    onclick='$("#geneDENotMappedTable").show();$("#geneDENotMappedTableShowBtn").hide();$("#geneDENotMappedTableHideBtn").show()'),
                  a("hide",
                    id="geneDENotMappedTableHideBtn",
                    href="javascript:void()",
                    onclick='$("#geneDENotMappedTable").hide();$("#geneDENotMappedTableShowBtn").show();$("#geneDENotMappedTableHideBtn").hide()',
                    style="display:none"),
                  tag("script", '$("#geneDENotMappedTable").hide()')
                  ) 
            } else NULL)
    })

    output$geneDENotMappedTable <- renderTable({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        notMapped <- notMappedGenes()
        if (length(notMapped) == 0) {
            return(NULL)
        }

        data <- data[order(pval)]
        data <- data[ID %in% notMapped]
        format(as.data.frame(head(data, n=20)), digits=3)
    })
    
    metDEInput <- reactive({
        if (loadExample()) {
            if (input$loadExampleMetDE) {
                return(example.met.de)
            }
            return(NULL)
        }


        if (is.null(input$metDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        if (!is(input$metDE, "data.frame")) {
            # Value was reset
            return(NULL)
        }

        loginfo("MetDE file:")
        loginfo(capture.output(str(input$metDE)))
        loginfo("reading met.de: %s", input$metDE$name)
        loginfo("     from file: %s", input$metDE$datapath)

        res <- read.table.smart.de.met(input$metDE$datapath)
        logdebug(capture.output(str(res)))
        if (!all(necessary.de.fields %in% names(res))) {
            loginfo("not all fields in DE file: %s", input$metDE$datapath)
            if (grepl("xlsx?$", input$metDE$name)) {
                stop("We do not support excel files yet, please, use tab-separated files instead")
            } else {
                stop(paste0("Metabolic differential expression data should contain at least these fields: ", 
                            paste(necessary.de.fields, collapse=", ")))
            }
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
    
    notMappedMets <- reactive({
        network <- getNetwork()
        metIT <- metIdsType()

        if (is.null(metIT)) {
            return(NULL)
        }

        if (metIT == network$met.ids) {
            return(NULL)
        }
        GAM:::lazyData("met.id.map")
        notMapped <- setdiff(metDEInput()$ID, met.id.map[[metIT]])
    })

    output$metDENotMapped <- renderUI({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        notMapped <- notMappedMets()
        network <- getNetwork()

        div(
            p(sprintf("Not mapped to %s: %s", network$met.ids, length(notMapped))),
            if (length(notMapped) > 0) { 
                p("Top unmapped metabolites:",
                  a("show",
                    id="metDENotMappedTableShowBtn",
                    href="javascript:void()",
                    onclick='$("#metDENotMappedTable").show();$("#metDENotMappedTableShowBtn").hide();$("#metDENotMappedTableHideBtn").show()'),
                  a("hide",
                    id="metDENotMappedTableHideBtn",
                    href="javascript:void()",
                    onclick='$("#metDENotMappedTable").hide();$("#metDENotMappedTableShowBtn").show();$("#metDENotMappedTableHideBtn").hide()',
                    style="display:none"),
                  tag("script", '$("#metDENotMappedTable").hide()')
                  ) 
            } else NULL)
    })

    output$metDENotMappedTable <- renderTable({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        notMapped <- notMappedMets()
        if (length(notMapped) == 0) {
            return(NULL)
        }

        data <- data[order(pval)]
        data <- data[ID %in% notMapped]
        format(as.data.frame(head(data, n=20)), digits=3)
    })
    
    output$updateEsParameters <- renderJs({        
        selected <- if (!is.null(metDEInput())) "edges" else "nodes"
        return(sprintf("$('#reactionsAs')[0].selectize.setValue('%s')", selected))
    })

    experimentTag <- reactive({
        geneData <- geneDEInput()
        metData <- metDEInput()
        name <- if (!is.null(geneData)) attr(geneData, "name") else attr(metData, "name")
        tag <- name
        tag <- gsub("\\.gz$", "", tag)
        tag <- gsub("\\.([ct]sv|txt)$", "", tag)
        tag
    })
# 
#     output$reactionsAsHolder <- renderUI({
#         gene.de <- geneDEInput()
# 
#         met.de <- metDEInput()
# 
#         selected <- if (!is.null(met.de)) "edges" else "nodes"
# 
#         selectInput("reactionsAs", 
#                     label="Interpret reactions as",
#                       c("edges"="edges", "nodes"="nodes"),
#                       selected=selected)
#     })


    esInput <- reactive({
        input$preprocess
        loginfo("Preprocessing")
        #input$runAll
        network <- isolate(getNetwork())
        gene.de <- isolate(geneDEInput())
        gene.ids <- isolate(geneIdsType())

        met.de <- isolate(metDEInput())
        met.ids <- isolate(metIdsType())
        tag <- isolate(experimentTag())

        if (is.null(gene.de) && is.null(met.de)) {
            return(NULL)
        }

        longProcessStart()

        tryCatch({
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
            attr(es, "tag") <- tag
            es
        }, finally=longProcessStop())
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
    
    output$FDRParameters <- reactive({
        res <- sprintf("%s;", input$resetFDRs)
        es <- NULL
        tryCatch({
            es <- isolate(esInput())
        }, error=function(e) {})
        
        if (!is.null(es)) {
            res <- paste0(res, generateFDRs(es))
        } else {
            res <- ""
        }
        res 
    })

    output$networkParameters <- reactive({
        es <- NULL
        tryCatch({
            es <- esInput()
        }, error=function(e) {})
        
        if (is.null(es)) {
            return("")
        }

        res <- paste0(
            makeJsAssignments(
                network.available = TRUE,
                network.hasReactionsAsNodes = !es$reactions.as.edges,
                network.hasReactionsAsEdges = es$reactions.as.edges,
                network.hasGenes = !is.null(es$fb.rxn),
                network.hasMets = !is.null(es$fb.met),
                network.usesRpairs = es$use.rpairs
            )
        )
        
        if (isolate(input$autoFindModule)) {
            res <- paste0(res, generateFDRs(es))
            res <- paste0(res, '$("#find").trigger("click");')
        }

        res <- paste0(res, '$("#find").removeAttr("disabled").addClass("btn-default");')
        res <- paste0(res, '$("#resetFDRs").removeAttr("disabled").addClass("btn-default");')
        res
    })

    output$enableMakeNetwork <- renderJs({
        res <- ""
        canRun <- FALSE
        tryCatch({
            geneDE <- geneDEInput()
            gIT <- geneIdsType()
            metDE <- metDEInput()
            mIT <- metIdsType()

            canRun <- !is.null(geneDE) || !is.null(metDE)
        }, error=function(e) {
            # if anything happened, not running
        })

        if (canRun) {
            res <- paste0(res, '$("#runStep1").removeAttr("disabled").addClass("btn-default");')
            res <- paste0(res, '$("#runAll").removeAttr("disabled").addClass("btn-default");')
        } else {
            res <- paste0(res, '$("#runStep1").attr("disabled", "disabled");')
            res <- paste0(res, '$("#runAll").attr("disabled", "disabled");')
        }
        res
    })
    
    output$showModulePanel <- renderJs({
        if (!is.null(esInput())) { return("mp = $('#module-panel'); mp[0].scrollIntoView();")
        }
        # return("mp = $('#module-panel'); mp.hide();")
        return("")
    })


    metFDR <- reactive({
        10^input$metLogFDR
    })

    geneFDR <- reactive({
        10^input$geneLogFDR
    })

    getSolver <- reactive({
        es <- esInput()
        if (is.null(es)) {
            return(NULL)
        }

        if (input$solveToOptimality) {
            h.solver
        } else if (es$reactions.as.edges && !is.null(es$fb.rxn)) {
            g.solver
        } else {
            h2.solver
        }
    })

    output$solverString <- reactive({
        es <- esInput()
        solver <- getSolver()
        if (!is.null(solver)) {
            sprintf("Solver: %s", attr(solver, "description"))
        } else {
            ""
        }
    })

    esScoredInput <- reactive({
        input$find
        met.fdr <- isolate(metFDR())
        gene.fdr <- isolate(geneFDR())
        absent.met.score <- isolate(input$absentMetScore)
        #absent.rxn.score <- isolate(input$absentRxnScore)
        #absent.rxn.score <- 0
        
        es <- isolate(esInput())
        
        if (is.null(es)) {
            return(NULL)
        }

        longProcessStart()
        loginfo(paste0(attr(es, "tag"),".mp", # min p-value
                       if (es$reactions.as.edges) ".re" else ".rn",
                       ".mf=", format(log10(met.fdr), digist=2),
                       ".rf=", format(log10(gene.fdr), digist=2),
                       ".ams=", absent.met.score
                       #, ".ars=", absent.rxn.score
                       ))

        res <- scoreNetwork(es,
                            met.fdr=met.fdr,
                            rxn.fdr=gene.fdr,
                            absent.met.score=absent.met.score,
                            #absent.rxn.score=absent.rxn.score,
                            met.score=-0.01,
                            rxn.score=-0.01)

        res$description.string <- paste0(attr(es, "tag"), 
                                         if (es$reactions.as.edges) ".re" else ".rn",
                                         ".mf=", format(log10(met.fdr), digist=2),
                                         ".rf=", format(log10(gene.fdr), digist=2),
                                         ".ams=", absent.met.score
                                         #, ".ars=", absent.rxn.score
                                         )
        res
    })
    
    rawModuleInput <- reactive({
        input$find
        
        esScored <- isolate(esScoredInput())
        
        if (is.null(esScored)) {
            return(NULL)
        }

        longProcessStart()
        tryCatch({
            solver <- isolate(getSolver())

            #res <- induced.subgraph(esScored$subnet.scored, V(esScored$subnet.scored)[adj(adj(1))])

            res <- findModule(esScored,
                        solver=solver)
            
            if (is.null(res) || length(V(res)) == 0) {
                stop("No module found")
            }
            res$description.string <- esScored$description.string
            res
        }, finally=longProcessStop())
    })
    
    moduleInput <- reactive({
        module <- rawModuleInput()
        if (is.null(module)) {
            return(NULL)
        }

        # for consistency
        module <- remove.vertex.attribute(module, "score")
        module <- remove.edge.attribute(module, "score")
        
        es <- isolate(esInput())
        
        if (es$reactions.as.edges) {
            if (isolate(input$useRpairs)) {
                if (input$addTransPairs) {
                    module <- addTransEdges(module, es)
                }
            }
        } else {
            if (input$addCommonMetabolites) {
                module <- addMetabolitesForReactions(module, es)
                module <- removeHangingNodes(module)
            }
        }

            
        module$description.string <- rawModuleInput()$description.string

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
    
    #output$module <- renderGraph({
    #    moduleInput()
    #})

    output$module <- renderUI({
        sf <- svgFile()

        if (!is.null(sf)) {
            HTML(readLines(sf), "<script>var panZoomModule = svgPanZoom('#module svg');</script>")
        } else {
            HTML("")
        }
    })
    
    output$downloadNetwork <- downloadHandler(
        filename = reactive({ paste0("network.", tolower(esInput()$network$organism), ".xgmml") }),
        content = function(file) {
            saveModuleToXgmml(esInput()$subnet, file=file, name=tolower(esInput()$network$organism))
        })
    
    output$downloadModule <- downloadHandler(
        filename = reactive({ paste0(moduleInput()$description.string, ".xgmml") }),
        content = function(file) {
            saveModuleToXgmml(moduleInput(), file=file, moduleInput()$description.string)
        })

    dotFile <- reactive({
        m <- moduleInput()
        if (is.null(m)) {
            return(NULL)
        }
        longProcessStart()
        res <- tempfile(pattern="module", fileext=".dot")
        saveModuleToDot(m, file=res, name=m$description.string)
        longProcessStop()
        res
    })

    svgFile <- reactive({
        df <- dotFile()
        if (is.null(df)) {
            return(NULL)
        }
        res <- paste0(df, ".svg")
        system2("neato", c("-Tsvg", 
                           "-o", res,
                           df), stderr=NULL)
        res
    })

    pdfFile <- reactive({
        df <- dotFile()
        if (is.null(df)) {
            return(NULL)
        }
        res <- paste0(df, ".pdf")
        system2("neato", c("-Tpdf", 
                           "-o", res,
                           df), stderr=NULL)
        res
    })

    output$downloadXlsx <- downloadHandler(
        filename = reactive({ paste0(moduleInput()$description.string, ".xlsx") }),
        content = function(file) {
            module <- moduleInput()
            es <- isolate(esScoredInput())
            stopifnot(require(xlsx))

            wb <- createWorkbook()

            vTable <- data.table(get.vertex.attributes(module))
            eTable <- data.table(get.edge.attributes(module, include.ends=T))
            if (es$reactions.as.edges) {
                metTable <- vTable
                metTable[, nodeType := NULL]
                rxnTable <- eTable
            } else {
                metTable <- vTable[nodeType == "met",]
                metTable[, nodeType := NULL]
                rxnTable <- vTable[nodeType == "rxn",]
                rxnTable[, nodeType := NULL]
            }

            metTable <- removeNAColumns(metTable)
            rxnTable <- removeNAColumns(rxnTable)

            addDataFrame(metTable, createSheet(wb, "metabolites"), row.names=F)
            addDataFrame(rxnTable, createSheet(wb, "reactions"), row.names=F)

            metInModule <- metTable$name
            geneInModule <- rxnTable$origin

            saveWorkbook(wb, file)
        })
    
    output$downloadPDF <- downloadHandler(
        filename = reactive({ paste0(moduleInput()$description.string, ".pdf") }),
        content = function(file) {
            file.copy(pdfFile(), file) 
        })

    output$downloadDot <- downloadHandler(

        filename = reactive({ paste0(moduleInput()$description.string, ".dot") }),
        content = function(file) {
            file.copy(dotFile(), file) 
        })
    
    output$downloadVizMap <- downloadHandler(
        filename = "GAM_VizMap.xml",
        content = function(file) {
            file.copy(
                from=system.file("GAM_VizMap.xml", package="GAM"),
                to=file)
        })

    output$GAMVersion <- renderUI({
        p(paste("GAM version:", sessionInfo()$otherPkgs$GAM$Revision))
    })

    output$geneDEExample <- renderUI({
        a("here", href=example.gene.de.path, target="_blank")
    })

    output$metDEExample <- renderUI({
        a("here", href=example.met.de.path, target="_blank")
    })
    

    output$downloadVizMapInHelp <- downloadHandler(
        filename = "GAM_VizMap.xml",
        content = function(file) {
            file.copy(
                from=system.file("GAM_VizMap.xml", package="GAM"),
                to=file)
        })
})
