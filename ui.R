library(markdown)

myActionButton <- function (inputId, label, icon = NULL, ...) 
{
    tags$button(id = inputId, type = "button", class = "btn action-button", 
        list(icon, label), ...)
}

mySidebarPanel <- function(...) {
    div(class="sidebar col-sm-3", tags$form(class="well", ...))
}

myMainPanel <- function(...) {
    div(class="mainPanel col-sm-9", ...)
}

workPanel <- tagList(
    fixedRow(
        column(12,
            div(
                class="alert alert-info",
                role="alert",
                HTML('Like Shiny GAM? Check out <a href="https://artyomovlab.wustl.edu/phantasus">Phantasus</a> where you can <a href="https://artyomovlab.wustl.edu/phantasus/phantasus-tutorial.html">do differential expression and submit the results to Shiny GAM</a>')
                )
            )
        ),
    fixedRow(
        mySidebarPanel(
            actionButton("resetInput", label="Reset all", 
                           onclick=paste(
                                         'resetFileInput("geneDE"); resetFileInput("metDE")',
                                         '$("#loadExampleGeneDE").attr("checked", false).trigger("change")',
                                         '$("#loadExampleMetDE").attr("checked", false).trigger("change")',
                                         '$("#preprocess").click(); $("#find").click()',
                                         sep=";")
                           ),
            checkboxInput(
                "loadExampleGeneDE",
                label="Example DE for genes",
                value=FALSE),
            checkboxInput(
                "loadExampleMetDE",
                label="Example DE for metabolites",
                value=FALSE),
            conditionalPanel("input.loadExampleGeneDE || input.loadExampleMetDE",
                p("Organism: Mouse")
                ),
            conditionalPanel("!input.loadExampleGeneDE && !input.loadExampleMetDE",
                selectInput(
                    "network",
                    label="Select an organism",
                    choices=c("Mouse"="mmu",
                      "Human"="hsa",
                      "Arabidopsis"="ath",
                      "Yeast"="sce"),
                    selected="mmu"
                    ),
                fileInput("geneDE", "File with DE for genes"),
                fileInput("metDE", "File with DE for metabolites")
                ),
            #uiOutput("reactionsAsHolder"),
            selectInput("reactionsAs", 
                        label="Interpret reactions as",
                        c("edges"="edges", "nodes"="nodes"),
                        selected="nodes"),
            conditionalPanel("input.reactionsAs == 'nodes'",
                checkboxInput(
                    "collapseReactions",
                    label="Collapse reactions",
                    value=TRUE)),
            conditionalPanel("input.reactionsAs == 'edges'",
                checkboxInput(
                    "useRpairs",
                    label="Use RPAIRs",
                    value=TRUE)),
            conditionalPanel("false",
                checkboxInput(
                    "autoFindModule",
                    label="Set FDRs & run step 2 automatically",
                    value=TRUE),
                actionButton("preprocess", label="Step 1: Make network")
            ),
            myActionButton("runAll",   label="Run step 1, autogenerate FDRs and run step 2", 
                          onclick='$("#autoFindModule")[0].checked=true; $("#autoFindModule").trigger("change"); $("#preprocess").trigger("click")',
                          disabled=""),
            div("or", style="text-align: center"),
            myActionButton("runStep1", label="Step 1: Make network", 
                          onclick='$("#autoFindModule")[0].checked=false; $("#autoFindModule").trigger("change"); $("#preprocess").trigger("click")',
                          disabled="")
            ),
        myMainPanel(
            div(class="DEBlock",
                h3("Differential expression for genes"),
                uiOutput("geneDESummary"),
                uiOutput("geneDETable"),
                uiOutput("geneDENotMapped"),
                uiOutput("geneDENotMappedTable")),
            div(class="DEBlock",
                h3("Differential expression for metabolites"),
                uiOutput("metDESummary"),
                uiOutput("metDETable"),
                uiOutput("metDENotMapped"),
                uiOutput("metDENotMappedTable")),
            div(class="bottom-buffer",
                h3("Network summary"),
                uiOutput("networkSummary"),
                conditionalPanel("network.available",
                    downloadButton("downloadNetwork", "Download XGMML"))
                )
            )
        ),
    div(id="module-panel", class="row top-buffer",
        div(class="sidebar col-sm-3",
            wellPanel(
                conditionalPanel("network.hasGenes",
                    numericInput("geneLogFDR",
                                 label=HTML("log<sub>10</sub> FDR for genes"),
                                 max=0, min=-100, value=-6, step=0.5)),
                conditionalPanel("network.hasMets",
                    numericInput("metLogFDR", 
                                 label=HTML("log<sub>10</sub> FDR for metabolites"),
                                 max=0, min=-100, value=-6, step=0.5),
                    numericInput("absentMetScore", 
                                 label="Score for absent metabolites",
                                 max=0, min=-100, value=-20, step=1)
                    ),
                myActionButton("resetFDRs", "Autogenerate FDRs", disabled=""),
                checkboxInput(
                    "solveToOptimality", 
                    label="Try to solve to optimality",
                    value=FALSE),
                conditionalPanel("network.available",
                    uiOutput("solverString")
                ),
                myActionButton("find", "Step 2: Find module", disabled=""),
                conditionalPanel("network.available",
                    conditionalPanel("network.hasReactionsAsEdges && network.usesRpairs",
                        checkboxInput(
                            "addTransPairs", 
                            label="Add trans- edges",
                            value=TRUE)),

                    conditionalPanel("network.hasReactionsAsNodes",
                        checkboxInput(
                            "addCommonMetabolites",
                            label="Add common metabolites",
                            value=FALSE)
                        )
                    )
                ),
            conditionalPanel("module.available",
                h3("Module summary"),
                uiOutput("moduleSummary"),
                downloadButton("downloadPDF", "PDF"),
                downloadButton("downloadModule", "XGMML"),
                downloadButton("downloadXlsx", "XLSX")
                ),
            div(id="legend",
                p(),
                div(img(src="img/log2FCscale.svg"))
                ),
            p(),
            downloadButton("downloadVizMap", "Cytoscape VizMap")
            ),
        myMainPanel(
            #div(id="module", class="graph-output")
            uiOutput("module")
            )
        )
    )

helpPanel <- fixedRow(
    mainPanel(id="helpPanel",
              includeMarkdown("help.markdown")))

aboutPanel <- fixedRow(
    mainPanel(includeMarkdown("about.markdown")))

shinyUI(
    fluidPage(
        tags$head(
        tags$script(src="svg-pan-zoom.min.js"),
        tags$script(src="d3.v3.min.js"),
        #tags$script(src="d3.v3.js"),
        tags$script(src="gam.js"),
        tags$link(rel="stylesheet", href="gam.css"),
        includeScript("ga.js"),
        tags$title("Shiny GAM")
        ),

    includeHTML("misc.xhtml"),
    div(id="updateEsParameters", class="js-output"),

    titlePanel("Shiny GAM: integrated analysis of genes and metabolites"),
    div(id="initialized", class="js-output"),

    fixedRow(
    column(12,
        tabsetPanel(
            tabPanel("Work", workPanel),
            tabPanel("Help", helpPanel),
            tabPanel("About", aboutPanel)
            )
        ))
))
