library(markdown)

workPanel <- tagList(
    fixedRow(
        column(12,
            div(
                class="alert alert-info",
                role="alert",
                HTML('If you have any feedback, please send it to <a href="mailto:asergushichev@path.wustl.edu">asergushichev@path.wustl.edu</a>')
                )
            )
        ),
    fixedRow(
        sidebarPanel(width=3,
            selectInput(
                "network",
                label="Select an organism",
                choices=c("Mouse"="mmu",
                  "Human"="hsa"),
                selected="mmu"
                ),
            fileInput("geneDE", "File with gene DE data"),
            fileInput("metDE", "File with metabolic DE data"),
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
            actionButton("preprocess", label="Make network")
            ),
        mainPanel(width=9,
            div(h3("Gene DE data"),
                uiOutput("geneDESummary"),
                uiOutput("geneDETable")),
            div(h3("Metabolite DE data"),
                uiOutput("metDESummary"),
                uiOutput("metDETable")),
            div(class="bottom-buffer",
                h3("Network summary"),
                uiOutput("networkSummary"),
                conditionalPanel("network.available",
                    downloadButton("downloadNetwork", "Download XGMML"))
                )
            )
        ),
    div(id="module-panel", class="row top-buffer",
        column(width=3,
        wellPanel(
            conditionalPanel("network.hasGenes",
                numericInput("geneLogFDR",
                             label=HTML("log<sub>10</sub> FDR for genes"),
                             max=0, min=-100, value=-6, step=1)),
            conditionalPanel("network.hasMets",
                numericInput("metLogFDR", 
                             label=HTML("log<sub>10</sub> FDR for metabolites"),
                             max=0, min=-100, value=-6, step=1),
                numericInput("absentMetScore", 
                             label="Score for absent metabolites",
                             max=0, min=-100, value=-1, step=1)
                ),
            actionButton("find", "Find module"),
            conditionalPanel("network.hasReactionsAsEdges && network.usesRpairs",
                checkboxInput(
                    "addTransPairs", 
                    label="Add trans- edges",
                    value=TRUE)),

            conditionalPanel("network.hasReactionsAsNodes",
                checkboxInput(
                    "addMetabolitesForReactions",
                    label="Add all reactnas for the reactions",
                    value=FALSE),
                checkboxInput(
                    "addInterconnections",
                    label="Add interconnections",
                    value=FALSE),
                checkboxInput(
                    "removeHangingNodes",
                    label="Remove hanging nodes",
                    value=FALSE),
                checkboxInput(
                    "simplifyReactionNodes",
                    label="Simplify reaction nodes",
                    value=FALSE)
                )
            ),
            conditionalPanel("module.available",
                h3("Module summary"),
                uiOutput("moduleSummary"),
                downloadButton("downloadModule")),
            div(id="legend",
                p(),
                div(style="color: red", "Red: log2FC > 0"),
                div(style="color: green", "Green: log2FC < 0"),
                div(style="color: #00acad", "Aqua: log2FC not available")
                )
            ),
        mainPanel(width=9,
            div(id="module", class="graph-output")
            )
        )
    )

helpPanel <- fixedRow(
    mainPanel(includeMarkdown("help.markdown")))

aboutPanel <- fixedRow(
    mainPanel(includeMarkdown("about.markdown")))

shinyUI(
    fixedPage(
        tags$head(
        tags$script(src="d3.v3.min.js"),
        tags$script(src="gam.js"),
        tags$link(rel="stylesheet", href="gam.css"),
        includeScript("ga.js"),
        tags$title("Shiny GAM")
        ),

    includeHTML("misc.xhtml"),
    div(id="updateEsParameters", class="js-output"),

    titlePanel("Shiny GAM: integrated analysis of genes and metabolites"),

    fixedRow(
    column(12,
        tabsetPanel(
            tabPanel("Work", workPanel),
            tabPanel("Help", helpPanel),
            tabPanel("About", aboutPanel)
            )
        ))
))
