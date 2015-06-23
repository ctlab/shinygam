### Help

This is GAM ("Genes And Metabolites"): a web-service for integrated network analysis of
transcriptional and steady-state metabolomic data focused on identification of
differential metabolic subnetworks most changing between two conditions of
interest.

The analysis consists of:

* Step 1: Creating a network of metabolites and reactions specific 
  to the provided data.
* Step 2: Finding a connected module that contains the most significant changes.

#### Constructing a network

To construct a network select an organism from the *Select an organism*
dropdown menu. Only reactions possible in the selected organism are used.
Currently supported organisms are Homo sapiens and Mus musculus. Please, 
send requests for additional organisms to Alexey Sergushichev at 
asergushichev @ path.wustl.edu.

<!-- :ToDo: add a screenshot? -->

The next step is to upload differential expression (DE) data
for genes and/or metabolites. GAM can be run using
either gene DE data or metabolite DE data or both datasets. Each DE dataset
must be in a separate CSV file (comma-, tab- and space- separated files are 
supported, archived files are supported too). The first line of each file must
contain a header with column names. Files should contain the following columns:

* "ID": RefSeq mRNA transcript ID, Entrez ID or symbol for genes 
  and HMDB or KEGG ID for metbolites.
* "pval": Differential expression p-value (non-adjusted).
* "log2FC": Base 2 logarithm of the fold-change.

**NB: Columns in your files can have somewhat different names and GAM will try
to guess which one to use.** If you favorite DE tool produce names that GAM
can not recognize, please, contact Alexey Sergushichev.

The "log2FC" column is optional, but we recommend to provide it if possible.
It is used for colors in a graph visualization. Any other columns will be copied to
a network as node or edge attributes. Example data can be downloaded:

* <div id="geneDEExample" style="display: inline-block" class="shiny-html-output"></div> for genes and
* <div id="metDEExample" style="display: inline-block" class="shiny-html-output"></div> for metabolites.

After files are uploaded, a file summary is displayed. Verify that the files
were parsed correctly.

![Summary of the example data](img/data_summary.png)

If present, DE for genes is converted to DE for reactions. This is done by 
considering for a reaction all the genes that code any
enzyme that takes part in this reaction. 
Gene with the minimal p-value is selected and its p-value is assigned as the
reaction p-value. All reactions without p-values are discarded as having no
expressed enzymes. We recommend to exclude genes with low expression prior
uploading.

Reactions in the network can be interpreted as either edges or nodes. Select
one of the interpretations.

**Reactions as edges.** If reactions are interpreted as edges, if a pair of
metabolites has multiple reaction connections between them, only the reaction
with the minimal p-value is kept. All the remaining reactions and metabolites
that take part in these reactions are combined into a network as edges and
nodes respectively. To only consider substrate-product pairs
that make up a main type KEGG reaction pair (RPAIR), as opposed to
all substrate-product pairs, select the *Use RPAIRs* option.

**Reactions as nodes.** If reactions are interpreted as nodes, then both
metabolites and reactions are added as nodes to the network. Edges are added
between a metabolite and a reaction if the metabolite takes part in this
reaction. We recommend to collapse groups of reactions that have at least one common
metabolite and the same most significant gene into single nodes with
*Collapse reactions* option. It removes artificial biases created by clusters
of similar reactions.

Generally, all of these options will lead to similar results. We recommend to
use the default option values as this makes the network simpler and the
analysis faster.

Click *Step 1: Make network* to finish this step.

#### Finding a module

After making a network, you can find a connected module that contains the most
significantly changed genes and metabolites. Internally, this is done by first
scoring nodes and edges based on their p-values such that positive scores
correspond to significant p-values and negative scores correspond to
insignificant changes. Then the problem of finding a connected subgraph with
maximum summary weight (maximum-weight connected subgraph, MWCS problem) is
solved.

The *FDR values*, *Score for absent metabolites* options control the size of the module. 
Increasing/decreasing the *FDR for reactions* (*FDR for metabolites*) value 
makes adding reactions (metabolites) to a module easier/harder.  We recommend
to start from the default values and then gradually change them depending on
the results. Default FDR values are generated so that the module will
a size of ~100 reactions (*Autogenerate FDRs* button).

There are three solvers available to solve MWCS instance: heinz, heinz 2 and gmwcs.
Heinz solver (Dittrich at al. 2008) is called when the option *Try to solve to optimality*
selected. Solving to optimality can take a long time, so we provide two other 
solvers that can be run for a predefined amount of time. These are Heinz 2 
(El-Kebir&Klau 2014, https://software.cwi.nl/cwisoftware/software/heinz), which is used when edges are not scored, and gmwcs 
(https://github.com/alexloboda/gmwcs-solver), which is used when edges are scored.


Click *Step 2: Find module* button to find a module in the network. The module will
be shown on the right panel.

![Example of a module](img/module.png)

#### Graph legend

We use the following scheme:

* Red nodes and edges are up-regulated (*log2FC* > 0).
* Green nodes and edges are down-regulated (*log2FC* < 0).
* Blue nodes and edges don't have *log2FC* values.
* Bigger size of nodes and width of edges means lower p-values.
* Dashed edges are trans-RPAIRs.
* Circle nodes represent metabolites, square nodes represent reactions.

#### Post-processing

There are couple of post-processing steps that are available. If in the network
reactions are edges and RPAIRs are used, then trans-connections can be added.

When reactions are nodes, metabolites that take part in at least two reactions
from the module can be added into the module.

#### Saving the module

You can download the module in an PDF or XGMML format by clicking the corresponding
button. XGMML file can be imported into Cytoscape. You can also download GAM's
VizMap style for Cytoscape 
<a id="downloadVizMapInHelp" class="shiny-download-link" href="" target="_blank">here</a>.
