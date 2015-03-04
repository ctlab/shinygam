### Help

This application is a web-interface to the GAM (Gene And Metabolite data analysis) R-package.  GAM
analyzes gene and/or metabolic differential expression (DE) data in a context
of biochemical reactions.


Analysis consists of creating a network describing connections between
metabolites and reactions specific to the data and finding a connected module
that contains the most significant changes.

#### Constructing a network

To construct a network select an organism from the *Select an organism*
dropdown menu.  Only reactions possible in the selected organism are used.
Currently supported organisms are Homo sapiens and Mus musculus. If your
organism is not supported, please contact Alexey Sergushichev at asergushichev@path.wustl.edu.

<!-- :ToDo: add a screenshot? -->

The next step is to upload gene and/or metabolite DE data. GAM can be run using
either gene DE data or metabolite DE data or both datasets. Each DE dataset
must be in a separate CSV file (comma-, tab- and space- separated files are 
supported, archived files are supported too). The first line of each file must
contain a header with column names. Files should contain the following columns:

* "ID": RefSeq mRNA transcript ID, Entrez ID or symbol for genes and
  and HMDB or KEGG ID for metbolites.
* "pval" for the p-value;
* "log2FC" for the base 2 logarithm of the fold-change.

**NB: Columns in your files can have a little different names and GAM will try
to guess which one to use.** If you favorite DE tool produce names that GAM
can not recoginze, please tell me and I will try to add it.

The "log2FC" column is optional, but we recommend to provide it if possible, as it's
used for colors in graph visualizing. Any other columns will be copied to a
network as node or edge attributes.  Example data for genes and metabolites can be
downloaded [here](/publications/supp_materials/GAM/).

After files are uploaded, a file summary is displayed. Verify that the files
were parsed correctly.

![Summary of the example data](img/data_summary.png)

If gene DE data was uploaded, the gene DE data is converted to reaction DE
data. This is done by considering for a reaction all the genes that code any
enzyme that takes part in this reaction. P-values of these genes are compared.
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
nodes respectively. To only consider cross-connections between substrates and
products in reactions that make up a KEGG reaction pair (RPAIR) as opposed to
all cross-connections, select the *Use RPAIRs* option, see [KEGG
REACTION](http://www.genome.jp/kegg/reaction/) for details.

**Reactions as nodes.** If reactions are interpreted as nodes, then both
metabolites and reactions are added as nodes to the network. Edges are added
between a metabolite and a reaction if the metabolite takes part in this
reaction. To collapse groups of reactions that have at least one common
metabolite and the same most significant gene into single nodes, select the
*Collapse reactions* option.

Generally, all of these options will lead to similar results. We recommend to
use the default option values as this makes the network simpler and the
analysis faster.

#### Finding a module

After making a network, you can find a connected module that contains the most
significantly changed genes and reactions. Internally, this is done by first
scoring nodes and edges based on their p-values such that positive scores
correspond to significant p-values and negative scores correspond to
insignificant changes. Then the problem of finding a connected subgraph with
maximum summary weight (maximum-weight connected subgraph, MWCS problem) is
solved. We use the [Heinz 2 solver](https://software.cwi.nl/cwisoftware/software/heinz).

The *FDR values*, *Score for absent metabolites*, and *Score for absent
reactions* options control the size of the module. 
Increasing/decreasing the *FDR for reactions* (*FDR for metabolites*) value 
makes adding reactions (metabolites) to a module harder/easier.  We recommend
to start from the default values and then gradually change them depending on
the results.

Click "Find module" button to find a module in the network. The module will
be shown on the right panel.

![Example of a module](img/module.png)

#### Graph legend

We use the following scheme:

* Red nodes and edges are up-regulated.
* Green nodes and edges are down-regulated.
* Blue nodes and edges don't have *log2FC* values.
* Bigger size of nodes and width of edges means lower p-values.
* Dashed edges are trans-RPAIRs.

#### Post-processing

There are couple of post-processing steps that are available. If in the network
reactions are edges and RPAIRs are used, then trans-connections can be added.
When reactions are nodes, the following operations can be applied:

* Adding all metabolites that are part of any reaction in the module.
* Adding all reactions that aren't present in the module but are directly
  connected to two metabolites from the module.
* Removing all metabolites with no data and only one connection in the module.
* Replacing reaction nodes with edges for reactions that connect only two
  metabolites in the module and these two metabolites are on different sides of
  the reaction.

#### Saving the module

You can download the module in an XGMML format by clicking the *Download XGMML*
button. This file can be imported into Cytoscape. You can also download GAM's
VizMap style for Cytoscape 
<a id="downloadVizMap" class="shiny-download-link" href="" target="_blank">here</a>.
