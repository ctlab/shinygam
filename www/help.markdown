
This application is a web-interface to the <a href="bioconductor/GAM"
target="_blank">GAM</a> (Gene And Metabolic data analysis) R-package.  GAM
analyzes gene and/or metabolic differential expression (DE) data in a context
of biochemical reactions.


Analysis consists of creating a network describing connections between
metabolites and reactions specific to the data and finding a connected module
that contains the most significant changes.

#### Constructing a network

To construct a network select an organism from the *Select an organism*
dropdown menu.  Only reactions possible in the selected organism are used.
Currently supported organisms are Homo sapiens and Mus musculus. If your
organism is not supported, please contact <a
href="mailto:asergushichev@wustl.edu">asergushichev@wustl.edu</a>.

<!-- :ToDo: add screenshot? -->

The next step is to upload gene and/or metabolite DE data. GAM can be run using
either gene DE data or metabolite DE data or both datasets. Each DE dataset
must be in a separate tab-delimited file. The first line of each file must
contain a header with column names. The columns of each file are:

* "ID" for the RefSeq mRNA transcript ID or Entrez ID (for the gene DE file)
  and the HMDB or KEGG ID (for the metabolite DE file);
* "pval" for the p-value;
* "logFC" for the natural logarithm of the fold-change.

The last column is optional, but we recommend to provide it if possible, as it's
used for colors in graph visualizing. Any other columns will be copied to a
network as node or edge attributes.  We recommend to exclude genes with low
expression prior uploading. Example data for genes and metabolites can be
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
all cross-connections, select the *Use RPAIRs* option (see [KEGG
REACTION](http://www.genome.jp/kegg/reaction/) for details).

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
solved.

This site supports three solvers:

* The [Heinz solver](http://www.mi.fu-berlin.de/w/LiSA/Heinz): it's the
  default solver, as it works for both interpretations of reactions.
* The MWCS solver: this solver can provide more interactive experience because
  when reaching time limit it outputs the best solution found so far.
* The heuristic search: this solver is from R-package BioNet. Added for
  comparison, as it's distributed freely.

As the last two solvers don't support edge scores they are not available when
reactions are interpreted as edges.

The *FDR values*, *Score for absent metabolites*, and *Score for absent
reactions* options control the size of the module. 
Increasing/decreasing the *FDR for reactions* (*FDR for metabolites*) value 
makes adding reactions (metabolites) to a module harder/easier.  We recommend
to start from the default values and then gradually change them depending on
the results.

Click "Find module" button to find a module in the network. The module will
be shown on the right panel.

![Example of a module](img/module.png)

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
[here](/publications/supp_materials/GAM/GAM_VizMap.xml).
