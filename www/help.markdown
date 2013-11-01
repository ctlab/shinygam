### Quick start

Reload the site, so that all parameters will be set to default.

Download example
data for <a href="/publications/supp_materials/GAM/gene.de.M1.M2.tsv">genes</a>
and <a href="/publications/supp_materials/GAM/met.de.M1.M2.tsv">metabolites</a>.

Upload downloaded files as TSV-files with gene and metabolite p-values respectively.
You'll see a summary of the data:

<div class="row">
<div class="span8">
<div class="example-panel">
<div class="row">
<div class="span4">
<h3>Genomic data</h3>
<div class="shiny-html-output shiny-bound-output"><ul>
<li> length :  14276 </li>
<li> ID type :  RefSeq </li>
</ul>
</div>
<div class="shiny-html-output shiny-bound-output"><!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<table class="data table table-bordered table-condensed">
<tbody><tr> <th>  </th> <th> ID </th> <th> pval </th> <th> logFC </th>  </tr>
<tr> <td align="right"> 1 </td> <td> NM_021274 </td> <td> 1.489938e-50 </td> <td> -12.292868 </td> </tr>
<tr> <td align="right"> 2 </td> <td> NM_013653 </td> <td> 1.098282e-46 </td> <td> -10.498761 </td> </tr>
<tr> <td align="right"> 3 </td> <td> NM_010501 </td> <td> 3.872073e-46 </td> <td> -12.144405 </td> </tr>
<tr> <td align="right"> 4 </td> <td> NM_001177351 </td> <td> 1.126691e-40 </td> <td>  -9.145330 </td> </tr>
<tr> <td align="right"> 5 </td> <td> NM_008332 </td> <td> 4.930544e-40 </td> <td>  -9.710362 </td> </tr>
<tr> <td align="right"> 6 </td> <td> NM_145137 </td> <td> 9.538526e-37 </td> <td>  12.364258 </td> </tr>
</tbody></table></div>
</div>
<div class="span4">
<h3>Metabolic data</h3>
<div class="shiny-html-output shiny-bound-output"><ul>
<li> length :  2119 </li>
<li> ID type :  HMDB </li>
</ul>
</div>
<div class="shiny-html-output shiny-bound-output"><!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<table class="data table table-bordered table-condensed">
<tbody><tr> <th>  </th> <th> ID </th> <th> pval </th> <th> logFC </th>  </tr>
<tr> <td align="right"> 1 </td> <td> HMDB00634 </td> <td> 5.077866e-34 </td> <td> -3.17008 </td> </tr>
<tr> <td align="right"> 2 </td> <td> HMDB00620 </td> <td> 5.077866e-34 </td> <td> -3.17008 </td> </tr>
<tr> <td align="right"> 3 </td> <td> HMDB02092 </td> <td> 5.077866e-34 </td> <td> -3.17008 </td> </tr>
<tr> <td align="right"> 4 </td> <td> HMDB00749 </td> <td> 5.077866e-34 </td> <td> -3.17008 </td> </tr>
<tr> <td align="right"> 5 </td> <td> HMDB10720 </td> <td> 2.273393e-32 </td> <td> -2.77394 </td> </tr>
<tr> <td align="right"> 6 </td> <td> HMDB03407 </td> <td> 2.273393e-32 </td> <td> -2.77394 </td> </tr>
</tbody></table></div>
</div>
</div>
</div>
</div>
</div>

Click "Make network" button and wait a few seconds till the network is constructed.
When it happens you'll be automatically scrolled down to a next panel.
You can scroll back up and check network summary or download it in XGMML:

<div class="row">
<div class="span3">
<div class="example-panel bottom-buffer">
<h3>Network summary</h3>
<div class="shiny-html-output shiny-bound-output"><ul>
<li> number of nodes :  1541 </li>
<li> number of edges :  1488 </li>
</ul>
</div>
<span class="btn shiny-download-link shiny-bound-output">Download XGMML</span><br>
</div>
</div>
</div>
<br>

Click "Find module" button to find a module in the network. You can preview it
on the right panel:

![Example of a module](img/module.png)


Click "Download XGMML" button in the "Module summary section" to download the
module in XGMML format for further analysis in Cytoscape:

<div class="row">
<div class="span3">
<div class="example-panel bottom-buffer">
<h3>Module summary</h3>
<div class="shiny-html-output shiny-bound-output"><ul>
<li> number of nodes :  34 </li>
<li> number of edges :  36 </li>
</ul>
</div>
<span class="btn shiny-download-link shiny-bound-output">Download XGMML</span>
</div>
</div>
</div>

### Details

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

The last column is optional, but we recommend to provide it if possible as it's
used for colors in graph visualizing. Any other columns will be copied to a
network as node or edge attributes.  We recommend to exclude genes with low
expression prior uploading. Example data for genes and metabolites can be
downloaded [here](/publications/supp_materials/GAM/).

After files are uploaded, a file summary is displayed. Verify that the files
were parsed correctly.

![Summary of the example data](img/data_summary.png)

If gene DE data was uploaded, the gene DE data is converted to reaction DE data.
This is done by considering for a reaction all the genes that code any enzyme
taking part in this reaction. P-values of these genes are compared. Gene with
a minimal p-value is selected and its p-value is assigned as the reaction
p-value.  All reactions without p-values are discarded as having no expressed
enzymes. We recommend to exclude genes with low expression prior uploading.

Reactions in the network can be interpreted as either edges or nodes.

**Reactions as edges** If reactions are interpreted as edges, if a pair of
metabolites has multiple reaction connections between them, only the reaction
with the minimal p-value is kept. All the remaining reactions and metabolites
that take part in these reactions are combined into a network as edges and
nodes respectively. To only consider cross-connections between substrates and
products in reactions that make up a KEGG reaction pair (RPAIR) as opposed to
all cross-connections, select the Use RPAIRs option (see [KEGG
REACTION](http://www.genome.jp/kegg/reaction/) for details).

Then you should select an interpretation of reactions in the network: either as
edges or as nodes.  For interpretation of reactions as edges the network is
constructed in a following way. For each pair of metabolites with multiple
reaction connections between them only one with a minimal p-value is kept.
Then all the remaining reactions and metabolites that take part in these
reactions are combined into a network as edges and nodes respectively. A
possible option is to consider not all cross-connections between substrates and
products for a reaction, but only those that make up a KEGG reaction pair
(RPAIR, see [KEGG REACTION](http://www.genome.jp/kegg/reaction/) for
details).

If reaction are interpreted as nodes then both metabolites and reactions are
added as nodes to the network. Edges are added between a metabolite and a
reaction if the metabolite takes part in this reaction.  Then, as an option,
groups of reactions having at least one common metabolite and the same most
significant gene can be collapsed into single nodes.

Generally, all this options lead to similar results. We recommend to use
default option values as this make a network simpler and analysis faster.

#### Finding a module

After making a network you can find a connected module that contains the most
significantly changed genes and reactions.  Internally this is done by firstly
scoring nodes and edges based on their p-values such that positive scores
correspond to significant p-values and negative scores correspond to
insignificant changes. Then the problem of finding a connected subgraph with
maximum summary weight (maximum-weight connected subgraph — MWCS — problem) is
solved.

This site supports three solvers:

* [Heinz solver](http://www.mi.fu-berlin.de/w/LiSA/Heinz) is a
  default solver, as it works for both interpretation of reactions.
* Using MWCS solver can provide more interactive experience because when
  reaching time limit it outputs the best solution found.
* Heuristic search is solver from R-package BioNet.

As the last two solvers don't support edge scores they are not available when
reactions are interpreted as edges.

By changing FDR values and scores for absent metabolites and reactions you can
control size of the module. We recommend to start from default values and then
gradually change them depending on the results.

There are couple of post-processing steps that are available. If in the network
reactions are edges and RPAIRs are used then trans- connections can be added.
When reactions are nodes following operations could be applied:

* Adding all metabolites that are part of any reaction in the module.
* Adding all reactions that aren't present in the module but are directly
  connected to two metabolites from the module.
* Removing all metabolites with no data and only one connection in the module.
* Replacing reaction nodes with edges for reactions that connect only two
  metabolites in the module and these two metabolites are on different sides of
  the reaction.

You can download the module in a XGMML format and, for example, explore it in
Cytoscape.  You can download GAM's VizMap style for Cytoscape
[here](/publications/supp_materials/GAM/GAM_VizMap.xml).
