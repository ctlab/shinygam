<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
    <script src="shared/jquery.js" type="text/javascript"></script>
    <script src="shared/shiny.js" type="text/javascript"></script>
    <link rel="stylesheet" type="text/css" href="shared/shiny.css"/>
    <link rel="stylesheet" type="text/css" href="shared/slider/css/jquery.slider.min.css"/>
    <script src="shared/slider/js/jquery.slider.min.js"></script>


    <link rel="stylesheet" type="text/css" href="shared/bootstrap/css/bootstrap.min.css"/>
    <script src="shared/bootstrap/js/bootstrap.min.js"></script>

    <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
    <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>
    <script src="gam.js" type="text/javascript"></script>
    <link rel="stylesheet" type="text/css" href="gam.css"/>
    <link rel="stylesheet" type="text/css" href="shared/bootstrap/css/bootstrap-responsive.min.css"/>
    <title>Shiny GAM</title>
  </head>
  <body>
<div>
</div>
      <div id="networkParameters" class="js-output">
          <script>
              network = {
                  available: false,
                  hasReactionsAsEdges: false,
                  hasReactionsAsNodes: false,
                  hasGenes: false,
                  usesRpairs: false
              };
          </script>
      </div>
      <div id="moduleParameters" class="js-output">
          <script>
              module = {
                  available: false
              };
          </script>
      </div>
    <div class="container">
      <div class="row">
        <div class="span12" style="padding: 10px 0px;">
          <h1>Shiny GAM</h1>
        </div>
      </div>
      <div class="row">
          <div class="span12">
              <div class="tabbable">
                  <ul class="nav nav-tabs" id="tab-panel-main">
                      <!--<li class="active"> <a href="#tab-home" data-toggle="tab">Home</a> </li>-->
                      <li class="active"> <a href="#tab-work" data-toggle="tab">Work</a> </li>
                      <li> <a href="#tab-help" data-toggle="tab">Help</a> </li>
                      <li> <a href="#tab-about" data-toggle="tab">About</a> </li>
                  </ul>
                  <div class="tab-content">
                      <!--
                      <div class="tab-pane active" id="tab-home">
                      </div>
                      -->
                      <div class="tab-pane active" id="tab-work">
                          <div class="row">
                              <div class="span3">
                                  <form class="well">
                                      <p class="form-header">Construct a network</p>
                                      <label class="control-label" for="network">
                                          Select an organism: 
                                          <a id="select-organism-tooltip" class="help-tooltip" href="#" 
                                              title="Select an organism that you want to analyse. If there is no one that you want contact the author (see About tab).">
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <select id="network" class="selectpicker">
                                          <option value="Mouse musculus" selected="selected">Mouse musculus</option>
                                          <option value="Homo sapiens">Homo sapiens</option>
                                      </select>
                                      <label>
                                          File with gene DE data: 
                                          <a id="upload-gene-data-tooltip" class="help-tooltip" href="#" 
                                              title='File with gene differential expression data shoud be tab-separated with header and contain columns "ID" and "pval".'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <input id="geneDE" type="file" accept="text/tab-separated-values"/>
                                      <label>
                                          File with metabolic DE data: 
                                          <a id="upload-met-data-tooltip" class="help-tooltip" href="#" 
                                              title='File with metabolic differential expression data shoud be tab-separated with header and contain columns "ID" and "pval".'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <input id="metDE" type="file" accept="text/tab-separated-values"/>
                                      <label class="control-label" for="reactionsAs">
                                          How to interpret reactions:
                                          <a id="reactionAs-tooltip" class="help-tooltip" href="#" 
                                              title='Reactions can be either edges connecting substrates with products or nodes connected to them. Reactions as edges makes analysis more flux-centric.'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <select id="reactionsAs">
                                          <option value="edges" selected="selected">As edges</option>
                                          <option value="nodes">As nodes</option>
                                      </select>
                                      <div data-display-if="input.reactionsAs == 'nodes'">
                                          <label class="checkbox" for="collapseReactions">
                                              <input id="collapseReactions" type="checkbox" checked="checked"/>
                                              <span>Collapse reactions
                                                  <a id="collapseReactions-tooltip" class="help-tooltip" href="#" 
                                                      title='Collapse reaction nodes with the same enzyme.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                      </div>
                                      <div data-display-if="input.reactionsAs != 'nodes'">
                                          <label class="checkbox" for="useRpairs">
                                              <input id="useRpairs" type="checkbox" checked="checked"/>
                                              <span>Use RPAIRs
                                                  <a id="useRpairs-tooltip" class="help-tooltip" href="#" 
                                                      title='Connect only those substrate-product pairs that make up a KEGG RPAIR.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                      </div>
                                      <button id="preprocess" type="button" class="btn action-button">Make network</button>
                                  </form>
                              </div>
                              <div class="span9">
                                  <div class="row">
                                      <div class="span4">
                                          <h3>Genomic data</h3>
                                          <div id="geneDESummary" class="shiny-html-output"></div>
                                          <div id="geneDETable" class="shiny-html-output"></div>
                                      </div>
                                      <div class="span4">
                                          <h3>Metabolic data</h3>
                                          <div id="metDESummary" class="shiny-html-output"></div>
                                          <div id="metDETable" class="shiny-html-output"></div>
                                      </div>
                                  </div>
                                  <div class="bottom-buffer">
                                      <h3>Network summary</h3>
                                      <div id="networkSummary" class="shiny-html-output"></div>
                                      <div data-display-if="network.available">
                                          <a id="downloadNetwork" class="btn shiny-download-link" href="" target="_blank">Download XGMML</a><br>
                                      </div>
                                  </div>
                              </div>
                          </div>
                          <div id="module-panel" class="row top-buffer">
                              <div class="span3">
                                  <form class="well">
                                      <p class="form-header">Discover modules</p>
                                      <div data-display-if="network.hasGenes">
                                          <label for="geneFDR">FDR for reactions
                                              <a id="geneFDR-tooltip" class="help-tooltip" href="#" 
                                                  title='Increasing/decreasing this value makes adding reactions to a module harder/easier.'>
                                                  <i class="icon-question-sign"></i>
                                              </a>
                                          </label>
                                          <input id="geneFDR" type="number" value="1e-06" min="1e-100" max="1" step="1e-200"/>
                                      </div>

                                      <label for="metFDR">FDR for metabolites
                                          <a id="metFDR-tooltip" class="help-tooltip" href="#" 
                                              title='Increasing/decreasing this value makes adding metabolites to a module harder/easier.'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <input id="metFDR" type="number" value="1e-06" min="1e-100" max="1" step="1e-200"/>
                                      <!--<input id="metMinusLogFDR" name="metMinusLogFDR" class="minusLogFDR" type="slider" value="5" data-skin="plastic"/>-->

                                      <label for="absentMetScore">Score for absent metabolites
                                          <a id="absentMetScore-tooltip" class="help-tooltip" href="#" 
                                              title='Increasing/decreasing this value makes adding metabolites with no data to a module harder/easier.'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                      <input id="absentMetScore" type="number" value="-20"/>

                                      <div data-display-if="!network.hasGenes">
                                          <label for="absentRxnScore">Score for absent reactions
                                              <a id="absentRxnScore-tooltip" class="help-tooltip" href="#" 
                                                  title='Increasing/decreasing this value makes adding reactions with no data to a module harder/easier (usable only when there is no gene expression data).'>
                                                  <i class="icon-question-sign"></i>
                                              </a>
                                          </label>
                                          <input id="absentRxnScore" type="number" value="-10"/>
                                      </div>

                                      <label class="control-label" for="solver">Select a solver:
                                          <a id="selectSolver-tooltip" class="help-tooltip" href="#" 
                                              title='Multiple MWCS solvers can be used to find modules. Check the help page for details.'>
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>

                                      <select id="solver">
                                          <option value="heinz" selected="selected">Heinz</option>
                                          <option value="mwcs">MWCS (under development)</option>
                                          <option value="fastHeinz">Heuristic search</option>
                                      </select>

                                      <div data-display-if="$('#solver').val() == 'mwcs'">
                                          <label for="mwcsTimeLimit">MWCS solver time limit (seconds):
                                          <a id="selectSolver-tooltip" class="help-tooltip" href="#" 
                                              title="It's an approximate maximum time MWCS solver is allowed to run. After hitting the time limit it will print the best solurion found so far.">
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                          <input id="mwcsTimeLimit" type="number" value="10" min="0" max="120"/>
                                      </div>

                                      <div data-display-if="$('#solver').val() == 'heinz'">
                                          <label for="heinzTimeLimit">Heinz time limit (seconds):
                                          <a id="selectSolver-tooltip" class="help-tooltip" href="#" 
                                              title="It's an approximate maximum time heinz solver is allowed to run. After hitting the time limit it will just stop searching outputing no module.">
                                              <i class="icon-question-sign"></i>
                                          </a>
                                      </label>
                                          <input id="heinzTimeLimit" type="number" value="60" min="0" max="240"/>
                                      </div>
                                      <button id="find" type="button" class="btn action-button">Find module</button>
                                      <div data-display-if="network.hasReactionsAsEdges && network.usesRpairs">
                                          <label class="checkbox" for="addTransPairs">
                                              <input id="addTransPairs" type="checkbox"/>
                                              <span>Add trans- edges
                                                  <a id="addTransPairs-tooltip" class="help-tooltip" href="#" 
                                                      title='Add trans- RPAIR edges.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                      </div>
                                      <div data-display-if="network.hasReactionsAsNodes">
                                          <label class="checkbox" for="addMetabolitesForReactions">
                                              <input id="addMetabolitesForReactions" type="checkbox"/>
                                              <span>Add all metabolites for reactions
                                                  <a id="addMetabolitesForReactions-tooltip" class="help-tooltip" href="#" 
                                                      title='Add all metabolites that are part of any reaction in the module.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                          <label class="checkbox" for="addInterconnections">
                                              <input id="addInterconnections" type="checkbox"/>
                                              <span>Add interconnections
                                                  <a id="addInterconnections-tooltip" class="help-tooltip" href="#" 
                                                      title="Addg all reactions that aren't present in the module but are directly connected to two metabolites from the module.">
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                          <label class="checkbox" for="removeHangingNodes">
                                              <input id="removeHangingNodes" type="checkbox"/>
                                              <span>Remove hanging nodes
                                                  <a id="removeHangingNodes-tooltip" class="help-tooltip" href="#" 
                                                      title='Remove all metabolites with no data and only one connection in the module.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                          <label class="checkbox" for="removeSimpleReactions">
                                              <input id="removeSimpleReactions" type="checkbox"/>
                                              <span>Remove simple reactions
                                                  <a id="removeSimpleReactions-tooltip" class="help-tooltip" href="#" 
                                                      title='Replacing reaction nodes with edges for reactions that connect only two metabolites in the module and these two metabolites are on different sides of the reaction.'>
                                                      <i class="icon-question-sign"></i>
                                                  </a>
                                              </span>
                                          </label>
                                      </div>
                                  </form>
                                  <h3>Module summary</h3>
                                  <div id="moduleSummary" class="shiny-html-output"></div>
                                  <div data-display-if="module.available">
                                      <a id="downloadModule" class="btn shiny-download-link" href="" target="_blank">Download XGMML</a>
                                  </div>
                              </div>
                              <div class="span9">
                                  <div id="module" class="graph-output"></div>
                              </div>
                          </div>
                      </div>
                      <div class="tab-pane" id="tab-help">
                          <div class="row">
                              <div class="span9 offset3">
                                  m4_include(help.xhtml)
                              </div>
                          </div>
                      </div>
                      <div class="tab-pane" id="tab-about">
                          <div class="row">
                              <div class="span9 offset3">
                                  <p> By Alexey Sergushichev (<a href="mailto:asergushichev@wustl.edu">asergushichev@wustl.edu</a>)</p>
                                  <p> This site must not be used for commercial purposes.<p>
                                  <p> You can find source code of this site <a href="https://bitbucket.org/assaron/gam">here</a>.</p>
                                  <div id="GAMVersion" class="shiny-html-output"></div>
                                  <p> Links <p>
                                  <p> Powered by Shiny. </p>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
          </div>
      </div>
    </div>
    <div id="showModulePanel" class="js-output"></div>
  </body>
</html>

m4_dnl vim: set syntax=html:
