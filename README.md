## Melianea RNA-Seq data analysis with AskoR 

**AskoR** is a pipeline for the analysis of gene expression data, available at [Github/AskoR](https://github.com/askomics/askoR).

* [User Guide](https://github.com/asusete/askoR/wiki/Pipeline-askoR:-User-Guide)
* [Parameters Table](https://github.com/asusete/askoR/wiki/Pipeline-askoR:-Parameters-Table)

### Melinaea data : 

You'll find the Melinaea data in the **data/input** folder. 

  - Samples.csv : description of all the samples with corresponding species, tissue and stage
  - Contrasts.txt : tabular file including all the contrasts performed by AskoR
  - MmEP.gos : Ids of Gene Ontology termes assigned to the contigs of the Melinea marsaeus transcriptome 



### Download and Run 
You can clone this repository, then run with Rscript or within Rstudio or your favorite R environment (including the AskoR dependencies) the  R script file (see _scriptR/AskoR_analysis_script.R_):  

### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.
