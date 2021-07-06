## Melinaea marsaeus RNA-Seq data analysis with AskoR 

**AskoR** is a pipeline for the analysis of gene expression data, available at [Github/AskoR](https://github.com/askomics/askoR).

* [User Guide](https://github.com/asusete/askoR/wiki/Pipeline-askoR:-User-Guide)
* [Parameters Table](https://github.com/asusete/askoR/wiki/Pipeline-askoR:-Parameters-Table)

### Melinaea data : 

You'll find the Melinaea data in the **data/input** folder. 

  - Samples.csv : description of all the samples with corresponding species, tissue and stage
  - Contrasts.txt : tabular file including all the contrasts performed by AskoR
  - MmEP.gos : Ids of Gene Ontology termes assigned to the contigs of the Melinea marsaeus transcriptome 



### Download and Run 
You can clone this repository, then run with Rscript or within Rstudio or your favorite R environment (including the AskoR dependencies) the  R script file (see _ScriptR/AskoR_analysis_script.R_).

### Results

Results will be exported to the _data/DE_ directory :
  - DE_CPM_NormCounts.txt : TMM normalized CPM for all genes and samples
  - DE_CPM_NormMeancounts.txt : TMM normalized CPM for all genes and conditions 
  - DE_summary.txt : Differential expression tests for all genes and contratst (0 not DE, 1 upregulated, -1 down regulated) with FDR  < 0.05 
  -  _data/DE/DataExplore_ : directory including the figures  of the filtering steps (barplot) and correlation analyses (Correlogram, heatmaps, clustering and MDS) 
  -  _data/DE/DEanalysis_/AskoTables : directory including all the results tables by contrasts (genes and tests values (fold changee,  pvalues, FDR)  
  -  _data/DE/DEanalysis_/AskoTables : figures of the DE tests by contrasts  (volcano plots and heatmaps of the highly DEGs)
  -  _data/DE/DEanalysis_/GOenrichment/OnContrasts : tables and Bubbles graphs of the GO enrichments
  -  _data/DE/DEanalysis_/UpsetGraphs/Global_upset/ : UpSet plot of the intersections of all the contrasts
  -  _data/DE/DEanalysis_/UpsetGraphs/Subset_upset/ : Upset plot of the intersections of the Mmp_1vsMmr_1, Mmp_2vsMmr_2 and Mmp_aavsMmr_aa contrasts 
### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.
