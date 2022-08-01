# snRNA-Seq pipeline
Functionalized pipeline to perform snRNA-Seq analyses based in Seurat scRNA-Seq pipeline.


## The pipeline include several scripts
- Creation the RNAseq seurat object
- QC
- cell filtering based in features/UMIC/perc Mitochondrial DNA
- Differential expression analysis (DEA) using both MAST and DESeq2 and producing as output the annotation of both genes showing more than 0.25 log2FC difference between conditions and ll the expressed genes even those without any log2FC difference without conditions to be able to visualize per cluster all DEGs values even the ones that are not DEG in a cluster to avoid having to use false 0 values.
- Differential functional analyses using topGO to get top 300 GO geneSets altered per condition; and using pathFindR to get the Go geneset, REACTOME, KEGG pathways altered between condition to further perform downstream analyses
- trajectory analyses using monocle3
- Gene correlation analyses of a vector of genes of interest against another specific gene per cluster using as independant observations the counts of expression per each cell per cluster i) when the gene counts is not cero at least in one cell of the cluster; ii) when the gene counts in all the cells consider for each pairwise comparison of genes per cluster is not 0 fro any fo the paired genes. In the plots generated and excel results you can see the number of cells totally considered per cluster and the number of cells with no 0 values really used to calculate the final ii) correlation.

Last update June  2022.
