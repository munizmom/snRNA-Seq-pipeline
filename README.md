# snRNA-Seq pipeline
Functionalized pipeline to perform snRNA-Seq analyses based in Seurat scRNA-Seq pipeline.


## The pipeline include several scripts
- Creation the RNAseq seurat object
- QC
- cell filtering based in features/UMIC/perc Mitochondrial DNA
- Differential expression analysis (DEA) using both MAST and DESeq2 and producing as output the annotation of both genes showing more than 0.25 log2FC difference between conditions and ll the expressed genes even those without any log2FC difference without conditions to be able to visualize per cluster all DEGs values even the ones that are not DEG in a cluster to avoid having to use false 0 values.
- Differential functional analyses using topGO to get top 300 GO geneSets altered per condition; and using pathFindR to get the Go geneset, REACTOME, KEGG pathways altered between conditio ¡n to further perform downstream analyses
- trajectory analyses using monocle3


Last update June  2022.
