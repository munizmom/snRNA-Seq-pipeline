###################################################################################
# Functional analyses
# a) GO enrichment based on TopGO 
###################################################################################

#libraries needed:
library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("topGO");
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");library("org.Hs.eg.db");
library("stringr");library("tidyr");
options(future.globals.maxSize = 100000 * 1024^2) #Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:
# install org.Hs.eg.db from Bioconductor if not already installed 



#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd <- "/home/omicData/scRNASeq/"
setwd(wd);

f_results <- "/results/";
f_Rdata <- "RData/";
f_tables <- "tables/";
f_plots <- "plots/";
f_markers_p_Cluster<- "perCluster/"
f_dfa_topGo <- "dfa_topGO/"
f_pathfindR <-"pathfindR/"
#########################################################  ##  #####################################################
dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_pathfindR), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_Rdata,f_list), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_Rdata,f_list), showWarnings = F);

#########################################################  ##  #####################################################
########################
# A) GO enrichment
########################

#scrnaseq data seurat object
load(file = paste0(wd,f_results,f_Rdata,f_allClusters,"scRNAseq_Seurat_object.RData"));#scRNASeq_seurat

Idents(scRNASeq_seurat) <- "RNA_snn_res.0.05"


clusterList_4_GO_per_cluster <- list();
for (i in 1:length(annotaDf$cluster)){
	clusterList_4_GO_per_cluster[[i]] <- subset(scRNASeq_seurat, idents = paste0(i-1));
	names(clusterList_4_GO_per_cluster)[i] <- paste0("cluster_",i-1);
	i <- i+1;
};

annotaDf <- data.frame(cluster=unique(scRNASeq_seurat$RNA_snn_res.0.05), cellType=unique(scRNASeq_seurat$seurat_cellType))



GOenrichment_perCluster.function <- function(annotaDf,scRNASeq_seurat,nameSubClusters,conditionName,analysis,ontology,nameFolderOutput,nameOutput){
	#enrichment using topGO:
	# https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
	# https://ucdavis-bioinformatics-training.github.io/2021-March-Single-Cell-RNA-Seq-Analysis/data_analysis/scRNA_Workshop-PART6_fixed


	if (analysis=="perCluster"){
		Idents(scRNASeq_seurat) <- "RNA_snn_res.0.05"

			df <- scRNASeq_seurat;
			expr <- as.matrix(GetAssayData(df));
			All_exprsList <- expr;
			bad <- which(rowSums(expr) == 0);
			exprF <- expr[-bad,];
			exprsList <- exprF;
			resList <- df;
			sum_exprsList <- rowSums(expr);
	
			# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
			n.gt.0 <- apply(exprF, 1, function(x)length(which(x > 0)));
			expressed.genes <- rownames(exprF)[which(n.gt.0/ncol(exprF) >= 0.5)];
			all.genes <- rownames(exprF);
			allGenesList <- all.genes;
			expressedF <-  as.data.frame(n.gt.0[which(names(n.gt.0 ) %in% expressed.genes)]);
			expressedF$Gene <- rownames(expressedF);
			expressedGenesList <- expressedF;
			

			# define geneList as 1 if gene is in expressed.genes, 0 otherwise
			geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
			names(geneList) <- all.genes
			
			geneListRes <- geneList;
			names(geneListRes) <- paste0(nameSubClusters,"_",i-1);
			
			# Create topGOdata object
			GOdata <- new("topGOdata",
					ontology = ontology, # use biological process ontology
					allGenes = geneList,
					geneSelectionFun = function(x)(x == 1),
			              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol");
			GOdataList <-GOdata;
			# Test for enrichment using Fisher's Exact Test
			resultFisher  <- runTest(GOdata, algorithm = "elim", statistic = "fisher");
			resultFisherList  <-resultFisher;
			
			GenTable <- GenTable(GOdata, Fisher = resultFisher, topNodes = 5000, numChar = 750);
			GenTable$BH <- p.adjust(GenTable$Fisher, method="BH");
			GenTableList  <-GenTable;
			goIDList <- GenTable[, "GO.ID"];
			

	} else {

			df <- subset(x=scRNASeq_seurat, subset= APOE_LA == conditionName );
			expr <- as.matrix(GetAssayData(df));
			All_exprsList <- expr;
			bad <- which(rowSums(expr) == 0)
			exprF <- expr[-bad,]
			exprsList <- exprF;
			resList <- df;
			sum_exprsList <- rowSums(expr);
				
			# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
			n.gt.0 <- apply(exprF, 1, function(x)length(which(x > 0)));
			expressed.genes <- rownames(exprF)[which(n.gt.0/ncol(exprF) >= 0.5)];
			all.genes <- rownames(exprF)
			allGenesList <- all.genes;
			expressedF <-  as.data.frame(n.gt.0[which(names(n.gt.0 ) %in% expressed.genes)]);
			expressedF$Gene <- rownames(expressedF);
			expressedGenesList <- expressedF;
	

			# define geneList as 1 if gene is in expressed.genes, 0 otherwise
			geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
			names(geneList) <- all.genes
			
			geneListRes <- geneList;
			
			# Create topGOdata object
			GOdata <- new("topGOdata",
					ontology = "BP", # use biological process ontology
					allGenes = geneList,
					geneSelectionFun = function(x)(x == 1),
			              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol");
			GOdataList <-GOdata;
			# Test for enrichment using Fisher's Exact Test
			resultFisher  <- runTest(GOdata, algorithm = "elim", statistic = "fisher");
			resultFisherList  <-resultFisher;
			
			GenTable <- GenTable(GOdata, Fisher = resultFisher, topNodes = 1000, numChar = 350);
			GenTable$BH <- p.adjust(GenTable$Fisher, method="BH");
			GenTableList  <-GenTable;

			goIDList <- GenTable[, "GO.ID"];

	};

	#save(GenTableList, file=paste0(wd,f_results,f_dfa_topGo,f_Rdata,f_list,"list_GO_",ontology,"_",nameSubClusters,"_",conditionName,".RData"));
	#save(resList, geneListRes, All_exprsList, exprsList, sum_exprsList, GOdataList,resultFisherList,goIDList, file=paste0(wd,f_results,f_dfa_topGo,f_Rdata,f_list,"GoSupp_",ontology,"_",nameSubClusters,"_",conditionName,".RData"));
	#save( expressedGenesList,All_exprsList, file=paste0(wd,f_results,f_dfa_topGo,f_Rdata,f_list,"ExpressedGenes_",ontology,"_",nameSubClusters,"_",conditionName,".RData"));


	excel_DFA_Results.function <- function(dfa,ontology,conditionName,nameFolderOutput,annotaDf,nameOutput){

		dfF <- dfa[[1]];
		dfF$cluster <- gsub("^.*_","",names(dfa)[1]);
		dfF$FunctionalTerm <-ontology;
		dfF <- left_join(dfF,annotaDf,by="cluster");

		assign(paste0("dfF"),dfF, envir=parent.frame());
		save(dfF, file=paste0(wd,f_dfa_topGo,f_Rdata,nameOutput,"_df_GO_",ontology,"_",conditionName,".RData"));
		openxlsx::write.xlsx(dfF, file = paste0(wd,f_results,f_dfa_topGo,f_tables,nameOutput,"_GO_",ontology,"_",conditionName, "Df.xlsx"),rowNames =FALSE, colNames =TRUE);
	};

	excel_DFA_Results.function(GenTableList,ontology,conditionName,nameFolderOutput,annotaDf,nameOutput);
	#assign(paste0("GO_",ontology,"_",conditionName),dfF,.GlobalEnv);
};		

GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond1_cluster","cond1","perCondition","BP",f_markers_p_Cluster,"allCellTypes");
GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond2_cluster","cond2","perCondition","BP",f_markers_p_Cluster,"allCellTypes");
GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond1_cluster","cond1","perCondition","CC",f_markers_p_Cluster,"allCellTypes");
GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond2_cluster","cond2","perCondition","CC",f_markers_p_Cluster,"allCellTypes");
GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond1_cluster","cond1","perCondition","MF",f_markers_p_Cluster,"allCellTypes");
GOenrichment_perCluster.function(annotaDf,scRNASeq_seurat,"cond2_cluster","cond2","perCondition","MF",f_markers_p_Cluster,"allCellTypes");




Pathway_enrichment_cond1_vs_cond2.function <- function(annotaDf,scRNASeq_seurat,conditionName,ontology,nameOutput,topNodes	){


	Idents(scRNASeq_seurat) <- "RNA_snn_res.0.05"

	df <- scRNASeq_seurat;
	expr <- as.matrix(GetAssayData(df));
	All_exprsList <- expr;
	bad <- which(rowSums(expr) == 0);
	exprF <- expr[-bad,];
	exprsList <- exprF;
	resList <- df;
	sum_exprsList <- rowSums(expr);

	# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
	n.gt.0 <- apply(exprF, 1, function(x)length(which(x > 0)));
	expressed.genes <- rownames(exprF)[which(n.gt.0/ncol(exprF) >= 0.5)];
	all.genes <- rownames(exprF);
	allGenesList <- all.genes;
	expressedF <-  as.data.frame(n.gt.0[which(names(n.gt.0 ) %in% expressed.genes)]);
	expressedF$Gene <- rownames(expressedF);
	expressedGenesList <- expressedF;


	# define geneList as 1 if gene is in expressed.genes, 0 otherwise
	geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
	names(geneList) <- all.genes

	geneListRes <- geneList;

	# Create topGOdata object
	GOdata <- new("topGOdata",
			ontology = ontology, # use biological process ontology
			allGenes = geneList,
			geneSelectionFun = function(x)(x == 1),
	              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol");
	GOdataList <-GOdata;
	# Test for enrichment using Fisher's Exact Test
	resultFisher  <- runTest(GOdata, algorithm = "elim", statistic = "fisher");
	resultFisherList  <-resultFisher;

	GenTable <- GenTable(GOdata, Fisher = resultFisher, topNodes = topNodes, numChar = 750);
	GenTable$BH <- p.adjust(GenTable$Fisher, method="BH");
	GenTableList  <-GenTable;
	goIDList <- GenTable[, "GO.ID"];

	allGO <- genesInTerm(GOdata); ##get all GOs and their genes from the topgo result.
	annotGenesin_dataset <- lapply(allGO,function(x) x[x %in% expressedGenesList[,2]] ); 
	tableGO.genes <- annotGenesin_dataset[goIDList]; 
	genes.df = as.data.frame(do.call(rbind, tableGO.genes)); ## Transform it to a data frame.
	tableGO.genes <- tableGO.genes[lengths(tableGO.genes) != 0];

	for (i in 1:length(tableGO.genes)){
	    if (i==1){
	 		df <- as.data.frame(tableGO.genes[[1]]);

			df$Term <- names(tableGO.genes)[1];		 		
	 		colnames(df) <- c("Gene","GO.ID")
	    } else{
	  	    df.tmp <- as.data.frame(tableGO.genes[[i]])
	  		df.tmp$Term <- names(tableGO.genes)[i];		 		
	  	    colnames(df.tmp) <- c("Gene","GO.ID");
			df <- bind_rows(df,df.tmp)
	  	}
	  i <-i+1;
	  }
	goID_name_conversion <- unique(GenTable[,c("GO.ID","Term")])
	assign(paste0("df"),df,.GlobalEnv);
	annotPathsGenes <-left_join(df,goID_name_conversion,by="GO.ID")
	annotPathsGenes$PathId <- paste0(annotPathsGenes$GO.ID, "___",annotPathsGenes$Term)
	annotPathsGenes_stats <-    unique(aggregate( Gene ~ PathId  , data = annotPathsGenes, paste, collapse = "// "));
	annotPathsGenes_stats<- separate(annotPathsGenes_stats, PathId, c("GO.ID" ,"Term"  ),sep="___");
	annotPathsGenes_stats$PathId <- paste0(annotPathsGenes_stats$GO.ID, " ",annotPathsGenes_stats$Term)
	annotPathsGenes_stats$nbGenes_per_Path <- str_count(annotPathsGenes_stats$Gene, "//") +1;
	colnames(annotPathsGenes_stats)[3] <- "All_genes_altered_per_Pathway"

	annotPathsGenes_statsF <- left_join(GenTable, annotPathsGenes_stats, by=c("GO.ID","Term"))
	openxlsx::write.xlsx(annotPathsGenes_statsF, file = paste0(wd,f_results,f_dfa_topGo,f_tables,nameOutput,"_GO_",ontology,"_",conditionName, "Df.xlsx"),rowNames =FALSE, colNames =TRUE);

	annotPathsGenes$PathId <- paste0(annotPathsGenes$GO.ID, " ",annotPathsGenes$Term)
	annotPathsGenesF <- left_join(GenTable, annotPathsGenes, by=c("GO.ID","Term"))
	annotPathsGenesF$ontology <- ontology;
	annotPathsGenesF$data <- nameOutput;
	openxlsx::write.xlsx(annotPathsGenesF, file = paste0(wd,f_results,f_dfa_topGo,f_tables,nameOutput,"_GO_",ontology,"_",conditionName, "4Network_Df.xlsx"),rowNames =FALSE, colNames =TRUE);
	save(annotPathsGenesF, file=paste0(wd,f_results,f_dfa_topGo,f_Rdata,f_list,"list_GO_",ontology,"_",nameOutput,"_",conditionName,".RData"));
	assign(paste0("GO_",ontology,"_",conditionName,"_",nameOutput),annotPathsGenesF,.GlobalEnv);
};
Pathway_enrichment_cond1_vs_cond2.function(annotaDf,scRNASeq_seurat,"cond1_vs_2","BP","allCellTypes",5000);
Pathway_enrichment_cond1_vs_cond2.function(annotaDf,scRNASeq_seurat,"cond1_vs_2","MF","allCellTypes",4000);
Pathway_enrichment_cond1_vs_cond2.function(annotaDf,scRNASeq_seurat,"cond1_vs_2","CC","allCellTypes",1900);
