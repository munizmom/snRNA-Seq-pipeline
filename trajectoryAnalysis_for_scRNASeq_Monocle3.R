###############################################################################################################
# script to perform trajectory analyses using monocle3 in seurat processed cluster objects :
######################
# 1) load all the clusters object called scRNASeq_seurat
#
# 2) load or produce each subcluster : inhNeu.clusters,excNeu.clusters, astrocytes.clusters, microglia.clusters,...
#
# 3) perform the trajectory analyses. then remove that object to not collapse R memory
# The annotation on each trajectory analysis is stored again in scRNASeq_seurat and used consecutively
# NOTE: at the end of the trajectory analysis the experiment.clusters object fully annotated is now stored in
# paste0(wd,f_results,f_Rdata,f_allClusters,f_monocle) not in paste0(wd,f_results,f_Rdata,f_allClusters anymore.
# But I am maintining both to be consistant with the steps.
# Script done and maintained by Mar Muniz Moreno. 
#
# to run in the cluster instructions
# seurat,signac and seurat wrapper use tools with a diff version not compatible with monocle
#thus install 2 conda environments, one specifically for monocle but allow int he cluster the use of both R libraries
###############################################################################################################

#####################
#load packages
#####################
library(Signac);library(Seurat);library(SeuratWrappers);library(monocle3);library(Matrix);
library(ggplot2);library(patchwork);library(dplyr);library(colorspace);library(ggplot2);
library(RColorBrewer);library(gtools);
#####################
#set seeed and wd
#####################

set.seed(22); # need to set it to get always the same random results and plots

wd <- "/home/omicData/scRNASeq/"
setwd(wd);

#####################
#set output directories
#####################
f_results <- "/results/";
f_Rdata <- "RData/";
f_qc <- "QC/";
f_violin <-"violinPLot/";
f_quantiles <-"quantiles/";
f_cellcycle<- "cellcycle/";
f_pca<- "pca/";
f_preIntegration <-"PreIntegration/";
f_after_integration <-"After_integration/";
f_ridgeplot <- "RidgePlot/";
f_foundClusters <- "FoundClusters/";
f_tables <- "tables/";
f_plots <- "plots/";
f_markers <- "markers/";
f_featurePlot <- "featurePlots/"
f_apoe <- "apoe/"
 f_umap_tsne <- "tsne_vs_umap/"
f_finalMarkers <- "finalMarkers/"
f_tree <- "tree/"
f_annotCluster <- "annotated_Cluster_cellType/"
f_markers_p_Cluster_cond <- "perCluster_per_condition/"
f_markers_p_Cluster_cond_Sex <- "perCluster_per_condition_PER_Sex/"
f_markers_p_Cluster <- "perCluster/"
f_LAgenes <-"LAregion_44_46kbs/"
f_ADgenes<- "ADgenes/"
f_dfa_topGo <- "dfa_topGO/"
f_subclusters <-"extracted_subClusters/"
f_allClusters <-"allClusters/"
f_dea_MAST <- "DEA_MAST/"
f_monocle <- "trajectoryAnalysis_monocle3/"
f_phase <- "phase/"
f_list <- "res_list_dfs/"


#####################
#create directories
#####################
dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata,f_dfa_topGo), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata,f_subclusters), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata,f_allClusters), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata,f_allClusters,f_monocle), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata,f_dea_MAST), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_violin), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_ridgeplot), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_quantiles), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_quantiles), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_cellcycle), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_after_integration), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_preIntegration,f_pca), showWarnings = F);
dir.create(file.path(wd, f_results,f_qc ,f_after_integration,f_pca), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata ,f_after_integration), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata ,f_preIntegration), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_umap_tsne), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_finalMarkers), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_tree), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_annotCluster), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_featurePlot,f_apoe), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_featurePlot,f_markers), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_featurePlot,f_markers), showWarnings = F);
dir.create(file.path(wd,f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond,f_LAgenes), showWarnings = F);
dir.create(file.path(wd,f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex,f_LAgenes), showWarnings = F);
dir.create(file.path(wd,f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond,f_phase), showWarnings = F);
dir.create(file.path(wd,f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex,f_phase), showWarnings = F);

dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_ADgenes), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_ADgenes,f_markers_p_Cluster_cond_Sex), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_ADgenes,f_markers_p_Cluster), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond, f_dfa_topGo), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex, f_dfa_topGo), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond, f_dfa_topGo,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond, f_dfa_topGo,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex, f_dfa_topGo,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex, f_dfa_topGo,f_plots), showWarnings = F);

dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond, f_dfa_topGo,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex, f_dfa_topGo,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond, f_dfa_topGo,f_Rdata,f_list), showWarnings = F);
dir.create(file.path(wd, f_results,f_foundClusters,f_markers,f_markers_p_Cluster_cond_Sex, f_dfa_topGo,f_Rdata,f_list), showWarnings = F);

############################################################################################



####################################
### load the pseudotime monocle trajectory function
####################################

pseudoTime_trajectory_analysis_monocle3.function <- function(allclusters,subcluster,nameOutput,Idents_name_rootCells,clusterIdentityMetadata,f_output){

	cellsNb_Subcluster <- length(Idents(subcluster)); #112781
	subcluster.cds <- as.cell_data_set(subcluster);
	subcluster.cds <- cluster_cells(cds = subcluster.cds, reduction_method = "UMAP",k = 20);
	subcluster.cds <- learn_graph(subcluster.cds, use_partition = TRUE);

	#establish 2 pseudotimes analysis, one based on previous knowledge, another based on statistical distances
	#results to be stored in:
	subcluster.cds1 <- subcluster.cds		
	subcluster.cds2 <- subcluster.cds		

	#1) seting us up the root cells based in our knowledge
	root_cells <- subset(x = subcluster, idents = c(Idents_name_rootCells), invert = FALSE);
	root_cells_Name <- names(Idents(root_cells))
	cellsNb_rootCluster <- length(Idents(root_cells)); #9020 astrocytes c1
	subcluster.cds1 <- order_cells(subcluster.cds1, reduction_method = "UMAP", root_cells = root_cells_Name);

	#2 setting up the root cells based on statistical distances, method provided in monocle3 doc
	#  helper function to identify the root principal points:
	# The function below does so by first grouping the cells according to which trajectory graph 
	# node they are nearest to. Then, it calculates what fraction of the cells at each node come 
	# from the earliest time point. Then it picks the node that is most heavily occupied by early 
	# cells and returns that as the root.
	#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
	# I modified a bit to produce the tables of results to keep track of which cluster/node is 
	# automatically selected as root as it was not aparent or easy to track.
	#normally is the same that the one you select manually by looking at the seurat clusters and dendrogram
	get_earliest_principal_node <- function(cds,clusterIdentityMetadata,nameOutput){
		cell_ids <- colData(cds)[, clusterIdentityMetadata];

		closest_vertex <-
		cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
		closest_vertex <- as.matrix(closest_vertex[colnames(cds), ]);
		root_pr_nodes <-
		igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
		(which.max(table(closest_vertex[cell_ids,]))))];
		assign(paste0("root_pr_nodes_",nameOutput),root_pr_nodes, envir=parent.frame());

		#saving the info of clusters node assignation
		clustersInfo <- as.data.frame(table(cell_ids))[-which(-as.data.frame(table(cell_ids))[,2]==0),]
		nodesInfo <- as.data.frame(table(closest_vertex[cell_ids,]));
		colnames(nodesInfo)[1] <- "NodeName";
		nodesInfo[,1] <- paste0("Y_",nodesInfo[,1]);
		root_Node <- data.frame(NodeName=root_pr_nodes,RootSelected="Yes");
		nodesClusters_info <- full_join(clustersInfo,nodesInfo,by="Freq");
		nodesClusters_info <- full_join(nodesClusters_info,root_Node,by="NodeName");
		colnames(nodesClusters_info)[2] <- "Nb of cells"
		write.csv(nodesClusters_info, file=paste0(wd,f_results,f_monocle,f_output, nameOutput,"_cluster_root_nodes_Info.csv"),row.names = FALSE);
		#Passing the programatically selected root node to order_cells() via the root_pr_nodeargument yields:

	root_pr_nodes;
	};
	subcluster.cds2 <- order_cells(subcluster.cds2, root_pr_nodes=get_earliest_principal_node(subcluster.cds2,clusterIdentityMetadata,nameOutput));
	#print(root_pr_nodes_AstrocytesTrajectory);

	# plot trajectories colored by pseudotime, 
	# NOTE: cannot on the cluster as the browser is not allowed by ssh.
	###############################################
	# plotTrajectory_pseudo <- plot_cells(
	#  cds = subcluster.cds,
	#  color_cells_by = "pseudotime",
	#  show_trajectory_graph = TRUE,
	#       graph_label_size=1.5);
	#
	# plotTrajectory_pseudo2 <-plot_cells(subcluster.cds2,
    #       color_cells_by = "pseudotime",
    #       label_cell_groups=FALSE,
    #       show_trajectory_graph = TRUE,
    #       graph_label_size=1.5);
    #
	#Extract the pseudotime values and add to the Seurat object
	allclusters <- AddMetaData(
		object = allclusters,
		metadata = subcluster.cds1@principal_graph_aux@listData$UMAP$pseudotime,
		col.name = paste0(gsub("Trajectory","",nameOutput), "_v_rootManual"));


	allclusters <- AddMetaData(
		object = allclusters,
		metadata = subcluster.cds2@principal_graph_aux@listData$UMAP$pseudotime,
		col.name = paste0(gsub("Trajectory","",nameOutput), "_v_rootAuto"));
	scRNASeq_seurat <-allclusters
	save(scRNASeq_seurat, file=paste0(wd,f_results,f_Rdata,f_allClusters,f_monocle,"experiment.clusters_monocle.RData"));
	save(scRNASeq_seurat, file=paste0(wd,f_results,f_Rdata,f_allClusters,f_monocle,nameOutput,"experiment.clusters_monocle.RData"));
	save(scRNASeq_seurat, file=paste0(wd,f_results,f_monocle,f_output,nameOutput,"_experiment.clusters.RData"));
	assign(paste0("scRNASeq_seurat"),allclusters,.GlobalEnv);

};



#####################
# load all the clusters object
#####################

Idents(scRNASeq_seurat) <- "integrated_snn_res.0.5"

excNeu.clusters <- subset(x = scRNASeq_seurat, idents = c("2", "3","9", "10","14", "15","17", "24","28", "34","35", "39", "41"), invert = FALSE)
astrocytes.clusters <-subset(x = scRNASeq_seurat, idents = c("1", "8", "7","4"), invert = FALSE)
inhNeu.clusters <-subset(x = scRNASeq_seurat, idents = c("12", "13","11", "37","19", "25","38", "40"), invert = FALSE)
microg.clusters <-subset(x = scRNASeq_seurat, idents = c("5", "42"), invert = FALSE)
oligosOPC.clusters <-subset(x = scRNASeq_seurat, idents = c("0", "6","16","18","20","21","23","26","27","29","31","33","36"), invert = FALSE)

excNeu.clusters <- RenameIdents(
  object = excNeu.clusters,
 '2' = 'Exc_neurons_L3_L5_FMN1_RORB_TOX_FMN1_CBLN2_NRG1_c2',
  '3' = 'Exc_neurons_L2_L3_PDZD2_CBLN2_LAMA2_PARD3_c3',
 '9' = 'Exc_neurons_L4_L6_RORB2_FOXP2_TOX_FAM19A1_LRRK1_SEMA3E_c9',
  '10' = 'Exc_neurons_L1_L3_CUX2_PDZD2_CBLN2_SLC38A11_MOXD1_LAMP5_c10',
  '14' = 'Exc_neurons_L4_PLCH1_L5_TOX_L6_PDZRN4_TSHZ2_c14',
  '15' = 'Exc_neurons_L2_L5_CUX2_RORB_TOX_FOXP2_COL5A2_c15',
  '17' = 'Exc_neurons_L6_TOX_THEMIS_FSTL5_c17',
  '24' = 'Exc_neurons_L6_PDZRN4_HS3ST4_CDH12_KIAA1217_CDH9_c24',
  '28' = 'Exc_neurons_L5_L6_FOXP2_TOX_PDZRN4_HS3ST4_KIAA1217_c28',
  '34' = 'Exc_neurons_L4_L6_RORB_TOX_PDZRN4_TSHZ2_CDH12_c34',
  '35' = 'Exc_neurons_L6_THEMIS_PDZRN4_HS3ST4_CPNE4_c35',
  '39' = 'Exc_neurons_L6_PDZRN4_CPNE4_HTR2C_c39',
  '41' = 'Exc_neurons_L5_L6_FOXP2_PDZRN4_SORCS1_ASIC2_c41')

excNeu.clusters
astrocytes.clusters <- RenameIdents(
  object = astrocytes.clusters,
  '1' = 'Astrocytes_c1',
 '4' = 'Astrocytes_c4',
  '7' = 'Astrocytes_c7',
  '8' = 'Astrocytes_reactive_c8');

microglia.clusters  <- RenameIdents(
  object = microg.clusters,
  '5' = 'Microglia_PTPRC_DOCK8_low_c5',
  '42' = 'Microglia_PTPRC_DOCK8_high_c42');


inhNeu.clusters <- RenameIdents(
  object = inhNeu.clusters,
  '11' = 'Inh_neurons_LHX6_SNCG_c11',
  '12' = 'Inh_neurons_LHX6_PVALB_c12',
  '13' = 'Inh_neurons_LHX6_SST_c13',
  '19' = 'Inh_neurons_ADARB2_VIP_c19',
  '25' = 'Inh_neurons_ADARB2_LAMP5_c25',
  '37' = 'Inh_neurons_TCF4_GRIA4_DDP10_KZN_c37',
  '38' = 'Inh_neurons_ADARB2_LAMP5_c38',
  '40' = 'Inh_neurons_ADARB2_MYO16_KCNT2_RELN_c40');

oligosOPC.clusters  <- RenameIdents(
  object = oligosOPC.clusters,
  '0' = 'Oligo_c0',
  '6' = 'OPCs_c6',
  '16' = 'Oligo_c16',
  '18' = 'Oligo_c18',
  '20' = 'Oligo_c20',
  '21' = 'Oligo_c21',
  '23' = 'Oligo_c23',
  '26' = 'Oligo_c26',
  '27' = 'Oligo_c27',
  '29' = 'Oligo_c29',
  '31' = 'Oligo_c31',
  '33' = 'Oligo_c33',
  '36' = 'Oligo_c36');

save(astrocytes.clusters,file=paste0(wd,f_results,f_Rdata,f_subclusters,"astrocytes.cluster.RData"))
save(inhNeu.clusters,file=paste0(wd,f_results,f_Rdata,f_subclusters,"inhNeu.cluster.RData"))
save(excNeu.clusters,file=paste0(wd,f_results,f_Rdata,f_subclusters,"excNeu.cluster.RData"))
save(microglia.clusters,file=paste0(wd,f_results,f_Rdata,f_subclusters,"microg.cluster.RData"))
save(oligosOPC.clusters,file=paste0(wd,f_results,f_Rdata,f_subclusters,"oligosOPC.cluster.RData"))

Idents(scRNASeq_seurat) <- "cellType_res05"



f_cond_1 <- "cond1/";
f_cond_2 <- "cond2/";
#####################
# perform each analysis consecutively
#####################
#1

#per LA trajectory
#########################################################################

nameOutput <-"AstrocytesTrajectory";
dir.create(file.path(wd,f_results,f_monocle,nameOutput,"/"), showWarnings = F);
Idents_name_rootCells<- "Astrocytes_c1";
clusterIdentityMetadata<- "cellType_res05";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,astrocytes.clusters,nameOutput,
	Idents_name_rootCells,clusterIdentityMetadata,paste0(nameOutput,"/"));


Idents(astrocytes.clusters) <- "perCluster_per_condition"
astrocytes.clusters_cond1 <- subset(x = astrocytes.clusters, idents = c("cond1"), invert = FALSE);
Idents(astrocytes.clusters_cond1) <- "cellType_res05"
astrocytes.clusters_cond2 <- subset(x = astrocytes.clusters, idents = c("cond2"), invert = FALSE);
Idents(astrocytes.clusters_cond2) <- "cellType_res05"

dir.create(file.path(wd,f_results,f_monocle,"AstrocytesTrajectory/"), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"AstrocytesTrajectory/",f_cond_1), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"AstrocytesTrajectory/",f_cond_2), showWarnings = F);
Idents_name_rootCells<- "Astrocytes_c1";
clusterIdentityMetadata<- "cellType_res05";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,astrocytes.clusters_cond2,"astrocytes_AF",
	Idents_name_rootCells,clusterIdentityMetadata, paste0("AstrocytesTrajectory/",f_cond_2));

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,astrocytes.clusters_cond1,"astrocytes_EU",
	Idents_name_rootCells,clusterIdentityMetadata,paste0("AstrocytesTrajectory/",f_cond_1));

rm(astrocytes.clusters);
rm(astrocytes.clusters_cond1);
rm(astrocytes.clusters_cond2);

#2

#per LA trajectory
#########################################################################

nameOutput <-"ExcNeuronTrajectory";
dir.create(file.path(wd,f_results,f_monocle,nameOutput,"/"), showWarnings = F);
clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Exc_neurons_L1_L3_CUX2_PDZD2_CBLN2_SLC38A11_MOXD1_LAMP5_c10";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,excNeu.clusters,nameOutput,
	Idents_name_rootCells,clusterIdentityMetadata,paste0(nameOutput,"/"));


Idents(excNeu.clusters) <- "perCluster_per_condition"
excNeu.clusters_cond1 <- subset(x = excNeu.clusters, idents = c("cond1"), invert = FALSE);
Idents(excNeu.clusters_cond1) <- "cellType_res05"
excNeu.clusters_cond2 <- subset(x = excNeu.clusters, idents = c("cond2"), invert = FALSE);
Idents(excNeu.clusters_cond2) <- "cellType_res05"

dir.create(file.path(wd,f_results,f_monocle,"ExcNeuronTrajectory/"), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"ExcNeuronTrajectory/",f_cond_1), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"ExcNeuronTrajectory/",f_cond_2), showWarnings = F);


clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Exc_neurons_L1_L3_CUX2_PDZD2_CBLN2_SLC38A11_MOXD1_LAMP5_c10";

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,excNeu.clusters_cond2,"Exc_neurons_AF",
	Idents_name_rootCells,clusterIdentityMetadata, paste0("ExcNeuronTrajectory/",f_cond_2));

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,excNeu.clusters_cond1,"Exc_neurons_EU",
	Idents_name_rootCells,clusterIdentityMetadata,paste0("ExcNeuronTrajectory/",f_cond_1));

rm(excNeu.clusters);
rm(excNeu.clusters_cond1);
rm(excNeu.clusters_cond2);


#3
nameOutput <-"InhNeuronTrajectory";
dir.create(file.path(wd,f_results,f_monocle,nameOutput,"/"), showWarnings = F);
clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Inh_neurons_LHX6_SNCG_c11";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,inhNeu.clusters,nameOutput,
	Idents_name_rootCells,clusterIdentityMetadata,paste0(nameOutput,"/"));


#per LA trajectory
#########################################################################
Idents(inhNeu.clusters) <- "perCluster_per_condition"
inhNeu.clusters_cond1 <- subset(x = inhNeu.clusters, idents = c("cond1"), invert = FALSE);
Idents(inhNeu.clusters_cond1) <- "cellType_res05"
inhNeu.clusters_cond2 <- subset(x = inhNeu.clusters, idents = c("cond2"), invert = FALSE);
Idents(inhNeu.clusters_cond2) <- "cellType_res05"

dir.create(file.path(wd,f_results,f_monocle,"InhNeuronTrajectory/"), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"InhNeuronTrajectory/",f_cond_1), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"InhNeuronTrajectory/",f_cond_2), showWarnings = F);

clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Inh_neurons_LHX6_SNCG_c11";

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,inhNeu.clusters_cond1,"inhNeu_EU",
	Idents_name_rootCells,clusterIdentityMetadata, paste0("InhNeuronTrajectory/",f_cond_1));

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,inhNeu.clusters_cond2,"inhNeu_AF",
	Idents_name_rootCells,clusterIdentityMetadata,paste0("InhNeuronTrajectory/",f_cond_2));

rm(inhNeu.clusters);
rm(inhNeu.clusters_cond1);
rm(inhNeu.clusters_cond2);


#4
nameOutput <-"MicrogliaTrajectory";
dir.create(file.path(wd,f_results,f_monocle,nameOutput,"/"), showWarnings = F);
clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Microglia_PTPRC_DOCK8_low_c5";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,microglia.clusters,nameOutput,
	Idents_name_rootCells,clusterIdentityMetadata,paste0(nameOutput,"/"));


#per LA trajectory
#########################################################################
Idents(microglia.clusters) <- "perCluster_per_condition"
microglia.clusters_cond1 <- subset(x = microglia.clusters, idents = c("cond1"), invert = FALSE);
Idents(microglia.clusters_cond1) <- "cellType_res05"
microglia.clusters_cond2 <- subset(x = microglia.clusters, idents = c("cond2"), invert = FALSE);
Idents(microglia.clusters_cond2) <- "cellType_res05"

dir.create(file.path(wd,f_results,f_monocle,"MicrogliaTrajectory/"), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"MicrogliaTrajectory/",f_cond_1), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"MicrogliaTrajectory/",f_cond_2), showWarnings = F);

clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "Microglia_PTPRC_DOCK8_low_c5";

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,microglia.clusters_cond1,"microglia_EU",
	Idents_name_rootCells,clusterIdentityMetadata, paste0("MicrogliaTrajectory/",f_cond_1));

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,microglia.clusters_cond2,"microglia_AF",
	Idents_name_rootCells,clusterIdentityMetadata,paste0("MicrogliaTrajectory/",f_cond_2));

rm(microglia.clusters);
rm(microglia.clusters_cond1);
rm(microglia.clusters_cond2);


#5
nameOutput <-"OligoTrajectory";
dir.create(file.path(wd,f_results,f_monocle,nameOutput,"/"), showWarnings = F);
clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "OPCs_c6";
pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,oligosOPC.clusters,nameOutput,
	Idents_name_rootCells,clusterIdentityMetadata,paste0(nameOutput,"/"));



#per LA trajectory
#########################################################################
Idents(oligosOPC.clusters) <- "perCluster_per_condition"
oligosOPC.cluster_cond1 <- subset(x = oligosOPC.clusters, idents = c("cond1"), invert = FALSE);
Idents(oligosOPC.cluster_cond1) <- "cellType_res05"
oligosOPC.cluster_cond2 <- subset(x = oligosOPC.clusters, idents = c("cond2"), invert = FALSE);
Idents(oligosOPC.cluster_cond2) <- "cellType_res05"

dir.create(file.path(wd,f_results,f_monocle,"OligoTrajectory/"), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"OligoTrajectory/",f_cond_1), showWarnings = F);
dir.create(file.path(wd,f_results,f_monocle,"OligoTrajectory/",f_cond_2), showWarnings = F);

clusterIdentityMetadata<- "cellType_res05";
Idents_name_rootCells <- "OPCs_c6";

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,oligosOPC.cluster_cond1,"oligosOPC_EU",
	Idents_name_rootCells,clusterIdentityMetadata, paste0("OligoTrajectory/",f_cond_1));

pseudoTime_trajectory_analysis_monocle3.function(scRNASeq_seurat,oligosOPC.cluster_cond2,"oligosOPC_AF",
	Idents_name_rootCells,clusterIdentityMetadata,paste0("OligoTrajectory/",f_cond_2));

rm(oligosOPC.clusters);
rm(oligosOPC.cluster_cond1);
rm(oligosOPC.cluster_cond2);



trajectoryPlots.function <- function(allClusters,columns2Plot,nameOutput,height, width,f_output,date){
	p <- FeaturePlot(scRNASeq_seurat,columns2Plot , pt.size = 0.1) & scale_color_viridis_c(); 

	png(file=paste0(wd,f_results,f_monocle,f_output,"Traj_",nameOutput,"_",date,".png"), height = height, width = width, units = "in", res = 1200);
	print(p);
	dev.off();
};

date <- "062222"
trajectoryPlots.function(scRNASeq_seurat,c("Astrocytes_v_rootManual", "Astrocytes_v_rootAuto"),"astrocytes",6,12,"AstrocytesTrajectory/",date);
trajectoryPlots.function(scRNASeq_seurat,c("ExcNeuron_v_rootManual" ,"ExcNeuron_v_rootAuto")  ,"ExcNeurons",6,12,"ExcNeuronTrajectory/",date);
trajectoryPlots.function(scRNASeq_seurat,c("InhNeuron_v_rootManual" ,"InhNeuron_v_rootAuto")  ,"InhNeurons",6,12,"InhNeuronTrajectory/",date);
trajectoryPlots.function(scRNASeq_seurat,c("Microglia_v_rootManual" ,"Microglia_v_rootAuto")  ,"Microglia",6,12,"MicrogliaTrajectory/",date);
trajectoryPlots.function(scRNASeq_seurat,c("Oligo_v_rootManual" ,"Oligo_v_rootAuto")  ,"Oligo",6,12,"OligoTrajectory/",date);
trajectoryPlots.function(scRNASeq_seurat,c("astrocytes_cond2_v_rootManual", "astrocytes_cond2_v_rootAuto"),"astrocytes_AF",6,12,paste0("AstrocytesTrajectory/",f_cond_2),date);
trajectoryPlots.function(scRNASeq_seurat,c("astrocytes_cond1_v_rootManual", "astrocytes_cond1_v_rootAuto"),"astrocytes_EU",6,12,paste0("AstrocytesTrajectory/",f_cond_1),date);

trajectoryPlots.function(scRNASeq_seurat,c("Exc_neurons_cond2_v_rootManual", "Exc_neurons_cond2_v_rootAuto"),"Exc_neurons_AF",6,12,paste0("ExcNeuronTrajectory/",f_cond_2),date);
trajectoryPlots.function(scRNASeq_seurat,c("Exc_neurons_cond1_v_rootManual", "Exc_neurons_cond1_v_rootAuto"),"Exc_neurons_EU",6,12,paste0("ExcNeuronTrajectory/",f_cond_1),date);


trajectoryPlots.function(scRNASeq_seurat,c("inhNeu_cond2_v_rootManual", "inhNeu_cond2_v_rootAuto"),"inhNeu_AF",6,12,paste0("InhNeuronTrajectory/",f_cond_2),date);
trajectoryPlots.function(scRNASeq_seurat,c("inhNeu_cond1_v_rootManual", "inhNeu_cond1_v_rootAuto"),"inhNeu_EU",6,12,paste0("InhNeuronTrajectory/",f_cond_1),date);


trajectoryPlots.function(scRNASeq_seurat,c("microglia_cond2_v_rootManual", "microglia_cond2_v_rootAuto"),"microglia_AF",6,12,paste0("MicrogliaTrajectory/",f_cond_2),date);
trajectoryPlots.function(scRNASeq_seurat,c("microglia_cond1_v_rootManual", "microglia_cond1_v_rootAuto"),"microglia_EU",6,12,paste0("MicrogliaTrajectory/",f_cond_1),date);


trajectoryPlots.function(scRNASeq_seurat,c("oligosOPC_cond2_v_rootManual", "oligosOPC_cond2_v_rootAuto"),"oligosOPC_AF",6,12,paste0("OligoTrajectory/",f_cond_2),date);
trajectoryPlots.function(scRNASeq_seurat,c("oligosOPC_cond1_v_rootManual", "oligosOPC_cond1_v_rootAuto"),"oligosOPC_EU",6,12,paste0("OligoTrajectory/",f_cond_1),date);




dimPlotUmap_combinedData_perCluster_cellType.function <- function(combined.log,resolution,width,height,nameOutput,date,metadataCondition){


		plot <- DimPlot(combined.log, reduction = "umap", label = FALSE, group.by="cellType_res05", split.by = metadataCondition,  raster = FALSE,ncol=2) + labs(title = paste0( resolution," Resolution")); #& theme(legend.position="bottom")

		png(file=paste0(wd,f_results,f_foundClusters,f_plots,f_annotCluster,"combined_",nameOutput,"_" ,resolution,"_UMAP_cellType_",metadataCondition,date,".png"), height = height, width = width, units = "in", res = 1200);
		print(plot);
		dev.off();



};
date<- "080622"
dimPlotUmap_combinedData_perCluster_cellType.function(scRNASeq_seurat,05,35,10,"integrated_snn_res.0.5",date,NULL);


#save.image(file=paste0(wd,f_results,f_Rdata,f_allClusters,f_monocle,"session_trajectoriesMonocle_",date,".RData"));



