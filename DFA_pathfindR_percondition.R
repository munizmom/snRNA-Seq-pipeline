library("pathfindR");library("readxl");library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("MAST");
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");library("org.Hs.eg.db");
library("stringr");library("tidyr");library("stringr");
options(future.globals.maxSize = 100000 * 1024^2) #Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:


Sys.setenv(RSTUDIO_PANDOC="/home/software/anaconda3/envs/r4-base/bin/pandoc")
# Sys.getenv("RSTUDIO_PANDOC")
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
f_dfa_topGo <- "dfa_topGO/"
f_list <- "res_list_dfs/"
f_pathfindR <-"pathfindR/"
f_customPaths <-"enrichment_MTDpaths/"
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

dir.create(file.path(wd, f_results, f_customPaths), showWarnings = F);
dir.create(file.path(wd, f_results, f_customPaths,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_customPaths,f_tables), showWarnings = F);

#########################################################  ##  #####################################################

dataInput<- as.data.frame(read_excel(paste0("/home/omicData/scRNASeq/rData/input/DEGs/pval005_log2fc_1_minus1/MAST_output_DEA.xlsx"),sheet =1, col_names =T)); 
degs_input <- dataInput[grep("Significant & ",dataInput$Group),];	


#degs_input <-degs_input[!(degs_input$cluster==0),] #oligo c0 only one gene pathway enrichmen gives error ofc


pathfindR_extraction.function <- function(dataDf,conditionName,nameOutput){
	for ( i in 1:length(unique(dataDf$cluster))){
		resPaths <- unique(dataDf[which(dataDf$cluster==unique(dataDf$cluster)[i]),c("Gene","avg_log2FC","p_val_adj")]);
		resPathsOriginal <- unique(dataDf[which(dataDf$cluster==unique(dataDf$cluster)[i]),]);
		colnames(resPaths) <-c("Gene.symbol","log2FC","adj.P.Val");
		resPaths<-resPaths[!is.na(resPaths$`adj.P.Val` ),]
		resPaths_kegg <- run_pathfindR(resPaths, output_dir = paste0(wd,f_results,f_pathfindR,nameOutput,"_",conditionName,"_KEGG_c",unique(dataDf$cluster)[i]),
			min_gset_size = 5, max_gset_size = 500, plot_enrichment_chart =FALSE,visualize_enriched_terms=FALSE)
		resPaths_GOs <- run_pathfindR(resPaths, output_dir = paste0(wd,f_results,f_pathfindR,nameOutput,"_",conditionName,"_GO_c",unique(dataDf$cluster)[i]),gene_sets = "GO-All",min_gset_size = 5, max_gset_size = 500, plot_enrichment_chart =FALSE,visualize_enriched_terms=FALSE)
		resPaths_reactome <- run_pathfindR(resPaths, output_dir = paste0(wd,f_results,f_pathfindR,nameOutput,"_",conditionName,"_REACTOME_c",unique(dataDf$cluster)[i]),gene_sets = "Reactome",min_gset_size = 5, max_gset_size = 500, plot_enrichment_chart =FALSE,visualize_enriched_terms=FALSE)

		resPathsF <- bind_rows(resPaths_kegg,resPaths_GOs,resPaths_reactome)
		openxlsx::write.xlsx(resPathsF, file = paste0(wd,f_results,f_pathfindR,f_tables,nameOutput,"_",conditionName,"_c",unique(dataDf$cluster)[i], "_4Network_Df.xlsx"),rowNames =FALSE, colNames =TRUE);
		save(resPathsF, file=paste0(wd,f_results,f_pathfindR,f_Rdata,f_list,"paths_",nameOutput,"_",conditionName,"_c",unique(dataDf$cluster)[i],".RData"));
		assign(paste0("paths_",conditionName,"_",nameOutput,"_c",unique(dataDf$cluster)[i]),resPathsF,.GlobalEnv);
		
		for (j in 1:nrow(resPathsF)){
			if (j==1){
				up <-  as.data.frame(str_split(resPathsF[j,8],", "));
				up$ID <- resPathsF[j,"ID"];
				up$GeneRegulation <- "Upregulated";
				down <-  as.data.frame(str_split(resPathsF[j,9],", "));
				down$ID <- resPathsF[j,"ID"];
				down$GeneRegulation <- "Downregulated";
				colnames(up)[1] <-"Gene";
				colnames(down)[1] <-"Gene";
				df <- bind_rows(up,down);
				dfF <- df;
			} else{
				up.tmp <-  as.data.frame(str_split(resPathsF[j,8],", "));
				up.tmp$ID <- resPathsF[j,"ID"];
				up.tmp$GeneRegulation <- "Upregulated";
				down.tmp <-  as.data.frame(str_split(resPathsF[j,9],", "));
				down.tmp$ID <- resPathsF[j,"ID"];
				down.tmp$GeneRegulation <- "Downregulated";
				colnames(up.tmp)[1] <-"Gene";
				colnames(down.tmp)[1] <-"Gene";
				df.tmp <- bind_rows(up.tmp,down.tmp);
				dfF <- bind_rows(dfF,df.tmp);

			};
			j<- j+1

		};
		network <- left_join(dfF,resPathsF,by="ID");
		network <- left_join(network,resPathsOriginal,by="Gene");

		assign(paste0("network_",conditionName,"_",nameOutput,"_c",unique(dataDf$cluster)[i]),network,.GlobalEnv);
		openxlsx::write.xlsx(network, file = paste0(wd,f_results,f_pathfindR,f_tables,nameOutput,"_",conditionName,"_c",unique(dataDf$cluster)[i], "_4Network_Df.xlsx"),rowNames =FALSE, colNames =TRUE);
		save(resPathsF,network, file=paste0(wd,f_results,f_pathfindR,f_Rdata,f_list,"paths_",nameOutput,"_",conditionName,"_c",unique(dataDf$cluster)[i],".RData"));
		i <- i+1;
	};


};
#run 1
conditionName <-"cond1_vs_cond2"
nameOutput <-"DEGs"
pathfindR_extraction.function(degs_input,"cond1_vs_cond2","DEGs");

resPaths <-mget(ls(pattern="paths_"));

for (ii in 1:(length(resPaths))){
	if(ii==1){
	  df <- resPaths[[ii]];
	  pathsDf=df;
	} else{
	  df.tmp <- resPaths[[ii]];
	  pathsDf<- bind_rows(pathsDf,df.tmp);

	};
	ii <- ii+1;
};
 rm(ii);
NetPaths <-mget(ls(pattern="network_"));

for (ii in 1:(length(NetPaths))){
	if(ii==1){
	  df <- NetPaths[[ii]];
	  NetpathsDf=df;
	} else{
	  df.tmp <- NetPaths[[ii]];
	  NetpathsDf<- bind_rows(NetpathsDf,df.tmp);

	};
	ii <- ii+1;
};


save(resPaths,NetpathsDf, file=paste0(wd,f_results,f_pathfindR,f_Rdata,"paths_all_clusters_DEGs_",conditionName,"_",nameOutput,".RData"));
save(resPaths,NetpathsDf, file=paste0(wd,f_results,f_pathfindR,f_Rdata,"paths_all_clusters_DEGs_",conditionName,"_",nameOutput,".RData"));
openxlsx::write.xlsx(pathsDf, file = paste0(wd,f_results,f_pathfindR,f_tables,conditionName, "_dfa_pathfindR_stats_allClusterss.xlsx"),rowNames =FALSE, colNames =TRUE);
openxlsx::write.xlsx(NetpathsDf, file = paste0(wd,f_results,f_pathfindR,f_tables, conditionName,"_dfa_pathfindR_4Net_allClusters.xlsx"),rowNames =FALSE, colNames =TRUE);





