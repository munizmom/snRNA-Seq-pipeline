

#########################################################################################
#########################################################################################
#correlation gene expression APOE with all the MIT altered genes per cell cluster
#########################################################################################
#########################################################################################
#########################################################################################
# script following the Mit_dysfunction_analysis.R script
# developed and maintained by Mar Muniz
#
#output folder
f_corr <-"correlation/thresholdDEG_1_minus1/";

################################
# scRNAseq data with cells coming per condition Mutant vs control for example
################################

####################################################################################################
#1) input data
####################################################################################################
#1 preparing input of genes of interest to test
# extract all genes normalized counts expression of all cells per cluster per condition
genesInterestcond_FC12<- unique(c(mit_percond_FC12[grep("Significant & ", mit_percond_FC12$Group),"Gene"],gene2Corr));
genesInterestcond<- unique(c(mit_percond[grep("Significant & ", mit_percond$Group),"Gene"],gene2Corr));

####################################################################################################
#2 cell data coming from the seurat onject of scRNASeq per cluster
####################################################################################################
load(file = paste0(wd,f_results,f_Rdata,f_allClusters,"scRNAseq_Seurat_object.RData"));#scRNASeq_seurat

res_logCounts_percond <- list();
Idents(scRNASeq_seurat) <- "cluster.percond"
for (i in 1:length(unique(scRNASeq_seurat$cluster.percond))){
  seuratObject <- subset(scRNASeq_seurat, idents = paste0(unique(scRNASeq_seurat$cluster.percond)[i]));
  res_logCounts_percond[[i]] <- as.matrix(GetAssayData(seuratObject), slot = "data");
  names(res_logCounts_percond)[i] <- unique(scRNASeq_seurat$cluster.percond)[i];
  i <- i+1;
};
res_logCounts_percond_1 <- res_logCounts_percond[c(1:43)]  #clusters with cells from cond1
res_logCounts_percond_2 <- res_logCounts_percond[c(44:86)] #clusters with cells from cond2

save(input_mitoXplorer_hsa,genesInterestcond, genesInterestcond_FC12,res_logCounts_percond_1,res_logCounts_percond_2, file=paste0(wd,f_results, f_corr,"cond_data_4corrAnalysis.RData"))


####################################################################################################
# 3) functions to calculate and plot the correlation of any pair of genes on a dataframe using each
#  cell member of a cluster as independant measument of the correlation pair
####################################################################################################
####################################################################################################


MTDcorrelation.function <- function(input_mitoXplorer_hsa,countsList, genesInterest, date,nameOutput ,gene2Corr){
  #the correlation calculated here will include also cells where the gene has a value of 0 counts of expression if
  # that gene is found expressed in at least 1 cell of the cluster. Thus to get the corr values removing cells where 
  # any of the genes whose corr we want to assess is 0 after this function we should run MTDcorrelation_not0_value.function

  f_corr <-"correlation/thresholdDEG_1_minus1/";
  f_per_group <-"percond_per_sex/"
  dir.create(file.path(wd,f_results,f_corr), showWarnings = F);
  dir.create(file.path(wd,f_results,f_corr,nameOutput,"/"), showWarnings = F);
  dir.create(file.path(wd,f_results,f_corr,nameOutput,"/",f_per_group), showWarnings = F);

  MTDcategories_perGene <- unique(input_mitoXplorer_hsa[,c("Gene", "Mitochondrial_process","gene_function")]);
  resList_MTDcorr <- list();
  resList_MTDnot0_corr <- list();
  resList_inputCorr_all <- list();
  resList_inputCorr_not0 <- list();


  for ( i in 1:length(names(countsList))){
  nameData <- gsub("__","_",gsub(".*_c","c",names(countsList)[i]));    
    a <-countsList[[i]]
    b <- a[which(rownames(a) %in% genesInterest ),]; #row:genes, cols:cells
    b.not0 <- b[,-which(b==0)];  #analysing just cells where none of the interesting genes have a value =0.
    not0 <- apply(b, 1, function(x)length(which(x > 0)));
        
    b.expressed <- rownames(b)[which(not0/ncol(b) >= 0.25)]; #genes expressed in cells totally more thna a 25 %
    b.expressed <- b[which(rownames(b) %in% b.expressed ),]; #row:genes, cols:cells


    resList_inputCorr_all[[i]] <- b;
    resList_inputCorr_not0[[i]] <- b.not0;

    names(resList_inputCorr_all)[i] <- names(countsList)[i];
    names(resList_inputCorr_not0)[i] <- names(countsList)[i];


    if (is.null(dim(b.not0))==FALSE && dim(b.not0)[2]>4  ){
        c.not0 <-  rcorr(t(b.not0), type="spearman")$r;
        c.not0_pval <- rcorr(t(b.not0), type="spearman")$P;
        c.not0_ncells <- rcorr(t(b.not0), type="spearman")$n;
        d.not0 <- c.not0[which(rownames(c.not0)==gene2Corr),];
        d.not0_pval <- c.not0_pval[which(rownames(c.not0_pval)==gene2Corr),];
        d.not0_ncells <- c.not0_ncells[which(rownames(c.not0_ncells)==gene2Corr),];
       d.not0_fdr <- p.adjust(d.not0_pval, method = "fdr", n = length(d.not0_pval));
        e.not0 <- data.frame(Gene=names(d.not0), corr=d.not0,pval=d.not0_pval,
         FDR=d.not0_fdr ,ncellsTested=d.not0_ncells,name=names(countsList)[i],
          Cluster=gsub("__.*","",names(countsList)[i]), 
          LA=gsub(".*__","",gsub("__M","",gsub("__F","",names(countsList)[i]))) ,
          Sex=gsub(".*__","",names(countsList)[i]));
         e.not0 <- left_join(e.not0,MTDcategories_perGene,by="Gene");
        resList_MTDnot0_corr[[i]] <-e.not0;
        names(resList_MTDnot0_corr)[i] <- names(countsList)[i];

        openxlsx::write.xlsx(e.not0, file = paste0(wd,f_results,f_corr,nameOutput ,"/",f_per_group,nameData,"_MTDnot0Cells_corr_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
        colnames(e.not0)[c(2:5)] <- paste0(colnames(e.not0)[c(2:5)], "_",names(countsList)[i] )

    };

    if (is.null(dim(b))==FALSE && dim(b)[2]>4  ){
      c <-  rcorr(t(b), type="spearman")$r;
      c_pval <- rcorr(t(b), type="spearman")$P;
      c_ncells <- rcorr(t(b), type="spearman")$n;
      d <- c[which(rownames(c)==gene2Corr),];
      d_pval <- c_pval[which(rownames(c_pval)==gene2Corr),];
      d_ncells <- c_ncells[which(rownames(c_ncells)==gene2Corr),];
      d_fdr <- p.adjust(d_pval, method = "fdr", n = length(d_pval));
      e <- data.frame(Gene=names(d), corr=d,pval=d_pval,
        FDR=d_fdr ,ncellsTested=d_ncells,name=names(countsList)[i],
        Cluster=gsub("__.*","",names(countsList)[i]), 
        LA=gsub(".*__","",gsub("__M","",gsub("__F","",names(countsList)[i]))) ,
        Sex=gsub(".*__","",names(countsList)[i]));
       e <- left_join(e,MTDcategories_perGene,by="Gene");

    resList_MTDcorr[[i]] <-e;
    names(resList_MTDcorr)[i] <- names(countsList)[i];

  resList <- list(corr_allCells=e);
  #openxlsx::write.xlsx(resList, file = paste0(wd,f_results,f_corr,nameOutput ,"/",f_per_group,nameData,"_MTDcorr_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
};


  i <-i+1;
  print(i)
  };
 
  resList_inputCorr_not0 <- resList_inputCorr_not0[lapply(resList_inputCorr_not0,length)>0]

  assign(paste0(nameOutput,"_resList_MTDnot0_corr"),resList_MTDnot0_corr,.GlobalEnv);
  assign(paste0(nameOutput,"_List_inputCorr_not0"),resList_inputCorr_not0,.GlobalEnv);

  for (ii in 1:(length(resList_MTDcorr))){
    if(ii==1){
      df <- resList_MTDcorr[[ii]];
      dfF=df;
    } else{
      df.tmp <- resList_MTDcorr[[ii]];
      dfF<- bind_rows(dfF,df.tmp);

    };
    ii <- ii+1;
  };

 
  for (iii in 1:(length(resList_MTDnot0_corr))){
    if(iii==1){
      df2 <- resList_MTDnot0_corr[[iii]];
      dfF2=df2;
    } else{
      df2.tmp <- resList_MTDnot0_corr[[iii]];
      dfF2<- bind_rows(dfF2,df2.tmp);
    };
    iii <- iii+1;
  };
  assign(paste0(nameOutput,"_",gene2Corr"_MTDcorr_not0Cells_df"),dfF2,.GlobalEnv);

  resListF <- list(corr_allCells=dfF,corr_cellsNot0=dfF2);
  #resListF <- list(corr_allCells=dfF);
  openxlsx::write.xlsx(resListF, file = paste0(wd,f_results,f_corr,nameOutput ,"/",nameData,"_MTDcorr_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
 
};


MTDcorrelation_not0_value.function <- function(nameOutput,inputList,date,gene2Corr){
    #NOTE:  the correlation calculated here will include only cells where the gene has a 
    # value of counts of expression>0 should be run after MTDcorrelation.function
  
 
  #function needed
  ############################################# ->
  lineal_corrFunction.function <- function(data,df1,nameOutput,inputColumnsFC,namePlot,nameSavePlot,f_cluster,ncells4corrPlot,cellsTotal,gene2Corr){

    f_corr <- "correlation/thresholdDEG_1_minus1/";
    f_linearCorr <- "linearCorr/";
    f_data <- paste0(nameOutput,"/");
    dir.create(file.path(wd,f_results,f_corr), showWarnings = F);
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr), showWarnings = F);
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr,f_data), showWarnings = F);
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr,f_data,f_cluster), showWarnings = F);


    pdf(paste0(wd,f_results, f_corr, f_linearCorr,f_data,f_cluster, "Scatter_corrPlot_",nameSavePlot,".pdf")) #,  width=13
    par(oma = c(0, 0, 2, 0))

    pairs.panels(round(df1[,c(1:2)],2), 
                 method = "spearman", # correlation method
                 hist.col = "pink2",
                 density = TRUE,cor=TRUE,smooth=TRUE,  # show density plots
                 ellipses = FALSE,
                 cex.cor=0.6,cex=0.7,digits = 2, # show correlation ellipses
                 stars=TRUE,breaks=10) #show significance starts
    title(paste0("Spearman corr ",namePlot, " nb cells not 0 values : ",ncells4corrPlot, "(",cellsTotal, ")"), outer=TRUE)


   pairs.panels(round(data[,inputColumnsFC],2), 
                 method = "spearman", # correlation method
                 hist.col = "pink2",
                 density = TRUE,cor=TRUE,smooth=FALSE,  # show density plots
                 ellipses = FALSE,
                 cex.cor=0.6,cex=0.7,digits = 2, # show correlation ellipses
                 stars=TRUE,breaks=10) #show significance starts
    title(paste0("Spearman corr ",namePlot, " nb total cells: ",cellsTotal), outer=TRUE)

    pairs.panels(round(df1[,c(1:2)],2), 
                 method = "pearson", # correlation method
                 hist.col = "pink2",
                 density = TRUE,cor=TRUE,smooth=TRUE,  # show density plots
                 ellipses = FALSE,
                 cex.cor=0.6,cex=0.7,digits = 2, # show correlation ellipses
                 stars=TRUE,breaks=10,
                 main=paste0("Pearson corr ",namePlot," nb cells not 0 values : ",ncells4corrPlot, "(",cellsTotal, ")")); #show significance starts
     dev.off();

  };

  ################################### <-

  resDfcluster <- list();
  for (j in 1:length(inputList)){
    data <-  inputList[[j]];
    data<- as.data.frame(t(data));
    data$cellName=rownames(data);
    rownames(data) <- 1:nrow(data)
    dataName <-str_sub(names(inputList)[[j]] , start= -30);
    
   f_corr <- "correlation/thresholdDEG_1_minus1/";
    f_linearCorr <- "linearCorr/";
    dir.create(file.path(wd,f_results,f_corr), showWarnings = F);
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr), showWarnings = F);
    f_cluster <- paste0(gsub(".*_c","",gsub("__.*$","",dataName)),"/");
    f_data <- paste0(nameOutput,"/");
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr,f_data), showWarnings = F);
    dir.create(file.path(wd,f_results,f_corr,f_linearCorr,f_data,f_cluster), showWarnings = F);
    f_linearCorr <- "linearCorr/";
    f_data <- paste0(nameOutput,"/");



    for (i in 1:(length(colnames(data))-1)) {
      print(paste0("cluster nb: ",j," and genePair: ", i));
      obs_apoe.col <-  grep(gene2Corr,colnames(data));
      cellsTotal <- dim(data)[1];
      inputColumnsFC<- c(obs_apoe.col,i);
      df1 <- data[,inputColumnsFC];
      df1 <- df1[-which(df1[,1]==0),]
      df1 <- df1[-which(df1[,2]==0),]
      
      if (is.null(dim(df1))==FALSE && dim(df1)[1]>3  ){
        ncells4corrPlot <- dim(df1)[1]
        nameSavePlot <- paste0(colnames(data)[obs_apoe.col],"_",colnames(data)[i]) 
        namePlot <- paste0(colnames(data)[obs_apoe.col]," : ",colnames(data)[i]) 
        lineal_corrFunction.function(data, df1,nameOutput,inputColumnsFC,namePlot,nameSavePlot,f_cluster,ncells4corrPlot,cellsTotal,gene2Corr);

        if (i ==1){
          a.tmp <-cor.ci(round(df1[,c(1:2)],2),n.iter = 100,  p = 0.05,overlap = FALSE,poly = FALSE, method = "spearman",cex=0.7, plot=TRUE,main="Spearman corr");
          resDf <- as.data.frame(a.tmp$ci)[,c(2,4)];
          rownames(resDf) <- paste0(colnames(df1[,c(1:2)])[1],"-",colnames(df1[,c(1:2)])[2]);
          resDf$corr <- a.tmp$rho[2];
          resDf$p_value <- cor_to_p(a.tmp$rho[2], ncells4corrPlot, method = "spearman")$p;
          resDf$FDR <- p.adjust(cor_to_p(a.tmp$rho[2], ncells4corrPlot, method = "spearman")$p, method = "fdr", n = ncells4corrPlot);
          resDf$GenePairs <- rownames(resDf);
          resDf$Total_nb_cells <-cellsTotal;
          resDf$ncells_non0_used_corrPlot <-ncells4corrPlot;
          resDf$cluster <-gsub(".*_c","",gsub("__.*$","",dataName));
          resDf$cellType_res05 <-gsub("__.*$","",dataName);
          resDf$condition <-gsub(".*__","",gsub("__.$","",dataName));
          resDf$sex <-gsub(".*__","",dataName);
          resDfF <- resDf;
        } else{
          a.tmp2 <-cor.ci(round(df1[,c(1:2)],2),n.iter = 100,  p = 0.05,overlap = FALSE,poly = FALSE, method = "spearman",cex=0.7, plot=TRUE,main="Spearman corr");
          resDf2 <- as.data.frame(a.tmp2$ci)[,c(2,4)];
          rownames(resDf2) <- paste0(colnames(df1[,c(1:2)])[1],"-",colnames(df1[,c(1:2)])[2]);
          resDf2$corr <- a.tmp2$rho[2];
          resDf2$p_value <- cor_to_p(a.tmp2$rho[2], ncells4corrPlot, method = "spearman")$p;
          resDf2$FDR <- p.adjust(cor_to_p(a.tmp2$rho[2], ncells4corrPlot, method = "spearman")$p, method = "fdr", n = ncells4corrPlot);
          resDf2$GenePairs <- rownames(resDf2);
          resDf2$Total_nb_cells <-cellsTotal;
          resDf2$ncells_non0_used_corrPlot <-ncells4corrPlot;
          resDf2$cluster <-gsub(".*_c","",gsub("__.*$","",dataName));
          resDf2$cellType_res05 <-gsub("__.*$","",dataName);
          resDf2$condition <-gsub(".*__","",gsub("__.$","",dataName));
          resDf2$sex <-gsub(".*__","",dataName);
          resDfF <- bind_rows(resDfF,resDf2);

        };
      } else{
        if(i==1){
          resDfF <- data.frame(low.e=0, up.e=0, corr=0,  GenePairs="test",Total_nb_cells=0, ncells_non0_used_corrPlot=0,cluster="test", cellType_res05="test", condition="test", sex="test")
        } else{
          resDfF<- resDfF; 
        };
      };
      i <- i+1;
    };
  resDfF <-resDfF[!(resDfF$sex=="test"),];
  openxlsx::write.xlsx(resDfF, file = paste0(wd,f_results,f_corr, f_linearCorr,f_data,f_cluster,nameOutput,"_",dataName,"_Stats_perC_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
  resDfcluster[[j]]<- resDfF;
  names(resDfcluster)[j] <- dataName;

  j<-j+1;

  };
  assign(paste0(nameOutput,"_CorrDfcluster"),resDfcluster,.GlobalEnv);

for (ii in 1:(length(resDfcluster))){
    if(ii==1){
      df <- resDfcluster[[ii]];
      dfF=df;
    } else{
      df.tmp <- resDfcluster[[ii]];
      dfF<- bind_rows(dfF,df.tmp);

    };
    ii <- ii+1;
  };
 
  openxlsx::write.xlsx(dfF, file = paste0(wd,f_results,f_corr, f_linearCorr,f_data,nameOutput,"_corSttsAllC_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);

};




#############################################################################################################

#############################################################################################################
#############################################################################################################
#  4 Run the code to analyse the MTD list of genes per condition and find the corr of
#  each of those genes with another gene of interes  code
#############################################################################################################
#############################################################################################################

gene2Corr <- "APOE" # We will extract the correlation fo all the genes in genesInterestcond versus APOE.

MTDcorrelation.function(input_mitoXplorer_hsa,res_logCounts_percond_1,genesInterestcond, date,"cond_DEGs_1_all",gene2Corr)
MTDcorrelation.function(input_mitoXplorer_hsa,res_logCounts_percond_2,genesInterestcond, date,"cond_DEGs_2_all",gene2Corr)
MTDcorrelation.function(input_mitoXplorer_hsa,res_logCounts_percond_2,genesInterestcond_FC12, date,"cond_FC12DEGs_2_all",gene2Corr)
MTDcorrelation.function(input_mitoXplorer_hsa,res_logCounts_percond_1,genesInterestcond_FC12, date,"cond_FC12DEGs_1_all",gene2Corr)

save(cond_DEGs_1_all_resList_MTDcorr, , cond_FC12DEGs_2_all_resList_MTDcorr,
 cond_FC12DEGs_1_all_resList_MTDcorr,cond_DEGs_1_all_resList_MTDnot0_corr,cond_DEGs_2_all_resList_MTDnot0_corr,
 cond_FC12DEGs_2_all_resList_MTDnot0_corr, cond_FC12DEGs_1_all_resList_MTDnot0_corr,file=paste0(wd,f_results, f_corr,"cond_DEGs_all.RData"))

save(cond_DEGs_1_all_APOE_MTDcorr_df, cond_DEGs_2_all_APOE_MTDcorr_df, cond_FC12DEGs_2_all_APOE_MTDcorr_df,
 cond_FC12DEGs_1_all_APOE_MTDcorr_df,cond_DEGs_1_all_APOE_MTDcorr_not0Cells_df, cond_DEGs_2_all_APOE_MTDcorr_not0Cells_df,
  cond_FC12DEGs_2_all_APOE_MTDcorr_not0Cells_df, cond_FC12DEGs_1_all_APOE_MTDcorr_not0Cells_df,
file=paste0(wd,f_results, f_corr,"cond_DEGs_all_df.RData"))


save(cond_DEGs_1_all_List_inputCorr_all, cond_DEGs_2_all_List_inputCorr_all, cond_FC12DEGs_2_all_List_inputCorr_all,
 cond_FC12DEGs_1_all_List_inputCorr_all,cond_DEGs_1_all_List_inputCorr_not0, cond_DEGs_2_all_List_inputCorr_not0,
  cond_FC12DEGs_2_all_List_inputCorr_not0, cond_FC12DEGs_1_all_List_inputCorr_not0,
file=paste0(wd,f_results, f_corr,"cond_DEGs_all_inputCorr.RData"))
load(file=paste0(wd,f_results, f_corr,"cond_DEGs_all_inputCorr.RData")):

#http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
MTDcorrelation_not0_value.function("cond_2_FC15",cond_DEGs_2_all_List_inputCorr_not0,"070122");
MTDcorrelation_not0_value.function("cond_2_FC12",cond_FC12DEGs_2_all_List_inputCorr_not0,"070122");
MTDcorrelation_not0_value.function("cond_1_FC15",cond_DEGs_1_all_List_inputCorr_not0,"070122");
MTDcorrelation_not0_value.function("cond_1_FC12",cond_FC12DEGs_1_all_List_inputCorr_not0,"070122");

save(cond_1_FC15_CorrDfcluster, cond_1_FC15_CorrDfcluster, 
  cond_2_FC15_CorrDfcluster,cond_2_FC15_CorrDfcluster,
file=paste0(wd,f_results, f_corr,"correlationList_not0measurements_percondalone.RData"))




#corr plot with each of the other genes
####################################



###

#####################################################################
# 5 correlations with whatever dataframe containing a columns of genes and Processes
#####################################################################

 input.df <- as.data.frame(read_excel(paste0("C:/Users/Documents/toTest_Corr_otherGenes.xlsx"),sheet =1, col_names =T))[,c("Gene stable ID", "Gene" , "gene_function", "Chromosome" ,"Process")]  
 load(file=paste0(wd,f_results, f_corr,"cond_data_4corrAnalysis.RData"))

correlation_anyDF_genes.function <- function(genesDfInput,countsList, date,nameOutput,gene2Corr ){
  f_corr <-"correlation/thresholdDEG_1_minus1/";
  dir.create(file.path(wd,f_results,f_corr), showWarnings = F);
  dir.create(file.path(wd,f_results,f_corr,nameOutput,"/"), showWarnings = F);

  MTDcategories_perGene <- unique(genesDfInput[,c("Gene", "Process","gene_function")]);
  resList_MTDcorr <- list();
  resList_MTDnot0_corr <- list();
  resList_inputCorr_all <- list();
  resList_inputCorr_not0 <- list();


  for ( i in 1:length(names(countsList))){
  nameData <- gsub("__","_",gsub(".*_c","c",names(countsList)[i]));    
    a <-countsList[[i]]
    b <- a[which(rownames(a) %in% genesDfInput$Gene ),]; #row:genes, cols:cells
    b.not0 <- b[,-which(b==0)];  #analysing just cells where none of the interesting genes have a value =0.
    not0 <- apply(b, 1, function(x)length(which(x > 0)));
        
    b.expressed <- rownames(b)[which(not0/ncol(b) >= 0.25)]; #genes expressed in cells totally more thna a 25 %
    b.expressed <- b[which(rownames(b) %in% b.expressed ),]; #row:genes, cols:cells


    resList_inputCorr_all[[i]] <- b;
    resList_inputCorr_not0[[i]] <- b.not0;

    names(resList_inputCorr_all)[i] <- names(countsList)[i];
    names(resList_inputCorr_not0)[i] <- names(countsList)[i];


    if (is.null(dim(b.not0))==FALSE && dim(b.not0)[2]>4  ){
        c.not0 <-  rcorr(t(b.not0), type="spearman")$r;
        c.not0_pval <- rcorr(t(b.not0), type="spearman")$P;
        c.not0_ncells <- rcorr(t(b.not0), type="spearman")$n;
        d.not0 <- c.not0[which(rownames(c.not0)==gene2Corr),];
        d.not0_pval <- c.not0_pval[which(rownames(c.not0_pval)==gene2Corr),];
        d.not0_ncells <- c.not0_ncells[which(rownames(c.not0_ncells)==gene2Corr),];
       d.not0_fdr <- p.adjust(d.not0_pval, method = "fdr", n = length(d.not0_pval));
        e.not0 <- data.frame(Gene=names(d.not0), corr=d.not0,pval=d.not0_pval,
         FDR=d.not0_fdr ,ncellsTested=d.not0_ncells,name=names(countsList)[i],
          Cluster=gsub("__.*","",names(countsList)[i]), 
          LA=gsub(".*__","",gsub("__M","",gsub("__F","",names(countsList)[i]))) ,
          Sex=gsub(".*__","",names(countsList)[i]));
         e.not0 <- left_join(e.not0,MTDcategories_perGene,by="Gene");
        resList_MTDnot0_corr[[i]] <-e.not0;
        names(resList_MTDnot0_corr)[i] <- names(countsList)[i];

        openxlsx::write.xlsx(e.not0, file = paste0(wd,f_results,f_corr,nameOutput ,"/",f_per_group,nameData,"_MTDnot0Cells_corr_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
        colnames(e.not0)[c(2:5)] <- paste0(colnames(e.not0)[c(2:5)], "_",names(countsList)[i] )

    };

    if (is.null(dim(b))==FALSE && dim(b)[2]>4  ){
      c <-  rcorr(t(b), type="spearman")$r;
      c_pval <- rcorr(t(b), type="spearman")$P;
      c_ncells <- rcorr(t(b), type="spearman")$n;
      d <- c[which(rownames(c)==gene2Corr),];
      d_pval <- c_pval[which(rownames(c_pval)==gene2Corr),];
      d_ncells <- c_ncells[which(rownames(c_ncells)==gene2Corr),];
      d_fdr <- p.adjust(d_pval, method = "fdr", n = length(d_pval));
      e <- data.frame(Gene=names(d), corr=d,pval=d_pval,
        FDR=d_fdr ,ncellsTested=d_ncells,name=names(countsList)[i],
        Cluster=gsub("__.*","",names(countsList)[i]), 
        LA=gsub(".*__","",gsub("__M","",gsub("__F","",names(countsList)[i]))) ,
        Sex=gsub(".*__","",names(countsList)[i]));
       e <- left_join(e,MTDcategories_perGene,by="Gene");

    resList_MTDcorr[[i]] <-e;
    names(resList_MTDcorr)[i] <- names(countsList)[i];

};


  i <-i+1;
  print(i)
  };
 
  resList_inputCorr_not0 <- resList_inputCorr_not0[lapply(resList_inputCorr_not0,length)>0]

  assign(paste0(nameOutput,"_List_inputCorr_not0"),resList_inputCorr_not0,.GlobalEnv);

  for (ii in 1:(length(resList_MTDcorr))){
    if(ii==1){
      df <- resList_MTDcorr[[ii]];
      dfF=df;
    } else{
      df.tmp <- resList_MTDcorr[[ii]];
      dfF<- bind_rows(dfF,df.tmp);

    };
    ii <- ii+1;
  };

 
  for (iii in 1:(length(resList_MTDnot0_corr))){
    if(iii==1){
      df2 <- resList_MTDnot0_corr[[iii]];
      dfF2=df2;
    } else{
      df2.tmp <- resList_MTDnot0_corr[[iii]];
      dfF2<- bind_rows(dfF2,df2.tmp);
    };
    iii <- iii+1;
  };
  assign(paste0(nameOutput,"_",gene2Corr,"_MTDcorr_not0Cells_df"),dfF2,.GlobalEnv);

  resListF <- list(corr_allCells=dfF,corr_cellsNot0=dfF2);
  openxlsx::write.xlsx(resListF, file = paste0(wd,f_results,f_corr,nameOutput ,"/",nameData,"_MTDcorr_",date, ".xlsx"),rowNames =FALSE, colNames =TRUE);
 
};


correlation_anyDF_genes.function(Test_corr_dev_others,res_logCounts_percond_2, date,"cond_Test_2_w0",gene2Corr);
correlation_anyDF_genes.function(Test_corr_dev_others,res_logCounts_percond_2, date,"cond_Test_2_w0_subset",gene2Corr);
correlation_anyDF_genes.function(Test_corr_dev_others,res_logCounts_percond_1, date,"cond_Test_1_w0",gene2Corr);

save(  cond_Test_2_w0_APOE_MTDcorr_not0Cells_df, cond_Test_1_w0_APOE_MTDcorr_not0Cells_df,
file=paste0(wd,f_results, f_corr,"Test_df_percond.RData"))

save(cond_Test_2_w0_List_inputCorr_not0, cond_Test_1_w0_List_inputCorr_not0,
file=paste0(wd,f_results, f_corr,"Test_inputCorr_percond.RData"))

load(file=paste0(wd,f_results, f_corr,"Test_inputCorr_percond.RData"))
MTDcorrelation_not0_value.function("Test_2_wo_0",cond_Test_2_w0_List_inputCorr_not0,date);
MTDcorrelation_not0_value.function("Test_1_wo_0",cond_Test_1_w0_List_inputCorr_not0,date);

save(Test_2_wo_0_CorrDfcluster, Test_1_wo_0_CorrDfcluster, 
file=paste0(wd,f_results, f_corr,"correlationList_not0measurements_Test_percond.RData"))

