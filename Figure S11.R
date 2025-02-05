library(monocle3)
library(VGAM) 
library(viridis)
library(stringr)
library(tibble) 
library(reticulate)
library(pheatmap)
library(dplyr)
library(magrittr)
library(ggplot2)
library(garnett)
library(ggsci)
library(ggpubr)
library(viridis)
library(leidenbase)
library(sf)
library(Rcpp)
library(magrittr)
library(scales)
library(ggsci)
library(reshape)
library(m3addon)


simple_theme <-  theme(text = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), legend.position = "none")  ### theme to remove legends and axis/font

simple_themeL <-  theme(text = element_text(size=25),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank())   #### theme show the legend with large font

simple_themeT <-  theme(text = element_text(size=25),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank(),legend.position="none")

#session info
writeLines(capture.output(sessionInfo()), "sessionInfo040724HashHSC.txt")

#location
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))

#Load cds
folders<-'/Users'
cds<-load_cellranger_data_h5(folders)

table(fData(cds)$feature_type)
prot<-cds[fData(cds)$feature_type=="Antibody Capture",]
gex<-cds[fData(cds)$feature_type=="Gene Expression",]
protm<-log10(exprs(prot)+0.001)
rownames(protm)<-paste0("Log_", rownames(protm))
colData(gex)<-cbind(colData(gex), t(protm))
cds<-gex

cds <- detect_genes(cds, min_expr=0.1)     
cds <- cds[,colData(cds)$num_genes_expressed > 1000]

# Remove low UMIs 
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
qplot(colData(cds)$n.umis, geom="density")
qplot(log10(colData(cds)$n.umis), geom="density")
summary(colData(cds)$n.umis)

cds <- cds[,Matrix::colSums(counts(cds)) > 3000]
cds <- cds[,Matrix::colSums(counts(cds)) < 50000] 

#m3addon
cds<-preprocess_cds(cds)
cds<-m3addon::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", verbose=T, cores=detectCores()/2)

#Remove cells with high content of mitochondrial genes
#List of all mitochondrial genes by Ensembl ID
mito_genes_updated <- c("ENSMUSG00000064351", "ENSMUSG00000064354", "ENSMUSG00000064370", "ENSMUSG00000064357")

num_mito <- Matrix::colSums(counts(cds[mito_genes_updated]))
cds$n.mito <- num_mito

colData(cds)$n.umis <- Matrix::colSums(counts(cds))
perc_mito <- 100*(cds$n.mito / cds$n.umis)
cds$perc_mito_umis <- perc_mito

cds_mito_filter <- cds[,pData(cds)$perc_mito_umis < 10]
cds<- cds_mito_filter

#cluster cells
set.seed(1)
cds<-cluster_cells(cds, resolution = 0.00001, random.seed=1)
plot_cells(cds, color_cells_by = "cluster",
           label_cell_groups=FALSE,cell_size = 1) 

#Annotation by Garnette
library(org.Mm.eg.db)

#Location 
RES_DIR <- file.path("/Users")

### check marker genes
marker_check <- check_markers(cds, marker_file_path,
                              db=org.Mm.eg.db,
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)

hspc_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker_file_path,
                                         db=org.Mm.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 500,
                                         marker_file_gene_id_type = "SYMBOL")


colData(cds)$garnett_cluster = monocle3::clusters(cds)

cds <- classify_cells(cds, hspc_classifier,
                      db = org.Mm.eg.db,
                      cluster_extend = TRUE,
                      cluster_extend_max_frac_unknown = 0.95,
                      cluster_extend_max_frac_incorrect = 0.25,
                      cds_gene_id_type = "ENSEMBL")


feature_genes <- get_feature_genes(hspc_classifier,
                                   node = "root",
                                   db = org.Mm.eg.db,convert_ids = TRUE)
head(feature_genes)

get_classifier_references(hspc_classifier)

### re-assign colors and order to cell_types
cell_type_color <- c( "Ptprc-Cdh5+"="darkgreen", "Ptprc+Cdh5-"="orange", "Unknown"="white")
colData(cds)$cell_type <- factor(colData(cds)$cell_type, levels = c("Ptprc-Cdh5+","Ptprct+Cdh5-","Unknown"))
colData(cds)$cluster_ext_type <- factor(colData(cds)$cluster_ext_type, levels = c("Ptprc-Cdh5+","Ptprc+Cdh5-","Unknown"))

### plot by cluster-extended cell typing                                                                                                                                                                 "AEC2-PreHE","EHT"))
plot_cells(cds,
           color_cells_by = "cluster_ext_type",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1, cell_size = 0.6) +
  scale_color_manual(values = cell_type_color) + simple_themeL


table(pData(cds)$cluster_ext_type)

cds <- cds[,(colData(cds)$cluster_ext_type == "Ptprc+Cdh5-"),]

#Demultiplex the sample
df<-melt(as.data.frame(colData(cds)[,colnames(colData(cds)) %in% c("Log_B0301"
                                                                   ,"Log_B0302"
                                                                   ,"Log_B0303"
                                                                   ,"Log_B0304"
                                                                   ,"Log_B0305"
                                                                   ,"Log_B0306"
                                                                   ,"Log_B0307")]))


thresh<-1.8 #higher
ggplot(df, aes(x=variable, y=value, fill=variable))+geom_violin()+geom_hline(yintercept = thresh)

cds$HTO1<-cds$Log_B0301>thresh 
cds$HTO2<-cds$Log_B0302>thresh 
cds$HTO3<-cds$Log_B0303>thresh 
cds$HTO4<-cds$Log_B0304>thresh 
cds$HTO5<-cds$Log_B0305>thresh 
cds$HTO6<-cds$Log_B0306>thresh 
cds$HTO7<-cds$Log_B0307>thresh 

#Colony1 : HTO1=TRUE, HTO2,3,4,5,6,7=FALSE
Colony1Preliminary <- cds[,(colData(cds)$HTO1 == "TRUE"),] 
Colony1Preliminary
Colony1a <- Colony1Preliminary[,(colData(Colony1Preliminary)$HTO2 == "FALSE"),] 
Colony1b <- Colony1a[,(colData(Colony1a)$HTO3 == "FALSE"),] 
Colony1c <- Colony1b[,(colData(Colony1b)$HTO4 == "FALSE"),] 
Colony1d <- Colony1c[,(colData(Colony1c)$HTO5 == "FALSE"),] 
Colony1e <- Colony1d[,(colData(Colony1d)$HTO6 == "FALSE"),]
Colony1 <- Colony1e[,(colData(Colony1e)$HTO7 == "FALSE"),]
Colony1
colData(Colony1)
table(Colony1$HTO1)

#Colony2 : HTO2=TRUE, HTO1,3,4,5,6,7=FALSE
Colony2Preliminary <- cds[,(colData(cds)$HTO2 == "TRUE"),] 
Colony2Preliminary
Colony2a <- Colony2Preliminary[,(colData(Colony2Preliminary)$HTO1 == "FALSE"),] 
Colony2b <- Colony2a[,(colData(Colony2a)$HTO3 == "FALSE"),] 
Colony2c <- Colony2b[,(colData(Colony2b)$HTO4 == "FALSE"),] 
Colony2d <- Colony2c[,(colData(Colony2c)$HTO5 == "FALSE"),] 
Colony2e <- Colony2d[,(colData(Colony2d)$HTO6 == "FALSE"),]
Colony2 <- Colony2e[,(colData(Colony2e)$HTO7 == "FALSE"),]
Colony2
colData(Colony2)
table(Colony2$HTO2)

#Colony3 : HTO3=TRUE, HTO1,2,4,5,6,7=FALSE
Colony3Preliminary <- cds[,(colData(cds)$HTO3 == "TRUE"),] 
Colony3Preliminary
Colony3a <- Colony3Preliminary[,(colData(Colony3Preliminary)$HTO1 == "FALSE"),] 
Colony3b <- Colony3a[,(colData(Colony3a)$HTO2 == "FALSE"),] 
Colony3c <- Colony3b[,(colData(Colony3b)$HTO4 == "FALSE"),] 
Colony3d <- Colony3c[,(colData(Colony3c)$HTO5 == "FALSE"),] 
Colony3e <- Colony3d[,(colData(Colony3d)$HTO6 == "FALSE"),]
Colony3 <- Colony3e[,(colData(Colony3e)$HTO7 == "FALSE"),]
Colony3
colData(Colony3)
table(Colony3$HTO3)

#Colony4 : HTO4=TRUE, HTO1,2,3,5,6,7=FALSE
Colony4Preliminary <- cds[,(colData(cds)$HTO4 == "TRUE"),] 
Colony4Preliminary
Colony4a <- Colony4Preliminary[,(colData(Colony4Preliminary)$HTO1 == "FALSE"),] 
Colony4b <- Colony4a[,(colData(Colony4a)$HTO2 == "FALSE"),] 
Colony4c <- Colony4b[,(colData(Colony4b)$HTO3 == "FALSE"),] 
Colony4d <- Colony4c[,(colData(Colony4c)$HTO5 == "FALSE"),] 
Colony4e <- Colony4d[,(colData(Colony4d)$HTO6 == "FALSE"),]
Colony4 <- Colony4e[,(colData(Colony4e)$HTO7 == "FALSE"),]
Colony4
colData(Colony4)
table(Colony4$HTO4)

#Colony5 : HTO5=TRUE, HTO1,2,3,4,6,7=FALSE
Colony5Preliminary <- cds[,(colData(cds)$HTO5 == "TRUE"),] 
Colony5Preliminary
Colony5a <- Colony5Preliminary[,(colData(Colony5Preliminary)$HTO1 == "FALSE"),] 
Colony5b <- Colony5a[,(colData(Colony5a)$HTO2 == "FALSE"),] 
Colony5c <- Colony5b[,(colData(Colony5b)$HTO3 == "FALSE"),] 
Colony5d <- Colony5c[,(colData(Colony5c)$HTO4 == "FALSE"),] 
Colony5e <- Colony5d[,(colData(Colony5d)$HTO6 == "FALSE"),]
Colony5 <- Colony5e[,(colData(Colony5e)$HTO7 == "FALSE"),]
Colony5
colData(Colony5)
table(Colony5$HTO5)

#Colony6 : HTO6=TRUE, HTO1,2,3,4,5,7=FALSE
Colony6Preliminary <- cds[,(colData(cds)$HTO6 == "TRUE"),] 
Colony6Preliminary
Colony6a <- Colony6Preliminary[,(colData(Colony6Preliminary)$HTO1 == "FALSE"),] 
Colony6b <- Colony6a[,(colData(Colony6a)$HTO2 == "FALSE"),] 
Colony6c <- Colony6b[,(colData(Colony6b)$HTO3 == "FALSE"),] 
Colony6d <- Colony6c[,(colData(Colony6c)$HTO4 == "FALSE"),] 
Colony6e <- Colony6d[,(colData(Colony6d)$HTO5 == "FALSE"),]
Colony6 <- Colony6e[,(colData(Colony6e)$HTO7 == "FALSE"),]
Colony6
colData(Colony6)
table(Colony6$HTO6)

#Colony7 : HTO7=TRUE, HTO1,2,3,4,5,6=FALSE
Colony7Preliminary <- cds[,(colData(cds)$HTO7 == "TRUE"),] 
Colony7Preliminary
Colony7a <- Colony7Preliminary[,(colData(Colony7Preliminary)$HTO1 == "FALSE"),] 
Colony7b <- Colony7a[,(colData(Colony7a)$HTO2 == "FALSE"),] 
Colony7c <- Colony7b[,(colData(Colony7b)$HTO3 == "FALSE"),] 
Colony7d <- Colony7c[,(colData(Colony7c)$HTO4 == "FALSE"),] 
Colony7e <- Colony7d[,(colData(Colony7d)$HTO5 == "FALSE"),]
Colony7 <- Colony7e[,(colData(Colony7e)$HTO6 == "FALSE"),]
Colony7
colData(Colony7)
table(Colony7$HTO7)

Colony1
Colony2
Colony3
Colony4
Colony5
Colony6
Colony7

#combine
Colony1$identifier<-("ESAM+HSC")
Colony2$identifier<-("ESAM+HSC")
Colony3$identifier<-("ESAM+HSC")
Colony4$identifier<-("ESAM+HSC")
Colony5$identifier<-("ESAM+HSC")
Colony6$identifier<-("ESAM-HSC")
Colony7$identifier<-("ESAM-HSC")

Colony1$identity<-("B301")
Colony2$identity<-("B302")
Colony3$identity<-("B303")
Colony4$identity<-("B304")
Colony5$identity<-("B305")
Colony6$identity<-("B306")
Colony7$identity<-("B307")

#combine
big_cds <- combine_cds(list(Colony1,Colony2,Colony3,Colony4,Colony5,Colony6,Colony7), keep_all_genes=TRUE)

##pre-process  (Here, num_dim=0)
big_cds <- preprocess_cds(big_cds)

#Align : remove batch effects in combination cds
big_cds = align_cds(big_cds, alignment_group = "identifier")

##Reduce the dimensions using UMAP
big_cds <- reduce_dimension(big_cds)

plot_cells(big_cds, color_cells_by = "identity", label_cell_groups = F, cell_size = 1) #+ simple_theme

group_type_color2 <- c("ESAM+HSC"= "red","ESAM-HSC"= "gray")

big_cds$identifier<-factor(big_cds$identifier, levels =  c("ESAM+HSC","ESAM-HSC"))

plot_cells(big_cds,
           color_cells_by = "identifier",
           label_cell_groups=FALSE,cell_size = 1) + scale_color_manual(values = group_type_color2)

### plot UMI and genes per cell for each sample
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
                        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
                        panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_text(size=15))

df1 <- data.frame("identity"=colData(big_cds)$identity, "n.umis"=colData(big_cds)$n.umis)

df1$identity <- factor(df1$identity, levels = c("B301","B302","B303","B304","B305","B306","B307"))

ggplot(df1, aes(x=identity, y=n.umis, fill=identity)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(1000,100000), labels=number) + simple_theme2 

df2 <- data.frame("identity"=colData(big_cds)$identity, "num_genes_expressed"=colData(big_cds)$num_genes_expressed)

df2$identity <- factor(df2$identity, levels = c("B301","B302","B303","B304","B305","B306","B307"))

ggplot(df2, aes(x=identity, y=num_genes_expressed, fill=identity)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 

#change the name
cds<-big_cds

#location
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))
#Save
saveRDS(cds, file.path(RES_DIR, "Figure S11.RDS")) 

set.seed(1)
cds <- cluster_cells(cds, resolution=1e-2, random_seed = 1)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size = 2, group_label_size=20)+simple_theme

#Figure S11B-C
plot_cells(cds, color_cells_by = "identity", label_cell_groups = F, cell_size = 2) 
plot_cells(cds, color_cells_by = "identity", label_cell_groups = F, cell_size = 2) + simple_theme
table(cds$identity)

group_type_color2 <- c("ESAM+HSC"= "red","ESAM-HSC"= "gray")
cds$identifier<-factor(cds$identifier, levels =  c("ESAM+HSC","ESAM-HSC"))
plot_cells(cds,
           color_cells_by = "identifier",
           label_cell_groups=FALSE,cell_size = 1.5) + scale_color_manual(values = group_type_color2)+ simple_theme
table(cds$identifier)

#add cluster to colData
cluster <- clusters(cds)
pData(cds)$cluster <- cluster   
colData(cds)
head(pData(cds))
table(pData(cds)$cluster)

#Figure S11D-H
#Dormancy score (Cell . 2017 May 18;169(5):807-823.e19.)
RES_DIR <- file.path("/Users")

estimate_score <- function(cds,gene_markers){
  cds_gene = cds[fData(cds)$gene_short_name %in% gene_markers,] 
  aggregate_gene = exprs(cds_gene)
  aggregate_gene = Matrix::t(Matrix::t(aggregate_gene) / pData(cds_gene)$Size_Factor)
  aggregate_gene = Matrix::colSums(aggregate_gene)
  pData(cds)$gene_score = log(aggregate_gene +1) ##not sure if this is really needed?
  return(cds)
}

gene_markers <- c(read.csv(file.path(RES_DIR, "Dormant.csv"))) ##kegg database list for hallmark oxphos genes, change it to whatever gene list
cds <- estimate_score(cds, gene_markers = gene_markers$Dormant)
plot_cells(cds, color_cells_by = 'gene_score', cell_size = 2, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

#violin plots
df <- data.frame(colData(cds)$identifier, colData(cds)$gene_score)
names(df) <- c("identifier", "gene_score")
write.csv(df, "~/gene_cluster.csv")

df <- data.frame(read.csv("~/gene_cluster.csv"))
df$identifier<-factor(df$identifier, levels = c("ESAM+HSC","ESAM-HSC"))
p <- ggplot(df, aes(x= identifier, y=gene_score, fill = identifier)) + geom_violin(trim = FALSE)  + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + theme_bw()+ scale_fill_manual(values = cell_type_color1)
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0.5,0.5,-0.5,-1), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=18))
p

#Wilcoxon test
p <- ggplot(df, aes(x= identifier, y=gene_score, fill = identifier)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.color=NA)  
p <- p + stat_compare_means() + stat_mean()  
p


#serial transplantaable score: Nature volume 583, pages585–589 (2020)
#commonly upregulated genes: Nature Cell Biology volume 26, pages181–193 (2024)
#Mouse Gene Set: HALLMARK_MYC_TARGETS_V2, MM2888

### Estimate proliferation index (Cell Rep. 2021 Sep 14;36(11):109675. doi: 10.1016/j.celrep.2021.109675.)
estimate_cell_cycle <- function(cds, g1s_markers, g2m_markers){
  cds_g1s = cds[fData(cds)$gene_short_name %in% g1s_markers,]
  aggregate_g1s_expression = exprs(cds_g1s)
  aggregate_g1s_expression = Matrix::t(Matrix::t(aggregate_g1s_expression) / pData(cds_g1s)$Size_Factor)
  aggregate_g1s_expression = Matrix::colSums(aggregate_g1s_expression)
  
  cds_g2m = cds[fData(cds)$gene_short_name %in% g2m_markers,]
  aggregate_g2m_expression = exprs(cds_g2m)
  aggregate_g2m_expression = Matrix::t(Matrix::t(aggregate_g2m_expression) / pData(cds_g2m)$Size_Factor)
  aggregate_g2m_expression = Matrix::colSums(aggregate_g2m_expression)
  
  pData(cds)$g1s_score = log(aggregate_g1s_expression+1)
  pData(cds)$g2m_score = log(aggregate_g2m_expression+1)
  pData(cds)$proliferation_index = log(aggregate_g1s_expression + aggregate_g2m_expression + 1)
  return(cds)
}
s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1",
             "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", 
             "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", 
             "Chaf1b", "Brip1", "E2f8")
g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo",
               "Cenpf", "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1",
               "Kif20b", "Hjurp", "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8",
               "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")

cds <- estimate_cell_cycle(cds, g1s_markers = s.genes, g2m_markers = g2m.genes)
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cds, color_cells_by = 'proliferation_index', cell_size = 2, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

#violin plots
df <- data.frame(colData(cds)$identifier, colData(cds)$proliferation_index)
names(df) <- c("identifier", "proliferation_index")
write.csv(df, "~/gene_cluster.csv")

df <- data.frame(read.csv("~/gene_cluster.csv"))
df$identifier<-factor(df$identifier, levels = c("ESAM+HSC","ESAM-HSC"))
p <- ggplot(df, aes(x= identifier, y=proliferation_index, fill = identifier)) + geom_violin(trim = FALSE)  + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + theme_bw()+ scale_fill_manual(values = cell_type_color1)
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0.5,0.5,-0.5,-0.3), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=18))
p

#Wilcoxon test
p <- ggplot(df, aes(x= identifier, y=proliferation_index, fill = identifier)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.color=NA)  
p <- p + stat_compare_means() + stat_mean()  
p

#Figure S11I
#q<0.05
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

set.seed(1) #set.seed and random seed set here to avoid variable output due to random number generator in function
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.00005, random_seed = 1)

### plot gene modules by cluster in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("cluster", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12)

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE, legend=F ,show_rownames = F, show_colnames = F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12)

### plot gene modules by colony type in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$identifier)  
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
colnames(agg_mat) <- stringr::str_c("", colnames(agg_mat))

pheatmap::pheatmap(agg_mat[c("2","4","1","3"),c("ESAM+HSC", "ESAM-HSC")]
                   ,cluster_rows=FALSE, cluster_cols=FALSE
                   ,scale="column", clustering_method="ward.D2"
                   ,fontsize=12, show_rownames = T, show_colnames = T, legend=T, cellheight = 50, cellwidth = 20) 

#pulling out genes in each module
mod1_genes<-gene_module_df%>%dplyr::filter(module==1)
names(mod1_genes)[1]<-"gene_id"
mod1_genes_names<-merge(mod1_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod1_genes_names<-mod1_genes_names[,c(1,9,10,12,14)]
write.csv(mod1_genes_names, "/Users/mod1.csv")

#use the code for the rest of modules
