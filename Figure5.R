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

writeLines(capture.output(sessionInfo()), "sessionInfoHashHSC.txt")

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
cdsf<-cds
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
cds<-cluster_cells(cds, resolution = 0.0001, random.seed=1)
plot_cells(cds, color_cells_by = "cluster", cell_size = 0.5)

#Annotation by Garnette
library(org.Mm.eg.db)

#Location 
RES_DIR <- file.path("/Users")

marker_file_path <- "/Users/CD45_selection.txt"  

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
cell_type_color <- c("Ptprc+Hbb-bt-" = "purple3", "Ptprc-Cdh5+"="darkgreen", "Hbb-bt+Cdh5-"="orange", "Unknown"="white")
colData(cds)$cell_type <- factor(colData(cds)$cell_type, levels = c("Ptprc+Hbb-bt-","Ptprc-Cdh5+","Hbb-bt+Cdh5-","Unknown"))
colData(cds)$cluster_ext_type <- factor(colData(cds)$cluster_ext_type, levels = c("Ptprc+Hbb-bt-","Ptprc-Cdh5+","Hbb-bt+Cdh5-","Unknown"))

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

cds <- cds[,(colData(cds)$cluster_ext_type == "Ptprc+Hbb-bt-"),] 


#Demultiplex the sample
df<-melt(as.data.frame(colData(cds)[,colnames(colData(cds)) %in% c("Log_A11_well"
                                                                   ,"Log_D4_well"
                                                                   ,"Log_E1_well"
                                                                   ,"Log_F10_well"
                                                                   ,"Log_G6_well"
                                                                   ,"Log_G8_well")]))


thresh<-1.5 #higher
ggplot(df, aes(x=variable, y=value, fill=variable))+geom_violin()+geom_hline(yintercept = thresh)

cds$HTO1<-cds$Log_A11_well>thresh 
cds$HTO2<-cds$Log_D4_well>thresh 
cds$HTO3<-cds$Log_E1_well>thresh 
cds$HTO4<-cds$Log_F10_well>thresh 
cds$HTO5<-cds$Log_G6_well>thresh 
cds$HTO6<-cds$Log_G8_well>thresh 

#A11 : HTO1=TRUE, HTO2,3,4,5,6=FALSE
A11Preliminary <- cds[,(colData(cds)$HTO1 == "TRUE"),] 
A11Preliminary
A11a <- A11Preliminary[,(colData(A11Preliminary)$HTO2 == "FALSE"),] 
A11b <- A11a[,(colData(A11a)$HTO3 == "FALSE"),] 
A11c <- A11b[,(colData(A11b)$HTO4 == "FALSE"),] 
A11d <- A11c[,(colData(A11c)$HTO5 == "FALSE"),] 
A11 <- A11d[,(colData(A11d)$HTO6 == "FALSE"),]
A11
colData(A11)
table(A11$HTO1)

#D4 : HTO2=TRUE, HTO1,3,4,5,6=FALSE
D4Preliminary <- cds[,(colData(cds)$HTO2 == "TRUE"),] 
D4Preliminary
D4a <- D4Preliminary[,(colData(D4Preliminary)$HTO1 == "FALSE"),] 
D4b <- D4a[,(colData(D4a)$HTO3 == "FALSE"),] 
D4c <- D4b[,(colData(D4b)$HTO4 == "FALSE"),] 
D4d <- D4c[,(colData(D4c)$HTO5 == "FALSE"),] 
D4 <- D4d[,(colData(D4d)$HTO6 == "FALSE"),] 
D4
colData(D4)
table(D4$HTO2)

#E1 : HTO3=TRUE, HTO1,2,4,5,6=FALSE
E1Preliminary <- cds[,(colData(cds)$HTO3 == "TRUE"),] 
E1Preliminary
E1a <- E1Preliminary[,(colData(E1Preliminary)$HTO1 == "FALSE"),] 
E1b <- E1a[,(colData(E1a)$HTO2 == "FALSE"),] 
E1c <- E1b[,(colData(E1b)$HTO4 == "FALSE"),] 
E1d <- E1c[,(colData(E1c)$HTO5 == "FALSE"),] 
E1 <- E1d[,(colData(E1d)$HTO6 == "FALSE"),] 
E1
colData(E1)
table(E1$HTO3)

#F10 : HTO4=TRUE, HTO1,2,3,5,6=FALSE
F10Preliminary <- cds[,(colData(cds)$HTO4 == "TRUE"),] 
F10Preliminary
F10a <- F10Preliminary[,(colData(F10Preliminary)$HTO1 == "FALSE"),] 
F10b <- F10a[,(colData(F10a)$HTO2 == "FALSE"),] 
F10c <- F10b[,(colData(F10b)$HTO3 == "FALSE"),] 
F10d <- F10c[,(colData(F10c)$HTO5 == "FALSE"),] 
F10 <- F10d[,(colData(F10d)$HTO6 == "FALSE"),] 
F10
colData(F10)
table(F10$HTO4)

#G6 : HTO5=TRUE, HTO1,2,3,4,6=FALSE
G6Preliminary <- cds[,(colData(cds)$HTO5 == "TRUE"),] 
G6Preliminary
G6a <- G6Preliminary[,(colData(G6Preliminary)$HTO1 == "FALSE"),] 
G6b <- G6a[,(colData(G6a)$HTO2 == "FALSE"),] 
G6c <- G6b[,(colData(G6b)$HTO3 == "FALSE"),] 
G6d <- G6c[,(colData(G6c)$HTO4 == "FALSE"),] 
G6 <- G6d[,(colData(G6d)$HTO6 == "FALSE"),] 
G6
colData(G6)
table(G6$HTO5)

#G8 : HTO6=TRUE, HTO1,2,3,4,5=FALSE
G8Preliminary <- cds[,(colData(cds)$HTO6 == "TRUE"),] 
G8Preliminary
G8a <- G8Preliminary[,(colData(G8Preliminary)$HTO1 == "FALSE"),] 
G8b <- G8a[,(colData(G8a)$HTO2 == "FALSE"),] 
G8c <- G8b[,(colData(G8b)$HTO3 == "FALSE"),] 
G8d <- G8c[,(colData(G8c)$HTO4 == "FALSE"),] 
G8 <- G8d[,(colData(G8d)$HTO5 == "FALSE"),] 
G8
colData(G8)
table(G8$HTO6)

#combine
A11$identifier<-("No engraft")
D4$identifier<-("Engraft")
E1$identifier<-("No engraft")
F10$identifier<-("No engraft")
G6$identifier<-("No engraft")
G8$identifier<-("No engraft")

A11$identity<-("A11")
D4$identity<-("D4")
E1$identity<-("E1")
F10$identity<-("F10")
G6$identity<-("G6")
G8$identity<-("G8")

#combine
big_cds <- combine_cds(list(A11, D4, E1, F10, G6, G8), keep_all_genes=TRUE)

##pre-process  
big_cds <- preprocess_cds(big_cds)

#Align : remove batch effects in combination cds
big_cds = align_cds(big_cds, alignment_group = "identifier")

##Reduce the dimensions using UMAP
big_cds <- reduce_dimension(big_cds)

plot_cells(big_cds, color_cells_by = "identity", label_cell_groups = F, cell_size = 1) #+ simple_theme

group_type_color2 <- c("Engraft"= "red","No engraft"= "gray")

big_cds$identifier<-factor(big_cds$identifier, levels =  c("Engraft","No engraft"))

plot_cells(big_cds,
           color_cells_by = "identifier",
           label_cell_groups=FALSE,cell_size = 1) + scale_color_manual(values = group_type_color2)

### plot UMI and genes per cell for each sample
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
                        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
                        panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_text(size=15))

df1 <- data.frame("identity"=colData(big_cds)$identity, "n.umis"=colData(big_cds)$n.umis)

df1$identity <- factor(df1$identity, levels = c("A11","D4","E1","F10","G6","G8"))

ggplot(df1, aes(x=identity, y=n.umis, fill=identity)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(1000,100000), labels=number) + simple_theme2 

df2 <- data.frame("identity"=colData(big_cds)$identity, "num_genes_expressed"=colData(big_cds)$num_genes_expressed)

df2$identity <- factor(df2$identity, levels = c("A11","D4","E1","F10","G6","G8"))

ggplot(df2, aes(x=identity, y=num_genes_expressed, fill=identity)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 


cds<-big_cds

#location
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))

#Save
#whole
saveRDS(cds, file.path(RES_DIR, "Figure5.RDS")) 
#Read
cds <- readRDS(file.path(RES_DIR, "Figure5.RDS"))


#Figure S8C
plot_cells(cds, genes=c("Esam"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Fgd5"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Hlf"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("H19"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Procr"), cell_size = 4, label_cell_groups = F) + simple_themeT
plot_cells(cds, genes=c("Mecom"), cell_size = 4, label_cell_groups = F) + simple_themeT
plot_cells(cds, genes=c("Cdkn1c"), cell_size = 4, label_cell_groups = F) + simple_themeT
plot_cells(cds, genes=c("Mllt3"), cell_size = 4, label_cell_groups = F) + simple_themeT

plot_cells(cds, genes=c("Gata1"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Klf1"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Elane"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Cd48"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Flt3"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Il7r"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Ccl5"),cell_size = 4, label_cell_groups = F)+ simple_themeT

#Figure 5C and D, Figure S8B
set.seed(1)
cds <- cluster_cells(cds, resolution=1.5e-2, random_seed = 1)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size = 2, group_label_size=20)+simple_theme

plot_cells(cds, color_cells_by = "identity", label_cell_groups = F, cell_size = 2) + simple_theme
table(cds$identity)

group_type_color2 <- c("Engraft"= "red","No engraft"= "gray")
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


#Figure5H, Figure S8C
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

set.seed(1) #set.seed and random seed set here to avoid variable output due to random number generator in function
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.00025, random_seed = 1)

### plot gene modules by cluster in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("cluster", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE, legend=F ,show_rownames = F, show_colnames = F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12)

#show annotation
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12)

### plot gene modules by colony type in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$identifier)  
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
colnames(agg_mat) <- stringr::str_c("", colnames(agg_mat))

pheatmap::pheatmap(agg_mat[c("1","4","3","2"),c("Engraft", "No engraft")]
                   ,cluster_rows=FALSE, cluster_cols=FALSE
                   ,scale="column", clustering_method="ward.D2"
                   ,fontsize=12, show_rownames = F, show_colnames = F, legend=F, cellheight = 100, cellwidth = 50) 
#show annotation
pheatmap::pheatmap(agg_mat[c("1","4","3","2"),c("Engraft", "No engraft")], cluster_rows=FALSE, cluster_cols=FALSE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12, legend=T, show_rownames = T, show_colnames = T, cellheight = 70, cellwidth = 50) 


#### plot gene modules in UMAP 
plot_cells(cds, cell_size = 2,
           genes=gene_module_df %>% filter(module %in% c(1,2,3,4)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           scale_to_range = TRUE,
           min_expr = 60) +scale_color_gradientn(colors=c("gray80", "darkgreen", "darkgreen"))+simple_theme

#Save genes in each module
mod1_genes<-gene_module_df%>%dplyr::filter(module==1)
names(mod1_genes)[1]<-"gene_id"
mod1_genes_names<-merge(mod1_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod1_genes_names<-mod1_genes_names[,c(1,9,10,12,14)]
write.csv(mod1_genes_names, "/Users/mod1.csv")

#Use above code for module2-4.

#Regression analysis in whole population
table(cds$identifier)
cds$identifier<-factor(cds$identifier, levels =  c("No engraft","Engraft"))
table(cds$identifier)

gene_fits <- fit_models(cds, model_formula_str = "~identifier")
fit_coefs <- coefficient_table(gene_fits)
identity_terms <- fit_coefs %>% filter(term != "(Intercept)")

identity_DEG <- identity_terms %>% filter (q_value < 0.05) %>%
  dplyr::select(id, gene_short_name, num_cells_expressed,term, estimate, std_err, test_val, normalized_effect, q_value, p_value)

write.csv(identity_DEG,"/Users/Data.csv")


#Figure5E-G, Figure S8D and E, Figure S9A, D and G
#color
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
cell_type_color1 <- c("No engraft"="gray","Engraft" = "red")

#Dormancy score : Cell, Volume 169, Issue 5, 18 May 2017, Pages 807-823.e19
RES_DIR <- file.path("/Users")

estimate_score <- function(cds,gene_markers){
  cds_gene = cds[fData(cds)$gene_short_name %in% gene_markers,] 
  aggregate_gene = exprs(cds_gene)
  aggregate_gene = Matrix::t(Matrix::t(aggregate_gene) / pData(cds_gene)$Size_Factor)
  aggregate_gene = Matrix::colSums(aggregate_gene)
  pData(cds)$gene_score = log(aggregate_gene +1) 
  return(cds)
}

gene_markers <- c(read.csv(file.path(RES_DIR, "Dormant.csv"))) 
cds <- estimate_score(cds, gene_markers = gene_markers$Dormant)
plot_cells(cds, color_cells_by = 'gene_score', cell_size = 2, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

#violin plots
df <- data.frame(colData(cds)$identifier, colData(cds)$gene_score)
names(df) <- c("identifier", "gene_score")
write.csv(df, "~/gene_cluster.csv")

df <- data.frame(read.csv("~/gene_cluster.csv"))
p <- ggplot(df, aes(x= identifier, y=gene_score, fill = identifier)) + geom_violin(trim = FALSE)  + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + theme_bw()+ scale_fill_manual(values = cell_type_color1)
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0.5,0.5,-0.5,0.5), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=18))
p

#Wilcoxon test
p <- ggplot(df, aes(x= identifier, y=gene_score, fill = identifier)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.color=NA)  
p <- p + stat_compare_means() + stat_mean()  
p

#Other gene-sets were tested similarly.

#serial transplantaable score: Nature volume 583, pages585–589 (2020)
#Low output score: Nature volume 583, pages585–589 (2020)
#Durable self-renewal, MolO HSC score,Cell Stem Cell Volume 16, Issue 6, 4 June 2015, Pages 712-724
#Diapause score : Dev Cell. 2015 Nov 9; 35(3): 366–382., Cancer Discov. 2021 Jun;11(6):1542-1561. 
#REPOSIG : EMBO Reports (2022)23:e55502https://doi.org/10.15252/embr.202255502
#High output score: Nature volume 583, pages585–589 (2020)
#Multilineage score: Nature volume 583, pages585–589 (2020)
#Activated HSC/MMP score : Cell, Volume 169, Issue 5, 18 May 2017, Pages 807-823.e19, Cell Stem Cell, Volume 29, Issue 1, 6 January 2022, Pages 131-148.e10
#Mouse Gene Set: HALLMARK_MYC_TARGETS_V1,	MM3887, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MYC_TARGETS_V1.html
#Mouse Gene Set: HALLMARK_MYC_TARGETS_V2, MM2888, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MYC_TARGETS_V2.html
#WP_TCA_CYCLE, MM15856, WP434, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_TCA_CYCLE.html
#Mouse Gene Set: WP_PURINE_METABOLISM,	WP2185, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_PURINE_METABOLISM.html
#HALLMARK_OXIDATIVE_PHOSPHORYLATION, MM3893, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_OXIDATIVE_PHOSPHORYLATION.html,  http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_OXIDATIVE_PHOSPHORYLATION
#WP_MRNA_PROCESSING, MM15946, WP31, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_MRNA_PROCESSING.html
#REACTOME_TRANSLATION, MM15420, R-MMU-72766, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_TRANSLATION.html
#WP_CHEMOKINE_SIGNALING_PATHWAY, MM15943, WP2292, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_CHEMOKINE_SIGNALING_PATHWAY.html
#Commonly upregulated dormancy genes : Nat Cell Biol. 2024 Feb;26(2):181-193. doi: 10.1038/s41556-023-01325-3.
#Serially engrafting FL-HSC genes(Table S5)


#Figure S9F
#E9&E9.5 data (Cell Rep. 2021 Sep 14;36(11):109675. doi: 10.1016/j.celrep.2021.109675.)
#https://github.com/FredHutch/dignum-etal-2021/blob/main/TD%20Manuscript%20Code%20vs042021.R
#Designate the location
RES_DIR <-file.path("/Users")
cds <- readRDS(file.path(RES_DIR, "Sample1_2_P.RDS"))
cds<-choose_cells(cds) 
#Test our original gene set
RES_DIR <- file.path("/Users")

estimate_score <- function(cds,gene_markers){
  cds_gene = cds[fData(cds)$gene_short_name %in% gene_markers,] 
  aggregate_gene = exprs(cds_gene)
  aggregate_gene = Matrix::t(Matrix::t(aggregate_gene) / pData(cds_gene)$Size_Factor)
  aggregate_gene = Matrix::colSums(aggregate_gene)
  pData(cds)$gene_score = log(aggregate_gene +1) 
  return(cds)
}

gene_markers <- c(read.csv(file.path(RES_DIR, "C2GENES.csv"))) 
cds <- estimate_score(cds, gene_markers = gene_markers$C2GENES)
plot_cells(cds, color_cells_by = 'gene_score', cell_size = 2, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol)+simple_theme


#Code from "https://github.com/FredHutch/dignum-etal-2021/blob/main/TD%20Manuscript%20Code%20vs042021.R"
#violin
pData(cds)$cxcr4 = Matrix::colSums(exprs(cds[fData(cds)$gene_short_name %in% c("Cxcr4")]))
df <- data.frame(monocle3::clusters(cds), (colData(cds)$cxcr4>0))
names(df) <- c("cluster", "cxcr4")
df2 <- df %>% 
  group_by(cluster,cxcr4) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))
df3<-df2 %>% filter(cxcr4 == "TRUE") 
df4<-df3 %>% filter(perc>0.05)
colData(cds)$cxcr4 = (monocle3::clusters(cds)  %in% df4$cluster)
cxcr4_color <- c("TRUE"= "royalblue", "FALSE"= "gray")

plot_cells(cds, color_cells_by = "cxcr4", label_cell_groups = F, cell_size = 1, show_trajectory_graph = FALSE) +
  scale_color_manual(values = cxcr4_color) + simple_themeL

df <- data.frame(colData(cds)$cxcr4,  colData(cds)$gene_score) 
names(df) <- c("cxcr4",  "gene_score")
df$cxcr4<-factor(df$cxcr4, levels = c("TRUE", "FALSE"))
p <- ggplot(df, aes(x= cxcr4, y=gene_score, fill = cxcr4)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + stat_compare_means() + stat_mean()   
p <- p + theme_bw() + scale_fill_manual(values = cxcr4_color)
p <- p + xlab(NULL) + ylab(NULL) 
p <- p + theme(plot.margin = unit(c(0.5,0.5,-0.5,0.5), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=15))
p

#Figure S9H
#Esam expression, violin
cell_type_color1 <- c("No engraft"="gray","Engraft" = "red")

hsc_genes<-c("Esam")
cds_hsc_genes<-cds[rowData(cds)$gene_short_name %in% hsc_genes,]
gene_fits<- fit_models(cds_hsc_genes, model_formula_str="~identifier")
fit_coefs<-coefficient_table(gene_fits) 
hsc_genes_terms <-fit_coefs %>% filter(term !="(Intercept)")
hsc_genes_terms %>% filter (q_value<0.05) %>% dplyr::select(gene_short_name, term, q_value, estimate,p_value)

#plot by violin
plot_genes_violin(cds_hsc_genes, group_cells_by= "identifier") 
plot_genes_violin(cds_hsc_genes, group_cells_by= "identifier") + scale_fill_manual(values = cell_type_color1) + xlab(NULL) + ylab(NULL) + theme(plot.margin = unit(c(0.5,0.5,-0.5,0.5), "cm"), text = element_text(size=18))  


#End
