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
library(ComplexHeatmap)


simple_theme <-  theme(text = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), legend.position = "none")  ### theme to remove legends and axis/font


simple_themeT <-  theme(text = element_text(size=25),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank(),legend.position="none")


writeLines(capture.output(sessionInfo()), "sessionInfoFreshFLHSC.txt")


#set resident directory 
RES_DIR <-file.path("/Users")

#load cds 
cds<-load_cellranger_data("/Users")


colData(cds)$n.umis <- Matrix::colSums(counts(cds))
summary(colData(cds)$n.umis)
qplot(colData(cds)$n.umis, geom="density")
qplot(log10(colData(cds)$n.umis), geom="density")
cds <- cds[,Matrix::colSums(counts(cds)) > 10000]

####identify expressed genes and remove cells with low genes per cell
cds <- detect_genes(cds, min_expr=0.1)     
cds <- cds[,colData(cds)$num_genes_expressed > 1000]
summary(colData(cds)$n.umis)
summary(pData(cds)$num_genes_expressed)

## Reduce dimension
cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

## Cluster the cells
cds <- cluster_cells(cds)
cluster <- clusters(cds)
plot_cells(cds)

#Remove cells with high content of mitochondrial genes
mito_genes_updated <- c("ENSMUSG00000064351", "ENSMUSG00000064354", "ENSMUSG00000064370", "ENSMUSG00000064357")

#Now calculate the number of umis associated with mitochondrial reads for each cell
num_mito <- Matrix::colSums(counts(cds[mito_genes_updated]))
cds$n.mito <- num_mito

#calculate n.umis
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
#calculate percent mito umis
perc_mito <- 100*(cds$n.mito / cds$n.umis)
cds$perc_mito_umis <- perc_mito

cds_mito_filter <- cds[,pData(cds)$perc_mito_umis < 5]

cds<- cds_mito_filter

#Designate the location
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))

#Save
#whole
saveRDS(cds, file.path(RES_DIR, "Figure3.RDS")) 
#Read
cds <- readRDS(file.path(RES_DIR, "Figure3.RDS"))

#Figure S2B-D
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
plot_cells(cds, genes=c("C1qb"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Elane"), cell_size = 4, label_cell_groups = F)+ simple_themeT

plot_cells(cds, genes=c("Cd48"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Flt3"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Il7r"),cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Cebpe"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Gypa"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Ms4a2"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Siglec1"), cell_size = 4, label_cell_groups = F)+ simple_themeT

plot_cells(cds, genes=c("Mki67"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Cdk1"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Top2a"), cell_size = 4, label_cell_groups = F)+ simple_themeT
plot_cells(cds, genes=c("Cdc20"), cell_size = 4, label_cell_groups = F)+ simple_themeT


### Estimate proliferation index, Fig.S2D (Cell Rep. 2021 Sep 14;36(11):109675. doi: 10.1016/j.celrep.2021.109675.)
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


##Figure3B-D, FigureS4E
#cluster cells
set.seed(1)
cds <- cluster_cells(cds, resolution=0.2, random_seed = 1)
cluster <- clusters(cds)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size = 2, group_label_size=20, label_cell_groups = T)+simple_theme

#add cluster to colData
cluster <- clusters(cds)
pData(cds)$cluster <- cluster  
table(pData(cds)$cluster)

#change the number of cluster into alphabetical order from the left on UMAP
colData(cds)
colData(cds)$group = as.character(clusters(cds))
colData(cds)


colData(cds)$group = dplyr::recode(colData(cds)$group,
                                   "1"= "K",
                                   "2"= "F",
                                   "3"= "H",
                                   "4" = "J",
                                   "5" = "A",
                                   "6" = "G",
                                   "7" = "D",
                                   "8" = "L",
                                   "9" = "B",
                                   "10" = "E",
                                   "11" = "C",
                                   "12" = "I")
colData(cds)
table(pData(cds)$group)

plot_cells(cds, color_cells_by = "group", cell_size=2, label_cell_groups=T, group_label_size=20) + simple_theme


#Gene-set scores
#color 
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#Dormancy score : Cell, Volume 169, Issue 5, 18 May 2017, Pages 807-823.e19
RES_DIR <- file.path("/Users")

estimate_score <- function(cds,gene_markers){
  cds_gene = cds[fData(cds)$gene_short_name %in% gene_markers,] 
  aggregate_gene = exprs(cds_gene)
  aggregate_gene = Matrix::t(Matrix::t(aggregate_gene) / pData(cds_gene)$Size_Factor)
  aggregate_gene = Matrix::colSums(aggregate_gene)
  pData(cds)$Dormancy = log(aggregate_gene +1)
  return(cds)
}

gene_markers <- c(read.csv(file.path(RES_DIR, "Dormant.csv")))
cds <- estimate_score(cds, gene_markers = gene_markers$Dormant)
plot_cells(cds, color_cells_by = 'Dormancy', cell_size = 2, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol)+simple_theme

### signature gene sets:
#serial transplantaable score: Nature volume 583, pages585–589 (2020)
#Diapause score : Dev Cell. 2015 Nov 9; 35(3): 366–382., Cancer Discov. 2021 Jun;11(6):1542-1561. 
#WP_CHEMOKINE_SIGNALING_PATHWAY, MM15943, WP2292, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_CHEMOKINE_SIGNALING_PATHWAY.html
#High output score: Nature volume 583, pages585–589 (2020)
#Multilineage score: Nature volume 583, pages585–589 (2020)
#Activated HSC/MMP score: : Cell, Volume 169, Issue 5, 18 May 2017, Pages 807-823.e19, Cell Stem Cell, Volume 29, Issue 1, 6 January 2022, Pages 131-148.e10
#Mouse Gene Set: HALLMARK_MYC_TARGETS_V1, 	MM3887, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MYC_TARGETS_V1.html
#Mouse Gene Set: HALLMARK_MYC_TARGETS_V2, MM2888, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MYC_TARGETS_V2.html
#WP_TCA_CYCLE, MM15856, WP434, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_TCA_CYCLE.html
#HALLMARK_OXIDATIVE_PHOSPHORYLATION, MM3893, http://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_OXIDATIVE_PHOSPHORYLATION.html,  http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_OXIDATIVE_PHOSPHORYLATION
#Mouse Gene Set: WP_PURINE_METABOLISM,	WP2185, https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/WP_PURINE_METABOLISM.html
#Serially engrafting FL-HSC genes(Table S5)

#Heatmap
avg_score_matrix<-function(cds, cell_group, scores){
  avg_score_list<-lapply(levels(factor(cds[[cell_group]])), function(x){
    sub<-cds@colData[cds[[cell_group]] == x, scores] %>% as.data.frame()
    avg<-colMeans(sub)
  })
  names(avg_score_list)<-levels(factor(cds[[cell_group]]))
  mat<-do.call(rbind, avg_score_list)
  mat<-t(scale(mat))
  mat
}

names(cds@colData)

##get scaled matrix
mat<-avg_score_matrix(cds, "group", c("Dormancy", "SerialTransplant","Diapause"
                                      ,"CHEMOKINE"
                                      ,"HighOutput", "Multilineage","ActivatedHSCMPP"
                                      ,"MYCV1","MYCV2"
                                      ,"TCACYCLE","OXOPHOS","PurineMetabolism"))

#Look at mat
mat

#make heatmap
pheatmap(mat[c("Dormancy", "SerialTransplant","Diapause"
               ,"CHEMOKINE"
               ,"HighOutput", "Multilineage","ActivatedHSCMPP"
               ,"MYCV1","MYCV2"
               ,"TCACYCLE","OXOPHOS","PurineMetabolism")
             ,c("A","B","C","D","E","F","G","H","I","J","K","L")]
         , legend=F, show_rownames = F, show_colnames = F, cluster_cols=F, cluster_rows = F
         , cellheight = 30
         , cellwidth = 45) 

#Heatmap Figure for Engraft FL-HSC colony genes
avg_score_heatmap<-function(cds, cell_group, scores){
  avg_score_list<-lapply(levels(factor(cds[[cell_group]])), function(x){
    sub<-cds@colData[cds[[cell_group]] == x, scores] %>% as.data.frame()
    avg<-colMeans(sub)
  })
  names(avg_score_list)<-levels(factor(cds[[cell_group]]))
  mat<-do.call(rbind, avg_score_list)
  mat<-t(scale(mat)) 
  pheatmap(mat, cluster_cols=F, legend=F, show_rownames = F, show_colnames = F,
           cellheight = 35,  
           cellwidth = 45)  
}

names(cds@colData)

avg_score_heatmap(cds, "group", 
                  c("GENES"))

#Figure3E, FigureS2E
#Gene module analysis
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
set.seed(1) 
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.001, random_seed = 1)


#each group
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$group)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("group", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   show_rownames = F, show_colnames = F, legend=F, cellheight=35, cellwidth = 40)

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12, show_rownames = T, show_colnames = T, legend=T, cellheight=20, cellwidth = 30)

#### plot gene modules in UMAP 
plot_cells(cds, cell_size = 2,
           genes=gene_module_df %>% filter(module %in% c(1,3,4,7)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           scale_to_range = TRUE,
           min_expr = 50) +scale_color_gradientn(colors=c("gray80", "darkgreen", "darkgreen"))+simple_theme

#Save genes in each module
mod1_genes<-gene_module_df%>%dplyr::filter(module==1)
names(mod1_genes)[1]<-"gene_id"
mod1_genes_names<-merge(mod1_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod1_genes_names<-mod1_genes_names[,c(1,9,10,11,12)]
write.csv(mod1_genes_names, "/Users/mod1.csv")

#use the code for module2-11

#End#
