#Figure S11

#read cds
cds1 <- load_cellranger_data("-")
cds2 <- load_cellranger_data("-")

cds1$identifier<-("HD1-HSC")
cds2$identifier<-("HD2-HSC")

#cds1
run_scrublet <- function(cds1) {
  
  require(reticulate)
  use_python("/usr/local/bin/python3")
  scr <- import('scrublet')
  plt <- import('matplotlib.pyplot') #IP added
  
  count.matrix <- r_to_py(t(exprs(cds1)),convert=TRUE)
  scrub <- scr$Scrublet(count.matrix)
  
  out <- scrub$scrub_doublets()
  
  pData(cds1)$scrub_kNN <- out[[1]]
  pData(cds1)$scrub_call <- out[[2]]
  cds1
}

cds1 <- run_scrublet(cds1)
colData(cds1)
table(cds1$scrub_call)
cds1 <- cds1[,(colData(cds1)$scrub_call == "FALSE"),]
table(cds1$scrub_call)

#cds2
run_scrublet <- function(cds2) {
  
  require(reticulate)
  use_python("/usr/local/bin/python3")
  scr <- import('scrublet')
  plt <- import('matplotlib.pyplot') #IP added
  
  count.matrix <- r_to_py(t(exprs(cds2)),convert=TRUE)
  scrub <- scr$Scrublet(count.matrix)
  
  out <- scrub$scrub_doublets()
  
  pData(cds2)$scrub_kNN <- out[[1]]
  pData(cds2)$scrub_call <- out[[2]]
  cds2
}

cds2 <- run_scrublet(cds2)
colData(cds2)
table(cds2$scrub_call)
cds2 <- cds2[,(colData(cds2)$scrub_call == "FALSE"),]
table(cds2$scrub_call)

#combine
cds <- combine_cds(list(cds1,cds2), keep_all_genes=TRUE)
# Remove low UMIs 
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
qplot(colData(cds)$n.umis, geom="density")
qplot(log10(colData(cds)$n.umis), geom="density")
summary(colData(cds)$n.umis)
cds <- cds[,Matrix::colSums(counts(cds)) > 3000]

#Remove cells with low genes per cell
cds <- detect_genes(cds, min_expr=0.1)     
cds <- cds[,colData(cds)$num_genes_expressed > 1000]
summary(colData(cds)$n.umis)
summary(pData(cds)$num_genes_expressed)

### plot UMI and genes per cell for each sample
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
                        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
                        panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_text(size=10))

df1 <- data.frame("identifier"=colData(cds)$identifier, "n.umis"=colData(cds)$n.umis)

df1$identifier <- factor(df1$identifier, levels = c("HD1-HSC","HD2-HSC"))

ggplot(df1, aes(x=identifier, y=n.umis, fill=identifier)) +geom_boxplot() +  
  scale_y_continuous(trans='log10', limits=c(100,200000), labels=number) + simple_theme2

df2 <- data.frame("identifier"=colData(cds)$identifier, "num_genes_expressed"=colData(cds)$num_genes_expressed)

df2$identifier <- factor(df2$identifier, levels = c("HD1-HSC","HD2-HSC"))

ggplot(df2, aes(x=identifier, y=num_genes_expressed, fill=identifier)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 

##pre-process 
cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)
#Align : remove batch effects 
cds = align_cds(cds, alignment_group = "identifier")

##Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

set.seed(1)
cds <- cluster_cells(cds, resolution=0.02,random_seed=1)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size = 1)

#Exclude cells with high mitochondrial reads
mito_genes <- rowData(cds)
mito_genes <- mito_genes[str_detect(mito_genes$gene_short_name, "MT-"),]
mito_genes_id <- mito_genes$id
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
num_mito <- Matrix::colSums(counts(cds[mito_genes_id]))
cds$n.mito <- num_mito

#calculate the percent of mitochondrial UMIs
perc_mito <- 100*(cds$n.mito / cds$n.umis)
cds$perc_mito_umis <- perc_mito
cds_mito_filter <- cds[,pData(cds)$perc_mito_umis < 10]
cds_mito_filter
cds<- cds_mito_filter
cds

#plot 
plot_cells(cds, genes=c("PTPRC"), cell_size = 1, label_cell_groups = F)
plot_cells(cds, genes=c("CDH5"), cell_size = 1, label_cell_groups = F)

#Exclude CDH5+PTPRC- 
cds<-choose_cells(cds)

#save cds
RES_DIR <- file.path("-")
saveRDS(cds, file.path(RES_DIR, "-")) 

#clustering
set.seed(1)
cds <- cluster_cells(cds, resolution=0.05, random_seed = 1)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size=1, group_label_size = 3, label_cell_groups = F)
plot_cells(cds, color_cells_by = "cluster", cell_size=1, label_cell_groups=T, group_label_size=10) + simple_theme

#add cluster to colData
cluster <- clusters(cds)
pData(cds)$cluster <- cluster  
table(pData(cds)$cluster)

#count
table(pData(cds)$identifier)
#A
cluster4 <- cds[,(colData(cds)$cluster == "4"),]
table(pData(cluster4)$identifier)
#B
cluster1 <- cds[,(colData(cds)$cluster == "1"),] 
table(pData(cluster1)$identifier)
#C
cluster2 <- cds[,(colData(cds)$cluster == "2"),]
table(pData(cluster2)$identifier)
#D
cluster3 <- cds[,(colData(cds)$cluster == "3"),]
table(pData(cluster3)$identifier)

#change the number of cluster into alphabetical order 
colData(cds)
colData(cds)$group = as.character(clusters(cds))
colData(cds)

colData(cds)$group = dplyr::recode(colData(cds)$group,
                                   "1"= "B",
                                   "2"= "C",
                                   "3"= "D",
                                   "4" = "A")
colData(cds)
table(pData(cds)$group)

plot_cells(cds, color_cells_by = "group", cell_size=2, label_cell_groups=T, group_label_size=20) + simple_theme

#Gene-set scores
#color 
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

##HSC_Cycle_Primed signature:Garcia-Prat et al Cell Stem Cell 2021.
#Cell Stem Cell. 2021 Oct 7;28(10):1838-1850.e10. doi: 10.1016/j.stem.2021.07.003.
RES_DIR <- file.path("-")
estimate_score <- function(cds,gene_markers){
  cds_gene = cds[fData(cds)$gene_short_name %in% gene_markers,]
  aggregate_gene = exprs(cds_gene)
  aggregate_gene = Matrix::t(Matrix::t(aggregate_gene) / pData(cds_gene)$Size_Factor)
  aggregate_gene = Matrix::colSums(aggregate_gene)
  pData(cds)$HSCPRIMED = log(aggregate_gene +1)
  return(cds)
}

gene_markers <- c(read.csv(file.path(RES_DIR, "-.csv")))
cds <- estimate_score(cds, gene_markers = gene_markers$HSCPRIMED)
plot_cells(cds, color_cells_by = 'HSCPRIMED', cell_size = 1, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

#Checked each gene-set score

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
mat<-avg_score_matrix(cds, "group", c("HSCPRIMED","Activated","MYCV1","MYCV2"
                                      ,"Quiescent","HSCSIGNATURE","LIPIDMETABOLISM","WEIJDEN"))


mat

pheatmap(mat[c("HSCPRIMED","Activated","MYCV1","MYCV2"
               ,"Quiescent","HSCSIGNATURE","LIPIDMETABOLISM","WEIJDEN")
             ,c("A","B","C","D")]
         , legend=F, show_rownames = F, show_colnames = F, cluster_cols=F, cluster_rows = F
         , cellheight = 40
         , cellwidth = 75) 

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

pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   show_rownames = F, show_colnames = F, legend=F, cellheight=45, cellwidth = 50)

#### plot gene modules in UMAP 
plot_cells(cds, cell_size = 1,
           genes=gene_module_df %>% filter(module %in% c(1,2,3,4,5,6,7,8)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           scale_to_range = TRUE,
           min_expr = 60) +scale_color_gradientn(colors=c("gray80", "darkgreen", "darkgreen"))+simple_theme

#Save genes in each module
mod1_genes<-gene_module_df%>%dplyr::filter(module==1)
names(mod1_genes)[1]<-"gene_id"
mod1_genes_names<-merge(mod1_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod1_genes_names<-mod1_genes_names[,c(1,9,10,11,12)]
write.csv(mod1_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod1.csv")

mod2_genes<-gene_module_df%>%dplyr::filter(module==2)
names(mod2_genes)[1]<-"gene_id"
mod2_genes_names<-merge(mod2_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod2_genes_names<-mod2_genes_names[,c(1,9,10,11,12)]
write.csv(mod2_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod2.csv")

mod3_genes<-gene_module_df%>%dplyr::filter(module==3)
names(mod3_genes)[1]<-"gene_id"
mod3_genes_names<-merge(mod3_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod3_genes_names<-mod3_genes_names[,c(1,9,10,11,12)]
write.csv(mod3_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod3.csv")

mod4_genes<-gene_module_df%>%dplyr::filter(module==4)
names(mod4_genes)[1]<-"gene_id"
mod4_genes_names<-merge(mod4_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod4_genes_names<-mod4_genes_names[,c(1,9,10,11,12)]
write.csv(mod4_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod4.csv")

mod5_genes<-gene_module_df%>%dplyr::filter(module==5)
names(mod5_genes)[1]<-"gene_id"
mod5_genes_names<-merge(mod5_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod5_genes_names<-mod5_genes_names[,c(1,9,10,11,12)]
write.csv(mod5_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod5.csv")

mod6_genes<-gene_module_df%>%dplyr::filter(module==6)
names(mod6_genes)[1]<-"gene_id"
mod6_genes_names<-merge(mod6_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod6_genes_names<-mod6_genes_names[,c(1,9,10,11,12)]
write.csv(mod6_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod6.csv")

mod7_genes<-gene_module_df%>%dplyr::filter(module==7)
names(mod7_genes)[1]<-"gene_id"
mod7_genes_names<-merge(mod7_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod7_genes_names<-mod7_genes_names[,c(1,9,10,11,12)]
write.csv(mod7_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod7.csv")

mod8_genes<-gene_module_df%>%dplyr::filter(module==8)
names(mod8_genes)[1]<-"gene_id"
mod8_genes_names<-merge(mod8_genes,pr_graph_test_res,by.x="gene_id",by.y="id")
mod8_genes_names<-mod8_genes_names[,c(1,9,10,11,12)]
write.csv(mod8_genes_names, "/Users/ti/Desktop/human scRNAseq/Aldinger-T21/HD-FL only analysis, 122424/onlyHSC/post/mod8.csv")