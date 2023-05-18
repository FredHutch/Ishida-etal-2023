#Figure7A-C
writeLines(capture.output(sessionInfo()), "sessionInfoNichenetr.txt")

#Load FL-AKT-EC
RES_DIR <-file.path("/Users")

cds<-load_cellranger_data("/Users/Data")

#Remove low UMI
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
summary(colData(cds)$n.umis)
qplot(colData(cds)$n.umis, geom="density")
qplot(log10(colData(cds)$n.umis), geom="density")
cds <- cds[,Matrix::colSums(counts(cds)) > 10000]

cds
summary(colData(cds)$n.umis)

cds <- detect_genes(cds, min_expr=0.1)     
cds <- cds[,colData(cds)$num_genes_expressed > 1000]
summary(colData(cds)$n.umis)
summary(pData(cds)$num_genes_expressed)


#Preprocess
cds <- preprocess_cds(cds)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)
cluster <- clusters(cds)
plot_cells(cds)


#List of all mitochondrial genes by Ensembl ID
mito_genes_updated <- c("ENSMUSG00000064351", "ENSMUSG00000064354", "ENSMUSG00000064370", "ENSMUSG00000064357")

num_mito <- Matrix::colSums(counts(cds[mito_genes_updated]))
cds$n.mito <- num_mito
cds

colData(cds)$n.umis <- Matrix::colSums(counts(cds))
perc_mito <- 100*(cds$n.mito / cds$n.umis)
cds$perc_mito_umis <- perc_mito

cds_mito_filter <- cds[,pData(cds)$perc_mito_umis < 10]

cds_mito_filter
cds<- cds_mito_filter

set.seed(1)
cds <- cluster_cells(cds, resolution=0.00001, random_seed = 1)
cluster <- clusters(cds)
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", cell_size = 0.5, group_label_size=5)

#Save
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))
saveRDS(cds, file.path(RES_DIR, "Figure7 FLAKTEC.RDS")) 



#Receptor-Ligand interaction analysis by Nichenetr
simple_theme <-  theme(text = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), legend.position = "none")  ### theme to remove legends and axis/font

simple_themeL <-  theme(text = element_text(size=15),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank())   #### theme show the legend with large font


library(nichenetr)
library(tidyverse)
library(patchwork)

#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
RES_DIR <- file.path("/Users")
ligand_target_matrix <- readRDS(file.path(RES_DIR, "ligand_target_matrix.rds"))

ligand_target_matrix[1:5,1:5]
dim(ligand_target_matrix)

#Convert the ligand-target model from human to mouse symbols.
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)

#Perform the NicheNet analysis
#Define expressed genes in sender and receiver cell populations
#receiver/target = Engraft FL-HSC colony (Colony#2)

#location
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))
#Read
cds <- readRDS(file.path(RES_DIR, "Figure4.RDS"))
colData(cds)
cds2 <- cds[,(colData(cds)$identifier == "Engraft"),] 
colData(cds2)

cds2 <- detect_genes(cds2, min_expr = 0.1)   
rowData(cds2)$use_for_rec_list <- rowData(cds2)$num_cells_expressed > (0.10 * ncol(cds2)) ### detect genes expressed in >10% of cells
cds2_exp <- (subset(rowData(cds2), rowData(cds2)$use_for_rec_list == TRUE))$gene_short_name  ###create list of expressed genes
length(cds2_exp)

#sender/niche=FL-AKT-EC
#location
RES_DIR <- file.path("/Users")
#Read
cds3 <- readRDS(file.path(RES_DIR, "Figure7 FLAKTEC.RDS"))
cds3 <- detect_genes(cds3, min_expr = 0.1)    
rowData(cds3)$use_for_rec_list <- rowData(cds3)$num_cells_expressed > (0.10 * ncol(cds3)) ### detect genes expressed in >10% of cells
cds3_exp <- (subset(rowData(cds3), rowData(cds3)$use_for_rec_list == TRUE))$gene_short_name  ###create list of genes
length(cds3_exp)

#Define the gene set of interest and a background of gene
RES_DIR <- file.path("/Users")

geneset_oi = read.csv(file.path(RES_DIR, "Table.csv"))%>% pull(gene_short_name)

#only consider genes also present in the NicheNet model
geneset_oi =geneset_oi %>% .[. %in% rownames(ligand_target_matrix)] 

head(geneset_oi)
length(geneset_oi)

#background of genes = All the HSC colonies (Colony#1-6)
cds <- detect_genes(cds, min_expr = 0.1)   
rowData(cds)$use_for_rec_list <- rowData(cds)$num_cells_expressed > (0.10 * ncol(cds)) 
cds_exp <- (subset(rowData(cds), rowData(cds)$use_for_rec_list == TRUE))$gene_short_name  
length(cds_exp)
background_expressed_genes = cds_exp  %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)


#Define a set of potential ligands
#Putative ligand-receptor links were gathered from NicheNet’s ligand-receptor data sources.
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

#for mouse, change the gene name
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
lr_network

#Ligand
ligands = lr_network %>% pull(from) %>% unique()
#Expressed Ligands in FL-AKT-EC
expressed_ligands = intersect(ligands,cds3_exp)
length(expressed_ligands)

#Receptor
receptors = lr_network %>% pull(to) %>% unique()

#Expressed Ligands in Engraft colony
expressed_receptors = intersect(receptors,cds2_exp)
length(expressed_receptors)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

#Perform NicheNet’s ligand activity analysis on the gene set of interest
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

#Ranking of ligand activities by pearson
ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="skyblue")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity2 (Engraft-colony)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands%>% intersect(expressed_ligands)

#Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()%>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","DE genes highly expressed in Engraft-colony", color = "#0457d2",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "#0457d2", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network + simple_themeL

#Ligand-receptor network inference for top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

#get the weights of the ligand-receptor interactions as used in the NicheNet model
#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds")) 
RES_DIR <- file.path("/Users")
weighted_networks <- readRDS(file.path(RES_DIR, "weighted_networks.rds"))

#change to mouse
weighted_networks$lr_sig =  weighted_networks$lr_sig %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks$gr =  weighted_networks$gr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

head(weighted_networks$lr_sig) 
head(weighted_networks$gr) 

lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized ligands","Receptors expressed in Engraft-colonies", color = "mediumvioletred",  x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network+simple_themeL

#Visualize top-predicted ligands with receptors
library(RColorBrewer)
library(cowplot)
library(ggpubr)

ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Ligands")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkred",x_axis = FALSE, y_axis = TRUE, legend_position = "left",  legend_title = "Pearson")
p_ligand_pearson+simple_themeL


#End