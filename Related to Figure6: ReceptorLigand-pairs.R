#Search for the ligand-receptor pairs
#Nature Communications volume 6, Article number: 7866 (2015)
#Nature Communications volume 13, Article number: 1584 (2022)

###define file path to output folders
RES_DIR <- file.path("/Users")

RecLigAll = read.csv(file.path(RES_DIR, "RecLigPairs.csv"))  # read csv file of all Receptor-Ligand pairs in database (mouse gene names)

#Receptor : Read HSC colonies
RES_DIR <- file.path(RES_DIR <- file.path("/Users"))
#Read
cds <- readRDS(file.path(RES_DIR, "Figure5.RDS"))

cds1 <- cds[,colData(cds)$identifier == "Engraft"]   
cds1 <- detect_genes(cds1, min_expr = 0.1)  
rowData(cds1)$use_for_rec_list <- rowData(cds1)$num_cells_expressed > (0.10 * ncol(cds1)) 
cds1_exp <- (subset(rowData(cds1), rowData(cds1)$use_for_rec_list == TRUE))$gene_short_name  
cds1_receptors <- intersect(RecLigAll$Receptor, cds1_exp)    

RES_DIR <- file.path("/Users")
write.csv(cds1_receptors, file.path(RES_DIR, "cds1_receptors.csv"))         

#Ligand : FL-AKT-EC
RES_DIR <- file.path("/Users")
#Read
cds2 <- readRDS(file.path(RES_DIR, "Figure6FLAKTEC.RDS"))


cds2 <- detect_genes(cds2, min_expr = 0.1)    
rowData(cds2)$use_for_rec_list <- rowData(cds2)$num_cells_expressed > (0.10 * ncol(cds2)) 
cds2_exp <- (subset(rowData(cds2), rowData(cds2)$use_for_rec_list == TRUE))$gene_short_name 
cds2_ligands <- intersect(RecLigAll$Ligand, cds2_exp)        

RES_DIR <- file.path("/Users")
write.csv(cds2_ligands, file.path(RES_DIR, "cds2_ligands.csv"))          

######  Identify R_L pairs
RES_DIR <- file.path("/Users")
cds1_receptors = read.csv(file.path(RES_DIR, "cds1_receptors.csv"))          
cds2_ligands = read.csv(file.path(RES_DIR, "cds2_ligands.csv"))         
RecLigPairs <- RecLigAll[(RecLigAll$Ligand %in% cds2_ligands$x & RecLigAll$Receptor %in% cds1_receptors$x),] 
write.csv(RecLigPairs, file.path(RES_DIR, "cds1_cds2_RecLigPairs.csv"))  

