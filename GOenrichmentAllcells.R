library(tidyr)
library(gdata) # requires ggplot2
library(ggrepel)

dge=read.csv('/Users/mprokhorenko/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine(2)/tbi/dentate/AllCellsdegswlabels2.csv')
sub0 = dge[dge$FDR<=0.05,]
Age = sub0[sub0$Contrast=="bothgroups.3movs18mo", ]
Age$negLogP = -log10(Age$PValue)
Age$posLogP = log10(Age$PValue)
TBI = sub0[sub0$Contrast=="TBIvssham.bothages", ]
TBI$negLogP = -log10(TBI$PValue)
TBI$posLogP = log10(TBI$PValue)
sub1 = rbind(TBI, Age)
sub1 = sub1[sub1$CellType!="Mixed",]
unique(sub1$CellType)


c("Micro" = "#661100", "Ependyma" = "#CC6677", 
  "CP" = "#AA4466", "CajalRhetzius" = "#882255",
  "EndoMural" = "#332288", "L3EC" = "#6C37C7", 
  "Nxph3" = "#51367A", "DG" = "#AA4499", 
  "CA1" = "#6699CC", "CA2CA3" = "#88CCEE", 
  "Astro" = "#44AA99", "Oligo" = "#47910E", 
  "OPC" = "#117733", "Lamp5" = "#045C4A", 
  "Sst" = "#999933", "Vip" = "#DDCC77")
my_cols2 = c("#44AA99","#332288", "#045C4A", "#999933", "#DDCC77",
             "#661100", "#6699CC", "#88CCEE", "#882255", "#AA4499", 
             "#51367A", "#6C37C7","#47910E", "#117733")


plt1 = ggplot(TBI) +
  aes(x = logFC, y = negLogP, colour = CellType) +
  geom_point() + 
  scale_color_manual(values= my_cols2) +
  theme_classic()


sub0 = dge[dge$CellType!="Mixed",]
Age = sub0[sub0$Contrast=="bothgroups.3movs18mo", ]
Age$negLogP = -log10(Age$PValue)
TBI = sub0[sub0$Contrast=="TBIvssham.bothages", ]
TBI$negLogP = -log10(TBI$PValue)
sub1 = rbind(TBI, Age)
length(unique(Age$Gene)) #20360
length(unique(TBI$Gene)) #20360
sharedDEGs = match(Age$gene, TBI$gene)
#sharedDEGs2 = match(sub2$Gene, sub1$Gene) #check that matching includes all
sub1$label = sub1$Gene
for (i in 1:length(sub1$label)){
  if (abs(sub1$logFC[i])<1){
    sub1$label[i] <- NA}
}
for (i in 1:length(sub1$label)){
  if (sub1$negLogP[i]<7){
    sub1$label[i] <- NA}
}
sub1$label <- sub("Rik", NA, sub1$label)
sub1$label[startsWith(sub1$label, 'Gm', trim=TRUE)] <- NA 
sub1$label[startsWith(sub1$label, 'AC', trim=TRUE)] <- NA 
write.table(sub1, sep = ",", "AllCellsdegswlabels3.csv")
temp <- read.table("AllCellsdegswlabels2.csv", header= TRUE, sep=",", row.names = 1)
temp = sub1
temp <- temp[order(temp$negLogP),]
plt1 = ggplot(temp) +
  aes(x = logFC, y = negLogP, colour = CellType) + #size = FDR
  geom_point() +
  scale_color_manual(values= my_cols2) +
  theme_classic() +
  facet_wrap(vars(Contrast)) +
  xlim(-4,6) + ylim(2,107) 

plt2 = plt1 + 
  geom_text_repel(inherit.aes = FALSE, 
                  aes(logFC, y = negLogP , label = label, size = 6),
                  max.overlaps = 18, segment.color = 'transparent') +
  theme_classic() +
  geom_vline(xintercept = 0, colour = "darkgray") 
png("allcellsDGE_volcanos27.png", width = 14, height = 5, units = 'in', res = 300)
plt2
dev.off()

plt1 = ggplot(temp) +
  aes(x = logFC, y = negLogP, colour = CellType) + #size = FDR
  geom_point() +
  scale_color_manual(values= my_cols2) +
  theme_classic() +
  facet_wrap(vars(Contrast)) +
  xlim(-4,6) + ylim(2,30) 
plt3 = plt1 + 
  geom_text_repel(inherit.aes = FALSE, 
                  aes(logFC, y = negLogP , label = label, size = 6),
                  max.overlaps = 10, segment.color = 'transparent') +
  theme_classic() +
  geom_vline(xintercept = 0, colour = "darkgray")
png("allcellsDGE_volcanos26.png", width = 14, height = 5, units = 'in', res = 300)
plt3
dev.off()

# prep data for GO analysis 
#dge=read.table('/Users/mprokhorenko/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine(2)/tbi/dentate/AllCellsdegswlabels2.csv',header = T,sep = '\t')
dge=read.csv('/Users/mprokhorenko/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine(2)/tbi/dentate/AllCellsdegswlabels2.csv')
cells=unique(dge$CellType)
ups = dge[dge$PValue<=0.05 & dge$logFC>0,]
downs = dge[dge$PValue<=0.05 & dge$logFC<0,]
degset = list()
cntrst = unique(dge$Contrast)
for (i in cntrst) {
  sub0 = ups[ups$Contrast==i, ]
  for (j in unique(dge$CellType)) {
    genes=as.character(sub0[sub0$CellType==j,8])
    degset[[paste("ups",i,j,sep='_')]]=genes
  }
}
for (i in cntrst) {
  sub0 = downs[downs$Contrast==i, ]
  for (j in unique(dge$CellType)) {
    genes=as.character(sub0[sub0$CellType==j,8])
    degset[[paste("downs",i,j,sep='_')]]=genes
  }
}
#L3EC done separately
cells="Neuron_Subiculum_Slc17a6" 
ups = dge[dge$PValue<=0.05 & dge$logFC>0,]
downs = dge[dge$PValue<=0.05 & dge$logFC<0,]
degset = list()
cntrst = unique(dge$Contrast)
for (i in cntrst) {
  sub0 = ups[ups$Contrast==i, ]
  genes=as.character(sub0[sub0$CellType==cells,8])
  degset[[paste("ups",i,cells,sep='_')]]=genes
}
for (i in cntrst) {
  sub0 = downs[downs$Contrast==i, ]
  genes=as.character(sub0[sub0$CellType==cells,8])
  degset[[paste("downs",i,cells,sep='_')]]=genes
}



#### go enrichment 
### 
library(enrichR)
options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
dbs <- listEnrichrDbs()

dbs[grep('KEGG',dbs$libraryName),]
dbs[grep('GO_Biological_Process',dbs$libraryName),]

cutoff=1
mydbs <- c("GO_Biological_Process_2023","KEGG_2019_Mouse")
golist2 = list()
for (i in 1:length(degset)){
  genes=degset[[i]]
  enriched <- enrichr(genes, mydbs)
  
  out=enriched$GO_Biological_Process_2023
  out=out[out$Adjusted.P.value<cutoff,]
  out$dbs = 'GO_BP_23'
  
  out2=enriched$KEGG_2019_Mouse
  out2=out2[out2$Adjusted.P.value<cutoff,]
  out2$dbs = 'KEGG_19'
  
  outs = rbind(out, out2)
  outs$cellcntrst = names(degset[i])
  golist2[[i]] = outs
}
golist2 = do.call( rbind , golist2)
golist1 = do.call( rbind , golist1)
golist2 = rbind(golist1, golist2)
#write.table(out,file=paste(i,'.enrichedKEGG.txt',sep=''),sep='\t',quote = F)
write.table(golist2,file=paste('tbi_pseudobulkDGEallcells_enriched.txt',sep=''),sep='\t',quote = F)

#######
#only use the GO Bio Processses
golist2 = golist2[golist2$dbs=="GO_BP_23",]
golist2 = golist2[golist2$Adjusted.P.value<0.2,]
# variable to split up and down regulated  
golist2$direction <- golist2$cellcntrst
golist2$direction[startsWith(golist2$direction, 'up', trim=TRUE)] <- 'up' 
golist2$direction[startsWith(golist2$direction, 'down', trim=TRUE)] <- 'down' 
# add celltype variable
golist2$cell <- golist2$cellcntrst
golist2$cell <- sub("Mixed", NA, golist2$cell)
golist2 = golist2 %>% drop_na(cell)
golist2$cell[grep("Astrocytes", golist2$cell)] <- "Astro"
golist2$cell[grep("Neuron_Dentate_C1ql2", golist2$cell)] <- "EN-DG"
golist2$cell[grep("Neuron_CA1_", golist2$cell)] <- "EN-CA1"
golist2$cell[grep("Entorhinal_Nxph3", golist2$cell)] <- "EN-Nxph3"
golist2$cell[grep("Microglia", golist2$cell)] <- "Micro"
golist2$cell[grep("Oligodendrocytes", golist2$cell)] <- "Oligo"
golist2$cell[grep("Neuron_CA2CA3", golist2$cell)] <- "EN-CA2CA3"
golist2$cell[grep("Polydendrocytes", golist2$cell)] <- "OPC"
golist2$cell[grep("Interneuron_Gad2-Sst", golist2$cell)] <- "IN-Sst"
golist2$cell[grep("Interneuron_Gad2-Vip", golist2$cell)] <- "IN-Vip"
golist2$cell[grep("Interneuron_Gad2-Lamp5", golist2$cell)] <- "IN-Lamp5"
golist2$cell[grep("Neuron_CajalRhetzius_Lhx1" , golist2$cell)] <- "EN-CajalRhetzius"
golist2$cell[grep("Endothelial-Mural-Fibroblast", golist2$cell)] <- "EndoMural"
golist2$cell[grep("Neuron_Subiculum_Slc17a6", golist2$cell)] <- "EN-L3EC"
#gsub(".*_scott80_.*", "incongruent", d)
# add contrast variable
golist2$contrast <- golist2$cellcntrst
golist2$contrast[grep("bothgroups", golist2$contrast)] <- "Age"
golist2$contrast[grep("sham.3movs18mo", golist2$contrast)] <- "Age-sham"
golist2$contrast[grep("TBI.3movs18mo", golist2$contrast)] <- "Age-TBI"
golist2$contrast[grep("bothages", golist2$contrast)] <- "TBI"
golist2$contrast[grep("TBIvssham.3mo", golist2$contrast)] <- "TBI-3mo"
golist2$contrast[grep("TBIvssham.18mo", golist2$contrast)] <- "TBI-18mo"
write.table(golist2,file=paste('tbi_pseudobulkDGEallcells_enriched_filtered_allPvals.txt',sep=''),sep='\t',quote = F)
grp_cell = paste( golist2$cell , golist2$contrast , sep='_' )
golist2$grp_cell = grp_cell
TermN = gsub(".*\\(", "", golist2$Term)
TermN = sub("\\).*", "", TermN)
golist2$TermN = TermN


length(unique(golist3$Term)) #488
#5691
grp_cell = paste( golist3$cell , golist3$contrast , sep='_' )
golist3$grp_cell = grp_cell
unique_terms <- unique(golist3$Term)
unique_grp_cell <- unique(golist3$grp_cell)
golist <- data.frame(Term = unique_terms)

# Add columns for each unique cellcntrst value
for (grp_cell in unique_grp_cell) {
  golist[, grp_cell] <- NA
  
  # Fill in Adjusted.P.values where matches are found
  for (i in 1:nrow(golist)) {
    term <- golist$Term[i]
    match_row <- golist3[golist3$Term == term & golist3$grp_cell == grp_cell, ]
    
    if (nrow(match_row) > 0) {
      golist[i, grp_cell] <- match_row$Adjusted.P.value[1]
    }
  }
}
golist2 = -log10(golist[,2:79])
golist2 = cbind(golist[,1], golist2)
golist2 = golist
golist2[is.na(golist2)] <- 1
# Calculate inter-rater reliability for each row: kapp-stat
library(irr)
sig = golist3[golist3$Adjusted.P.value < 0.05 , ]
u = unique(sig$grp_cell)
TermN = gsub(".*\\(", "", sig$Term)
TermN = sub("\\).*", "", TermN)
sig$TermN = TermN


TermN = gsub(".*\\(", "", golist2$Term)
TermN = sub("\\).*", "", TermN)
golist2$TermN = TermN
sets = strsplit(unique(sig$TermN) , split = ',' )

l = length(sets)
k = matrix( NA , nrow = l , ncol = l )
for( i in 1:l ) {
  cat( i , '\n' )
  for( j in 1:l ) {
    set1 = sig[sig$TermN==sets[[i]],]#sig[sig$grp_cell==sets[[i]],]#sig[sig$Term==sets[[i]],15]#sig[sig$cellcontrast==sets[[i]],1]#sets[[i]]
    set2 = sig[sig$TermN==sets[[j]],]#sig[sig$grp_cell==sets[[j]],]#sig[sig$Term==sets[[i]],15]#sig[sig$cellcontrast==sets[[j]],1]#sets[[j]]
    r1 = u %in% set1$grp_cell#set1$Term#sig[sig$Term==set1,] sig[sig$cellcontrast==set1,]
    r2 = u %in% set2$grp_cell#set2$Term#sig[sig$Term==set2,] sig[sig$cellcontrast==set1,]
    k2 = kappa2( ratings = cbind(r1,r2) )
    k[i,j] = k2$value
  }
}
diag(k) <- NA
k[upper.tri(k)] <- NA
k2 = k
k[k < 0] <- NA
col_sums <- colSums(k, na.rm = TRUE)
names(col_sums) <- sets
topgos = order(col_sums, decreasing = T)
topgos = topgos[1:100]
topgos = sets[topgos]

#order(colSums(k, na.rm = TRUE), decreasing = TRUE)
#1   86  199  663  664  674 1342  200  250  497   46
#sig[c(1,86,199,644,674,1342,200,250,497,46),1]
# Function to calculate ICC for a single row
calculate_icc <- function(row) {
  values <- as.numeric(row[-1])  # Exclude the 'Term' column
  values <- values[!is.na(values)]  # Remove NA values
  
  if (length(values) >= 2) {
    icc_result <- icc(matrix(values, nrow = 1), model="t", type="a", unit="a")
    return(icc_result$value)
  } else {
    return(NA)
  }
}
calculate_icc <- function(row) {
  values <- as.numeric(row[-1])  # Exclude the 'Term' column
  values <- values[!is.na(values)]  # Remove NA values
  
  if (length(values) >= 2) {
    icc_result <- icc(matrix(values, nrow = 1))
    return(icc_result$value)
  } else {
    return(NA)
  }
}
# Apply the function to each row and add the results as a new column
#ICC <- apply(golist2, 1, calculate_icc)
#kappam.fleiss(diagnoses[, 1:3])

golist4 = golist2[golist2$TermN %in% topgos,]


fill1 = c("#BBE7C8FF", "#60CEACFF","#40498EFF", "#F6B893FF", "#BD1655FF", "#F47F58FF")

golist4$logP = -log10(golist4$Adjusted.P.value)

golist= golist4
# Create color mapping based on contrast type and logP
golist$col <- NA
tbi_contrasts <- c("TBI", "TBI-3mo", "TBI-18mo")
age_contrasts <- c("Age", "Age-sham", "Age-TBI")

library(viridis)
golist$col[golist$contrast %in% tbi_contrasts] <- viridis::viridis(
  n = sum(golist$contrast %in% tbi_contrasts),
  option = "rocket",
  direction = -1
)[rank(golist$logP[golist$contrast %in% tbi_contrasts])]

golist$col[golist$contrast %in% age_contrasts] <- viridis::viridis(
  n = sum(golist$contrast %in% age_contrasts), 
  option = "mako",
  direction = -1
)[rank(golist$logP[golist$contrast %in% age_contrasts])]

# Create dircon column combining direction and contrast
golist$dircon <- paste(golist$direction, substr(golist$contrast, 1, 3), sep = "-")
golist$majcontrast <- substr(golist$contrast, 1, 3)
saveRDS(golist, "tbi_DGEallcells_golist_irrTop100.rds")
write.table(golist,file=paste('tbi_DGEallcells_golist_irrTop100.txt',sep=''),sep='\t',quote = F)

golist4 = golist[golist$majcontrast == "Age",]
golist5 = golist[golist$majcontrast == "TBI",]

golist6 <- read.delim("/Users/mprokhorenko/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/tbi/Monocle/tbi_allcells_golist_irrTop100_edited.txt")

# select top 12 GOs
#[1] "Axonogenesis (GO:0007409)"                                              
#[2] "Generation Of Neurons (GO:0048699)"                                     
#[3] "Negative Regulation Of Intracellular Signal Transduction (GO:1902532)"  
#[4] "Negative Regulation Of Macrophage Activation (GO:0043031)"              
#[5] "Negative Regulation Of Neuroinflammatory Response (GO:0150079)"         
#[6] "Positive Regulation Of Cytokine Production (GO:0001819)"                
#[7] "Positive Regulation Of Cytokine-Mediated Signaling Pathway (GO:0001961)"
#[8] "Positive Regulation Of Excitatory Postsynaptic Potential (GO:2000463)"  
#[9] "Positive Regulation Of Glial Cell Migration (GO:1903977)"               
#[10] "Positive Regulation Of Phagocytosis (GO:0050766)"                       
#[11] "Positive Regulation Of Synaptic Transmission (GO:0050806)"              
#[12] "Regulation Of Primary Metabolic Process (GO:0080090)"    
topgos = unique(golist6$Term)
golist6 = golist[golist$Term %in% topgos,]
golist4 = golist6[golist6$majcontrast == "Age",]
golist5 = golist6[golist6$majcontrast == "TBI",]


p <- ggplot(golist4) +
  aes(x = grp_cell, y = TermN, fill = col) +
  geom_tile() +
  scale_fill_identity() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
  facet_wrap(vars(direction))
png("allcellsDGE_irrGO_age6.png", width = 9, height = 4, units = 'in', res = 300)
p
dev.off()
p <- ggplot(golist5) +
  aes(x = grp_cell, y = TermN, fill = col) +
  geom_tile() +
  scale_fill_identity() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
  facet_wrap(vars(direction))
png("allcellsDGE_irrGO_tbi6.png", width = 9, height = 4, units = 'in', res = 300)
p
dev.off()

####
# print color scales
# viridis scale rocket for TBI
range(golist$logP[golist$contrast %in% tbi_contrasts])
z=matrix(1:100,nrow=1)
x=1
yy=seq(1.828820e-06, 1.165416e+01,len=100) 
image(x,yy,z,col=rocket(100),axes=FALSE,xlab="",ylab="")
axis(2)
#viridis scale mako for AGE
range(golist$logP[golist$contrast %in% age_contrasts])
#   1.713839e-06 5.553079e+00
z=matrix(1:100,nrow=1)
x=1
yy=seq(1.713839e-06, 5.553079e+00, len=100)
image(x,yy,z,col=mako(100),axes=FALSE,xlab="",ylab="")



range(golist$logP)
#1.713839e-06 1.165416e+01
summary(rank(golist$logP[golist$contrast %in% age_contrasts]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#38    1160    2317    2320    3473    4640 
summary(rank(golist$logP[golist$contrast %in% tbi_contrasts]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#21.5  1175.0  2350.0  2350.0  3525.0  4699.0 