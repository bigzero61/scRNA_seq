############Single-cell RNAseq analysis (GSE171169)################
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(viridis)
#pre-process using seurat
setwd("D:/scRNA/GSE171169_RAW")
counts.list <- list()
all_folders <- dir()
for(i in 1:length(all_folders)){
  counts.list <- Read10X(data.dir=all_folders[i])
  object <- CreateSeuratObject(counts=counts.list,project = all_folders[i])
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  print(VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  print(FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  assign(all_folders[i],object)
}
rm(list=c("counts.list","object","i","all_folders"))
mcao_5d <- merge(`5d_N1`,`5d_N2`,add.cell.ids=c("n1","n2"),project="mcao_5d")
rm(list=c("5d_N1","5d_N2"))
mcao_5d <- subset(mcao_5d, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
mcao_5d <- NormalizeData(mcao_5d, normalization.method = "LogNormalize", scale.factor = 10000)
mcao_5d <- FindVariableFeatures(mcao_5d, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mcao_5d), 10)
LabelPoints(plot = VariableFeaturePlot(mcao_5d), points = top10, repel = TRUE)
all.genes <- rownames(mcao_5d)
mcao_5d <- ScaleData(mcao_5d, features = all.genes)
mcao_5d <- RunPCA(mcao_5d)
ElbowPlot(mcao_5d, ndims = ncol(Embeddings(mcao_5d, "pca")))
mcao_5d <- RunUMAP(mcao_5d, reduction = "pca", dims = 1:30)
mcao_5d <- RunTSNE(mcao_5d,dims = 1:30)
mcao_5d <- FindNeighbors(mcao_5d, reduction = "pca", dims = 1:30)
mcao_5d <- FindClusters(mcao_5d, resolution = c(0.4,0.6,0.8,1.0))
Idents(mcao_5d) <- mcao_5d$RNA_snn_res.0.4
mcao_5d.markers <- FindAllMarkers(mcao_5d, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
top_markers <- mcao_5d.markers %>%
  group_by(cluster) %>%
slice_max(n = 5, order_by = avg_log2FC)
mcao_5d <- RenameIdents(mcao_5d,"0"="MG 1","1"="MG 2","5"="MG 3","2"="MΦ 1","3"="MΦ 2","4"="MΦ 3","9"="MΦ 4","6"="NK cell",
                        "7"="Proliferating monocytes","8"="T cell","10"="B cell","12"="Neutrophil","11"="MΦ 1","13"="Unknown","14"="Unknown","15"="Unknown")
mcao_5d@meta.data$cell_type <- mcao_5d@active.ident
mcao_5d <- subset(mcao_5d,invert=TRUE,subset = cell_type=="Unknown")
color_tsne <- viridis(13, alpha = 0.7, begin = 0, end = 1, direction = 1, option = "turbo")
DoHeatmap(mcao_5d, features = top_markers$gene,group.colors = color_tsne,label = FALSE)+
  scale_fill_viridis(alpha = 1,direction = 1,end = 1,option = "viridis",na.value="white")
TSNEPlot(mcao_5d,label=FALSE,cols=color_tsne,pt.size=0.7)+NoLegend()
p <- FeaturePlot(mcao_5d,features = c("P2ry12","Tmem119","Gpr34","Cx3cr1"),label=FALSE,reduction = "tsne",ncol = 1,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+scale_colour_viridis(alpha=1,direction=-1,begin=0.5,option="mako")
}
cowplot::plot_grid(plotlist = p,ncol = 2)

#######microglia subclusters
microglia <- subset(mcao_5d,subset = cell_type==c("MG 1","MG 2","MG 3"))
homeo <- c("Cx3cr1","P2ry12","P2ry13","Gpr34","Tmem119","Selplg","Olfml3")
p1 <- DotPlot(microglia,features = homeo,dot.scale = 8)+
  scale_colour_viridis(alpha = 0.8,direction=-1,begin=0.5,option = "mako")+theme(legend.position = "top",axis.title.y=element_blank(),
                                                                                 axis.text.y=element_blank(),
                                                                                 axis.ticks.y=element_blank(),
                                                                                 axis.line.y = element_blank(),
                                                                                 axis.title.x = element_blank(),
                                                                                 axis.text.x = element_text(angle=45,size=20,vjust = 0.6),
                                                                                 #legend.text = element_text(size = 12),
                                                                                 #legend.title = element_text(size = 16),
                                                                                 #axis.text.y=element_text(size = 20),
                                                                                 legend.key.width= unit(0.7, 'cm'),
                                                                                 legend.title = element_blank())

pro <- c("Itgax","Lyz2","Il1b","S100a8","Cd74","Cd86","Ccr2","Nfkb1","Mrc1","Cd163","Tgfbi")
p2 <- DotPlot(microglia,features = pro,dot.scale = 8)+
  scale_colour_viridis(alpha = 0.8,direction=-1,begin=0.5,option = "rocket")+theme(legend.position = "top",axis.title.y=element_blank(),
                                                                                   axis.text.y=element_blank(),
                                                                                   axis.ticks.y=element_blank(),
                                                                                   axis.line.y = element_blank(),
                                                                                   axis.title.x = element_blank(),
                                                                                   axis.text.x = element_text(angle=45,size=20,vjust = 0.6),
                                                                                   legend.key.width= unit(0.7, 'cm'),
                                                                                   legend.justification = "right",
                                                                                   #legend.text = element_text(size = 12),
                                                                                   legend.title = element_blank())
cowplot::plot_grid(p1,p2,ncol=2,rel_widths = c(7,11))

microglia_all_markers <- FindAllMarkers(microglia,only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
microglia_top_markers <- microglia_all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
DoHeatmap(microglia, features = microglia_top_markers$gene,group.colors = color_tsne[1:3],label = TRUE,angle=0,hjust = 1)+
  scale_fill_viridis(alpha = 1,direction = 1,end = 1,option = "viridis",na.value="white")+NoLegend()

####### MG 3 enrichment analysis
####### KEGG over-representatation analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
microglia.3 <- FindMarkers(microglia,ident.1 = "MG 3") %>% filter(p_val_adj<0.05)
write.csv(microglia.3,"MG3deg.csv")
genedf <- bitr(rownames(microglia.3),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
microglia.3$symbol <- rownames(microglia.3)
microglia.3 <- merge(x=microglia.3,y=genedf,by.x="symbol",by.y="SYMBOL")
genelist <- microglia.3$avg_log2FC
names(genelist) <- microglia.3$ENTREZID
genelist = sort(genelist,decreasing = T)
kk <- enrichKEGG(gene         = microglia.3$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

####### GSVA analysis
library(GSVA)
library(msigdbr)
library(pheatmap)
m <- as.matrix(microglia@assays$RNA@counts)
cluster <- microglia@meta.data[,"cell_type"]
mouse_KEGG <- msigdbr(species = "Mus musculus",
                      category = "C2", 
                      subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)
mouse_KEGG_set <- mouse_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
gsva_KEGG <- gsva(expr = m, 
                  gset.idx.list = mouse_KEGG_set,
                  kcdf="Poisson", #查看帮助函数选择合适的kcdf方法 
                  parallel.sz = 5)
kegg <- as.data.frame(gsva_KEGG) 
m_data <- kegg[c("KEGG_MAPK_SIGNALING_PATHWAY","KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                 "KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                 "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_P53_SIGNALING_PATHWAY"
                 #,"KEGG_PPAR_SIGNALING_PATHWAY","KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_MAPK_SIGNALING_PATHWAY","KEGG_MTOR_SIGNALING_PATHWAY","KEGG_NOTCH_SIGNALING_PATHWAY"
                 ),]
annotation_col <- data.frame(cluster)
rownames(annotation_col) <- colnames(gsva_KEGG)
m_data <- t(m_data)
m_data <- as.data.frame(m_data)
m_data <- cbind(m_data,annotation_col)
m_data <- arrange(m_data,cluster)
m_data <- m_data[,-grep("cluster",colnames(m_data))]
m_data <- as.matrix(t(m_data))
annotation_color <- list(cluster=c(`MG 1`="#30123BB3",`MG 2`="#434FBBB3",`MG 3`="#4686FBB3"))
kk2 <- pheatmap(kegg,
                cluster_rows = F,
                cluster_cols = F,
                annotation_col =annotation_col,
                annotation_legend=F, 
                show_rownames = T,
                show_colnames = F,
                annotation_colors = annotation_color,
                color =colorRampPalette(rev(brewer.pal(n = 7, name =
                                                         "RdBu")))(100),
                cellwidth = 0.5, cellheight = 13,
                fontsize = 10,
                legend = T)
#annotation_legend=FALSE)
jakstat_gene <- strsplit(x=kk_selected["JAK-STAT signaling pathway","geneID"],"/")[[1]]
m_data <- as.data.frame(microglia@assays$RNA@data)
jakstat <- bitr(jakstat_gene,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
m_data <- m_data[jakstat$SYMBOL,]
m_data <- as.matrix(m_data)
#column ordered by cluster
m_data <- t(m_data)
m_data <- as.data.frame(m_data)
m_data <- cbind(m_data,annotation_col)
m_data <- arrange(m_data,cluster)
m_data <- m_data[,-grep("cluster",colnames(m_data))]
m_data <- as.matrix(t(m_data))
#row ordered by fold change
m_data <- as.data.frame(m_data)
rownames(microglia.3) <- microglia.3$symbol
microglia.3_selected <- microglia.3[rownames(m_data),]
microglia.3_selected <- arrange(microglia.3_selected,-avg_log2FC)
mmm <- as.data.frame(microglia.3_selected[,3])
rownames(mmm) <- rownames(microglia.3_selected)
m_data$gene <- rownames(m_data)
mmm$gene <- rownames(mmm)

m_data <- merge(x=m_data,y=mmm,by.x="gene",by.y="gene")
rownames(m_data) <- m_data$gene
m_data <- m_data[,-1]

m_data <- arrange(m_data,-`microglia.3_selected[, 3]`)
m_data <- m_data[,-743]
m_data <- as.matrix(m_data)
geneHM <- pheatmap(m_data,
                   cluster_rows = F,
                   cluster_cols = F,
                   annotation_col =annotation_col,
                   annotation_legend=F, 
                   show_rownames = T,
                   show_colnames = F,
                   annotation_colors = annotation_color,
                   color =colorRampPalette(rev(brewer.pal(n = 7, name =
                                                            "RdBu")))(100),
                   cellwidth = 0.5, cellheight = 13,
                   fontsize = 10,
                   legend = TRUE)
#### bar plot showing kegg over-representaton analysis of selected terms
kk_selected <- as.data.frame(kk) %>% filter(Description=c())
ggplot(data = kk_selected,aes(x=Description,y=Count))+
  geom_bar(stat = "identity",aes(fill=-log10(kk_selected$p.adjust)))+
  coord_flip()+theme_void()+
  geom_text(mapping=aes(y=Count+1.5),label=kk_selected$Count,size=8)+
  scale_fill_viridis(option = "rocket",direction = -1,alpha = 1,begin = 0.5)+
  ggtitle("Gene Count")+theme(plot.title=element_text(hjust = 0.5))

######## Il13ra1
m <- as.matrix(microglia@assays$RNA@data)
microglia$target <- m["Il13ra1",]
test <- microglia
test$target[which(test$target>0)] <- "pos"
test$target[which(test$target==0)] <- "non"
Idents(test) <- "target"
pos.non <- test %>% 
  FindMarkers(ident.1 = "pos",ident.2 = "non",logfc.threshold=0) %>% 
  filter(p_val_adj<0.05)
write.csv(pos.non,"il13ra1_degs.csv")
genedf <- bitr(rownames(pos.non),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
pos.non$symbol <- rownames(pos.non)
pos.non <- merge(x=pos.non,y=genedf,by.x="symbol",by.y="SYMBOL")
genelist <- pos.non$avg_log2FC
names(genelist) <- pos.non$ENTREZID
genelist = sort(genelist,decreasing = T)
gsekk_target <- gseKEGG(geneList     = genelist,
                        organism     = 'mmu',
                        minGSSize    = 5,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
go2 <- gseGO(geneList     = genelist,
                     OrgDb        = org.Mm.eg.db,
                     ont          = "BP",
                     minGSSize    = 100,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)
table <- gsekk_target@result
gseaplot2(gsekk_target,geneSetID = c("mmu04621","mmu04620","mmu04115","mmu04062","mmu04010"),pvalue_table = TRUE)
dotplot(gsekk_target, showCategory=30) + ggtitle("Dotplot for GSEA (Il13ra1+ vs Il13ra1-)")
ridgeplot(gsekk_target)


#######fgsea
mouse_KEGG <- msigdbr(species = "Mus musculus",
                      category = "C2", 
                      subcategory = "CP:KEGG")
mdb_kegg = mouse_KEGG [grep("^KEGG",mouse_KEGG$gs_name),]
fgsea_sets<- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)
pos.non$genes = rownames(pos.non)
cluster0.genes<- pos.non %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
ranks<- deframe(cluster0.genes)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
