#ANALISIS TRANSKRIPTOMIK HIV
#Dataset: GSE108296 (HIV Controller vs cART-treated patients)
#Platform: Microarray (Illumina GPL10558)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) dan Enrichment 

#PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE)
#1. Install BiocManager (manajer paket Bioconductor)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
#2. Install paket Bioconductor (GEOquery & limma)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update =
                       FALSE)
#Install annotation package sesuai platform
#GPL10558= Illumina HumanHT-12 V4.0 expression beadchip
BiocManager::install("illuminaHumanv4.db", ask = FALSE, update = FALSE)
#3. Install paket CRAN untuk visualisasi dan manipulasi data
install.packages(c("pheatmap", "ggplot2", "dplyr"))
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
#4. Memanggil library
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(umap)
library(illuminaHumanv4.db)

#PENGAMBILAN DATA DARI GEO
gset <- getGEO("GSE108296", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
#PRE-PROCESSING DATA EKSPRESI
ex <- exprs(gset)
qx <- quantile(ex,c(0,0.25,0.5,0.75,0.99,1),na.rm=TRUE)
#LogTransform 
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#METADATA & GROUPING
metadata <- pData(gset)
group_raw <- metadata$characteristics_ch1.1
group <- ifelse(
  grepl("Controller", group_raw, ignore.case = TRUE),
  "Controller",
  "cART"
)
gset$group <- factor(group)
table(gset$group)

#DESIGN MATRIX (KERANGKA STATISTIK)
#model.matrix():membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + group, data=gset)
#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)
contrast_matrix <- makeContrasts(
  Controller_vs_cART = Controller - cART,
  levels = design
)

#ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
#lmFit(): membangun model linear untuk setiap gen
fit <- lmFit(ex, design)
#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)
#eBayes():empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)
#topTable():mengambil hasil akhir DEG
deg_results <- topTable(fit2,
                        adjust="fdr",
                        number=Inf)
deg_filtered <- deg_results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

#ANOTASI NAMA GEN
#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(deg_filtered)
#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL","GENENAME"),
  keytype = "PROBEID"
)
#Gabungkan dengan hasil limma
deg_filtered$PROBEID <- rownames(deg_filtered)
deg_annotated <- merge(
  deg_filtered,
  gene_annotation,
  by="PROBEID",
  all.x=TRUE
)

#BOXPLOT DISTRIBUSI NILAI EKSPRESI
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  cex.axis = 0.6,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

#UMAP (VISUALISASI DIMENSI RENDAH)
umap_input <- t(ex)
#Jalankan UMAP
umap_result <- umap(umap_input)
#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#VISUALISASI VOLCANO PLOT
volcano_data <- deg_results
#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[
  volcano_data$logFC > 1 &
    volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[
  volcano_data$logFC < -1 &
    volcano_data$adj.P.Val < 0.05] <- "DOWN"
#Visualiasi
ggplot(volcano_data,
       aes(logFC,
           -log10(adj.P.Val),
           color=status))+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c(
    "DOWN"="blue",
    "NO"="grey",
    "UP"="red"
  ))+
  geom_vline(xintercept=c(-1,1),
             linetype="dashed")+
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed")+
  theme_minimal()+
  ggtitle("Volcano Plot Controller vs cART")

#VISUALISASI HEATMAP
#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
top50 <- deg_filtered %>%
  arrange(adj.P.Val) %>%
  head(50)
#Ambil ID probe dari top 50
probe_ids <- rownames(top50)
#Terjemahkan ID Probe Illumina ke Simbol Gen resmi
gene_symbols <- mapIds(illuminaHumanv4.db, 
                       keys = probe_ids, 
                       column = "SYMBOL", 
                       keytype = "PROBEID", 
                       multiVals = "first")
#Jika ada probe yang tidak memiliki simbol gen, tetap gunakan ID aslinya 
gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]
#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[probe_ids, ]
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]
rownames(mat_heatmap) <- gene_symbols
#Visualisasi heatmap
pheatmap(mat_heatmap,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Top 50 DEG")

#ANALISIS GENE ONTOLOGY (GO) & KEGG
# Install package enrichment
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), 
                     ask = FALSE, update = FALSE)
# Load library
library(clusterProfiler)
library(org.Hs.eg.db)
# Ambil gene symbol unik dan tidak NA
deg_symbols <- unique(deg_annotated$SYMBOL)
deg_genes <- deg_symbols[!is.na(deg_symbols)]
length(deg_genes)  
#Konversi symbol menjadi entrez id
gene_df <- bitr(
  deg_symbols,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)
entrez_ids <- gene_df$ENTREZID
#Analisis GO
ego <- enrichGO(
  gene=entrez_ids,
  OrgDb=org.Hs.eg.db,
  ont="ALL",
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  qvalueCutoff=0.05,
  readable=TRUE
)
#Visualisasi barplot dan dotplot
barplot(ego, showCategory = 15, label_format = 30) + 
  theme(axis.text.y = element_text(size = 8))
dotplot(ego, showCategory = 15, label_format = 30) + 
  theme(axis.text.y = element_text(size = 8))

#Analisis KEGG
ekegg <- enrichKEGG(
  gene=entrez_ids,
  organism="hsa",
  pvalueCutoff=0.05
)
# Konversi agar gene symbol terbaca
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(ekegg)
#Visualisasi dotplot & barplot
barplot(ekegg,showCategory = 15)
dotplot(ekegg,showCategory = 15)

#MENYIMPAN HASIL
# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(deg_annotated,"DEG_Controller_vs_cART.csv")
write.csv(as.data.frame(ego),"GO_Controller_vs_cART.csv")
write.csv(as.data.frame(ekegg),"KEGG_Controller_vs_cART.csv")
message("Analisis selesai. File hasil telah disimpan.")