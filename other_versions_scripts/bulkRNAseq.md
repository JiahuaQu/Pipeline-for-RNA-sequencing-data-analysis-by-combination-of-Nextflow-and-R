Ref: https://github.com/JiahuaQu/Pipeline-for-RNA-sequencing-data-analysis-by-combination-of-Nextflow-and-R/blob/master/RNA-seq%20pipeline.md

# Upstream

human UCSF h38 genome

nf-core-rnaseq_3.10.1

```bash
ln -s /research/sharedresources/immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/ /home/jqu/project/OLAH/analysis/

cd 	/research/sharedresources/immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/

mkdir 01-Nextflow
cd 01-Nextflow
module load mamba/1.5.8
mamba create -n NextFlow-nf-core-Bulk_RNAseq
mamba activate NextFlow-nf-core-Bulk_RNAseq
python -V     # Python 3.10.14

# install the latest version of nextflow,
# rather than loading the old version of module.
mamba install nextflow -y
mamba list | grep nextflow   
# nextflow                  24.10.5              hdfd78af_0    bioconda

# install nf-core
mamba search nf-core
mamba install nf-core -y
mamba list | grep nf-core
# nf-core                   3.2.0              pyhdfd78af_0    bioconda

nf-core download nf-core/rnaseq
# select version：3.10.1 
# OR: 3.18.0
# Singularity
# Copy
# Not zip
```

```markdown
INFO     Saving 'nf-core/rnaseq'                                                      
          Pipeline revision: '3.18.0'                                                 
          Use containers: 'singularity'                                               
          Container library: 'quay.io'                                                
          Using $NXF_SINGULARITY_CACHEDIR': /home/jqu/NXF_SINGULARITY_CACHEDIR'       
          Output directory: 'nf-core-rnaseq_3.18.0'                                   
          Include default institutional configuration: 'False'                        
INFO     Downloading workflow files from GitHub                                       
INFO     Processing workflow revision 3.18.0, found 42 container images in total.
```

```markdown
INFO     Saving 'nf-core/rnaseq'                                                      
          Pipeline revision: '3.10.1'                                                 
          Use containers: 'singularity'                                               
          Container library: 'quay.io'                                                
          Using $NXF_SINGULARITY_CACHEDIR': /home/jqu/NXF_SINGULARITY_CACHEDIR'       
          Output directory: 'nf-core-rnaseq_3.10.1'                                   
          Include default institutional configuration: 'False'                        
INFO     Downloading workflow files from GitHub                                       
INFO     Processing workflow revision 3.10.1, found 31 container images in total.
```



```bash
# reference genome
# should download human genome
cd /research/sharedresources/immunoinformatics/common/jqu/reference_genome/UCSC/AWS_iGenomes/human/hg38

module load aws-cli/2.22.12

nohup aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/ ./references/Homo_sapiens/UCSC/hg38/Annotation/Genes/ --exclude "*" --include "genes.gtf" &
# /research/sharedresources/immunoinformatics/common/jqu/reference_genome/UCSC/AWS_iGenomes/human/hg38/references/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf

nohup aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ ./references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ &
# /research/sharedresources/immunoinformatics/common/jqu/reference_genome/UCSC/AWS_iGenomes/human/hg38/references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa


# fastq
# /research/groups/crawfgrp/projects/crawfgrp_hartwell/common/illumina/crawfgrp_842475_RNAseq_total_stranded_illumina-1/
mkdir fastq/
ln -s /research/groups/crawfgrp/projects/crawfgrp_hartwell/common/illumina/crawfgrp_842475_RNAseq_total_stranded_illumina-1/*/*.fastq.gz fastq/

# sample sheet
/research/sharedresources/immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/01-Nextflow

wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py

chmod a+x fastq_dir_to_samplesheet.py

./fastq_dir_to_samplesheet.py \
-r1 _R1_001.fastq.gz \
-r2 _R2_001.fastq.gz \
fastq \
samplesheet.csv

# Process Multiple runs of the same sample
awk 'BEGIN {FS=OFS=","} NR==1 {print; next} { $1 = substr($1, 1, length($1)-8); print }' samplesheet.csv > new_samplesheet.csv
# Rename sample name in excel

mkdir output

vim nextflow.config
```

```bash
// Global default params, used in configs
params {
 
    // Input options
    input = 'new_samplesheet.csv'
 
    // References
    fasta = '/research/sharedresources/immunoinformatics/common/jqu/reference_genome/UCSC/AWS_iGenomes/human/hg38/references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
    gtf = '/research/sharedresources/immunoinformatics/common/jqu/reference_genome/UCSC/AWS_iGenomes/human/hg38/references/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
	gencode = false
	save_reference = true
	     
	// QC
	skip_bigwig = true
	skip_stringtie = true
	skip_preseq = true
	skip_qualimap = true
	deseq2_vst = true
	
	// Boilerplate options
	outdir = 'output'
	email = "jiahua.qu@stjude.org"
	 
	// Max resource options
	// Defaults only, expecting to be overwritten
	max_memory = '500.GB'
	max_cpus = 32
	max_time = '240.h'
}
```

```bash
cd /research/sharedresources/immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/01-Nextflow
vim my_script.sh
```

```bash
#!/bin/bash
#BSUB -P OLAH
#BSUB -J nextflow
#BSUB -q superdome
#BSUB -n 32
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error

module load conda3/202402
source activate NextFlow-nf-core-Bulk_RNAseq

cd /research/sharedresources/immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/01-Nextflow

# set path
log="./my.log"
workflow="./nf-core-rnaseq_3.10.1/3_10_1"
profile=singularity
config="./nextflow.config"

# run pipeline
nextflow -log $log run $workflow -profile $profile -resume -c $config

conda deactivate
```

# Midstream and Downstream

```R
# path_script <- "/home/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq-script/R/"
# path_save <- "/mnt/iminfo/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/"
path_save <- "Z:/ResearchHome/SharedResources/Immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/"

### QC
setwd(path_save)
if (!dir.exists("02-QC")){
  dir.create("02-QC")
} 
setwd("./02-QC")
library(data.table)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
gene_counts <- readRDS("../01-Nextflow/output/star_salmon/salmon.merged.gene_counts.rds")  

# extract and process count matrix
counts <- gene_counts@assays@data$counts

temp <- counts[rownames(counts) == "OLAH",]
fwrite(temp, "raw_counts-OLAH.csv", row.names = T)

# round non-integer to integer for loading into DESeq2
counts_round <- round(counts, digits = 0)

# remove rows whose values are all zero
library(dplyr)
counts_round <- counts_round %>%
  filter(rowSums(.) > 0)
temp <- counts_round[rownames(counts_round) == "OLAH",]
fwrite(temp, "counts_round-OLAH.csv", row.names = T)

# meta_data
meta_data <- data.frame(
  sample=colnames(counts_round),
  condition=factor(c("Ctrl","Ctrl","Ctrl","Ctrl","OLAH","OLAH","OLAH","OLAH"), levels = c("Ctrl", "OLAH"))
)
saveRDS(meta_data, "meta_data.rds", compress = T)

#BiocManager::install("DESeq2")
library(DESeq2)
# Load in Salmon output
dds_new <- DESeqDataSetFromMatrix(countData = counts_round,
                                  colData = meta_data,
                                  design= ~ condition)

# Normalization
dds_new <- estimateSizeFactors(dds_new)
dds_rlog <- rlog(dds_new)
dds_vst <- varianceStabilizingTransformation(dds_new)
assay(dds_new, "rlog") <- assay(dds_rlog)
assay(dds_new, "vst") <- assay(dds_vst)
saveRDS(dds_new, "dds_new.rds", compress = T)
saveRDS(dds_rlog, "dds_rlog.rds", compress = T)
saveRDS(dds_vst, "dds_vst.rds", compress = T)

### create a folder
if (!dir.exists("QCplots")){
  dir.create("QCplots")
}

### 1. boxplot
#library(DESeq2)
## 1.1. counts_vst
counts_vst <- dds_new@assays@data$vst

library(tidyverse)
library(tibble)
dat <- counts_vst %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  gather(key = "sample", value = "counts_vst", -gene) %>%
  left_join(meta_data, by="sample")

#library(ggplot2)
p_1 <- ggplot(data=dat,aes(x=sample,y=counts_vst,fill=condition)) +
  geom_boxplot()+
  labs(x="Sample",y="VST Normalized Read Counts")+
  ggtitle("VST Normalized Read Count by Sample") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "QCplots/boxplot-counts_vst.pdf", plot = p_1, width = 10, height = 5)
ggsave(filename = "QCplots/boxplot-counts_vst.png", plot = p_1, width = 10, height = 5, units = "in", dpi = 300)

## 1.2. counts_rlog
counts_rlog <- dds_new@assays@data$rlog
#library(tibble)
dat <- counts_rlog %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  gather(key = "sample", value = "counts_rlog", -gene) %>%
  left_join(meta_data, by="sample")

#library(ggplot2)
p_2 <- ggplot(data=dat,aes(x=sample,y=counts_rlog,fill=condition)) +
  geom_boxplot()+
  labs(x="Sample",y="rlog Normalized Read Counts")+
  ggtitle("rlog Normalized Read Count by Sample") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "QCplots/boxplot-counts_rlog.pdf", plot = p_2, width = 10, height = 5)
ggsave(filename = "QCplots/boxplot-counts_rlog.png", plot = p_2, width = 10, height = 5, units = "in", dpi = 300)


### 2. sample correlation heatmap
sampleDists <- dist(t(counts_vst))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
mydata_col = data.frame(condition=meta_data$condition, row.names = meta_data$sample)

library(pheatmap)
p_3 <- pheatmap(
  sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors,
  annotation_col = mydata_col,
  main="Euclidean distance between samples"
)

### make my customized function to save heatmap
save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_3, "QCplots/sample_heatmap.pdf", 5, 4.1)

save_pheatmap_png <- function(x, filename, width, height, dpi = 300) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  
  # 计算像素尺寸
  width_px <- width * dpi
  height_px <- height * dpi
  
  # 使用 png 设备保存图像
  png(filename, width = width_px, height = height_px, res = dpi, units = "px")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# 调用示例
save_pheatmap_png(p_3, "QCplots/sample_heatmap.png", 5, 4.1)



### 3. gene expression heatmap of top 1000 variable genes
counts_vst_var <- as.data.frame(counts_vst)
counts_vst_var$var <- apply(counts_vst,1,var)
counts_vst_var <- counts_vst_var %>%
  arrange(desc(var))
saveRDS(counts_vst_var, "counts_vst_var.rds", compress=T)
counts_vst_var <- counts_vst_var %>%
  slice(1:1000) %>%
  as.matrix()

p_4 <- pheatmap(counts_vst_var[,1:8],
                scale="row",
                color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
                annotation_col = mydata_col,
                main="Gene Expression (VST) of Top 1000 Variable Genes",
                show_rownames = F)
save_pheatmap_pdf(p_4, "QCplots/gene_heatmap.pdf", 5, 6)
save_pheatmap_png(p_4, "QCplots/gene_heatmap.png", 5, 6)

### 4. PCA plots
# make my customized function to plot PCA
plotPCA_vst <- function (object, assay, ntop = 1000) {
  rv <- rowVars(assay(object, assay))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  df <- cbind( as.data.frame(colData(object)), pca$x)
  # order points so extreme samples are more likely to get label
  ord <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
  df <- df[ord,]
  attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
  return(df)
}

pca.data <- plotPCA_vst(dds_new, assay="vst", ntop=1000)
saveRDS(pca.data, "pca.data.rds", compress = T)

# getwd()
# setwd("/mnt/iminfo/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/02-QC")
# pca.data <- readRDS("pca.data.rds")

for (i in c(1:4)) {
  for (j in c((i+1):5)) {
    percentVar <- round(attr(pca.data, "percentVar")$percentVar)
    p <- ggplot(pca.data, aes_string(x=paste0("PC",i), y=paste0("PC",j), color="condition")) +
      geom_point(size=3) +
      geom_text(aes(label = sample), vjust=-0.5, hjust=0.5) +
      xlab(paste0("PC",i,": ",percentVar[i],"% variance")) +
      ylab(paste0("PC",j,": ",percentVar[j],"% variance")) +
      theme_bw()
    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".pdf"), plot = p, width=5, height=4)
    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".png"), plot = p, width=5, height=4, units = "in", dpi = 300)
  }
}

### create folders
setwd("../")
if (!dir.exists("03-DE")){
  dir.create("03-DE")
} 
setwd("./03-DE")

### read in gene information obtained from Nextflow
bg_ID <- read_delim("../01-Nextflow/output/star_salmon/salmon_tx2gene.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE)
names(bg_ID) <- c("transcript", "gene", "gene_name")
identical(bg_ID$gene, bg_ID$gene_name)   # TRUE
bg_ID <- bg_ID[,-1] %>% unique()

# library(clusterProfiler, lib.loc = "/home/jqu/R/workbench/4.3.2")
# library(org.Mm.eg.db, lib.loc = "/home/jqu/R/workbench/4.3.2")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#keytypes(org.Mm.eg.db)
bg_ENTREZID <- bitr(bg_ID$gene,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
bg_ID <- left_join(bg_ID, bg_ENTREZID, by=c("gene"="SYMBOL"))
saveRDS(bg_ID, "bg_ID.rds", compress = T)

### create folders
if (!dir.exists("DEfiles")){
  dir.create("DEfiles")
} 

if (!dir.exists("MAplots")){
  dir.create("MAplots")
} 

if (!dir.exists("Volcanoplots")){
  dir.create("Volcanoplots")
} 

if (!dir.exists("Heatmaps")){
  dir.create("Heatmaps")
} 

### generate FC and FDR via DESeq2
library(DESeq2)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
# read in dds object and metadata
dds_new <- readRDS("../02-QC/dds_new.rds")
meta_data <- readRDS("../02-QC/meta_data.rds")

### run DESeq to get result objects in loop
library(ggplot2)
library(readr)
library(tibble)
library(pheatmap)
counts_vst <- dds_new@assays@data$vst
dds_new_res <- DESeq(dds_new)
resultsNames(dds_new_res)
# "Intercept"              "condition_OLAH_vs_Ctrl"
saveRDS(dds_new_res, "dds_new_res.rds", compress = T)

### 1. DE results
res <- results(dds_new_res)
saveRDS(res, paste0("DEfiles/DE-",sampleB,"vs",sampleA,".rds"), compress = T)

# join gene_name
res_df <- as.data.frame(res)
res_symbol <- res_df %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  left_join(bg_ID, by="gene")

### 2. MA plots
pdf(file = paste0("MAplots/MAplot-",sampleB,"vs",sampleA,".pdf"),
    width=7, height=5)
MAplot <- plotMA(object=res, alpha=0.05,
                 main=paste0(sampleB,"vs",sampleA),
                 colNonSig="black",
                 colSig="red",
                 colLine="blue")
dev.off()

### 3. volcano plots
vol <- EnhancedVolcano(res_symbol,
                       lab = NA,
                       x = 'log2FoldChange',
                       y = 'padj',
                       ylab = bquote(~-Log[10] ~ italic(FDR)),
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = paste0(sampleB,"vs",sampleA),
                       legendLabels=c('Not DE &\nabsolute log2FC < 1',
                                      'Not DE &\nabsolute log2FC > 1',
                                      'FDR < 0.05 &\nabsolute log2FC < 1',
                                      'FDR < 0.05 &\nabsolute log2FC > 1'),
                       legendPosition = 'right')
ggsave(filename = paste0("Volcanoplots/Volcanoplot-",sampleB,"vs",sampleA,".pdf"), plot=vol, width = 8, height = 7)
ggsave(filename = paste0("Volcanoplots/Volcanoplot-",sampleB,"vs",sampleA,".png"), plot=vol, width = 8, height = 7, units = "in", dpi=300)

vol_lab <- EnhancedVolcano(res_symbol,
                           lab = res_symbol$gene_name,
                           x = 'log2FoldChange',
                           y = 'padj',
                           ylab = bquote(~-Log[10] ~ italic(FDR)),
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           title = paste0(sampleB,"vs",sampleA),
                           legendLabels=c('Not DE &\nabsolute log2FC < 1',
                                          'Not DE &\nabsolute log2FC > 1',
                                          'FDR < 0.05 &\nabsolute log2FC < 1',
                                          'FDR < 0.05 &\nabsolute log2FC > 1'),
                           legendPosition = 'right')
ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-",sampleB,"vs",sampleA,".pdf"), plot=vol_lab, width = 8, height = 7)
ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-",sampleB,"vs",sampleA,".png"), plot=vol_lab, width = 8, height = 7, units = "in", dpi=300)









### Excluse two outlier samples
# 
colnames(counts)
# "3259029_JCC489_Control-GFP_DOX-4"
# "3259033_JCC489_OLAH-GFP_DOX-4" 
### QC
setwd(path_save)

if (!dir.exists("04-exclude_outliers")){
  dir.create("04-exclude_outliers")
} 

setwd("./04-exclude_outliers")

if (!dir.exists("02-QC")){
  dir.create("02-QC")
} 

setwd("./02-QC")
library(data.table)

gene_counts <- readRDS("../01-Nextflow/output/star_salmon/salmon.merged.gene_counts.rds")  
counts <- gene_counts@assays@data$counts
counts <- counts[,-c(4,8)]
saveRDS(counts, "counts-exclude_outliers.rds", compress = T)

# round non-integer to integer for loading into DESeq2
counts_round <- round(counts, digits = 0)

# remove rows whose values are all zero
library(dplyr)
counts_round <- counts_round %>%
  filter(rowSums(.) > 0)
temp <- counts_round[rownames(counts_round) == "OLAH",]
fwrite(temp, "counts_round-OLAH.csv", row.names = T)

# meta_data
meta_data <- data.frame(
  sample=colnames(counts_round),
  condition=factor(c("Ctrl","Ctrl","Ctrl","OLAH","OLAH","OLAH"), levels = c("Ctrl", "OLAH"))
)
saveRDS(meta_data, "meta_data.rds", compress = T)

#BiocManager::install("DESeq2")
library(DESeq2)
# Load in Salmon output
dds_new <- DESeqDataSetFromMatrix(countData = counts_round,
                                  colData = meta_data,
                                  design= ~ condition)

# Normalization
dds_new <- estimateSizeFactors(dds_new)
dds_rlog <- rlog(dds_new)
dds_vst <- varianceStabilizingTransformation(dds_new)
assay(dds_new, "rlog") <- assay(dds_rlog)
assay(dds_new, "vst") <- assay(dds_vst)
saveRDS(dds_new, "dds_new.rds", compress = T)
saveRDS(dds_rlog, "dds_rlog.rds", compress = T)
saveRDS(dds_vst, "dds_vst.rds", compress = T)

### create a folder
if (!dir.exists("QCplots")){
  dir.create("QCplots")
}

### 1. boxplot
#library(DESeq2)
## 1.1. counts_vst
counts_vst <- dds_new@assays@data$vst

library(tidyverse)
library(tibble)
dat <- counts_vst %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  gather(key = "sample", value = "counts_vst", -gene) %>%
  left_join(meta_data, by="sample")

#library(ggplot2)
p_1 <- ggplot(data=dat,aes(x=sample,y=counts_vst,fill=condition)) +
  geom_boxplot()+
  labs(x="Sample",y="VST Normalized Read Counts")+
  ggtitle("VST Normalized Read Count by Sample") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "QCplots/boxplot-counts_vst.pdf", plot = p_1, width = 10, height = 5)
ggsave(filename = "QCplots/boxplot-counts_vst.png", plot = p_1, width = 10, height = 5, units = "in", dpi = 300)

## 1.2. counts_rlog
counts_rlog <- dds_new@assays@data$rlog
#library(tibble)
dat <- counts_rlog %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  gather(key = "sample", value = "counts_rlog", -gene) %>%
  left_join(meta_data, by="sample")

#library(ggplot2)
p_2 <- ggplot(data=dat,aes(x=sample,y=counts_rlog,fill=condition)) +
  geom_boxplot()+
  labs(x="Sample",y="rlog Normalized Read Counts")+
  ggtitle("rlog Normalized Read Count by Sample") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "QCplots/boxplot-counts_rlog.pdf", plot = p_2, width = 10, height = 5)
ggsave(filename = "QCplots/boxplot-counts_rlog.png", plot = p_2, width = 10, height = 5, units = "in", dpi = 300)


### 2. sample correlation heatmap
sampleDists <- dist(t(counts_vst))
sampleDistMatrix <- as.matrix(sampleDists)

library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
mydata_col = data.frame(condition=meta_data$condition, row.names = meta_data$sample)

library(pheatmap)
p_3 <- pheatmap(
  sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors,
  annotation_col = mydata_col,
  main="Euclidean distance between samples"
)

### make my customized function to save heatmap
save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_3, "QCplots/sample_heatmap.pdf", 5, 4.1)

save_pheatmap_png <- function(x, filename, width, height, dpi = 300) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  
  # 计算像素尺寸
  width_px <- width * dpi
  height_px <- height * dpi
  
  # 使用 png 设备保存图像
  png(filename, width = width_px, height = height_px, res = dpi, units = "px")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# 调用示例
save_pheatmap_png(p_3, "QCplots/sample_heatmap.png", 5, 4.1)



### 3. gene expression heatmap of top 1000 variable genes
counts_vst_var <- as.data.frame(counts_vst)
counts_vst_var$var <- apply(counts_vst,1,var)
counts_vst_var <- counts_vst_var %>%
  arrange(desc(var))
saveRDS(counts_vst_var, "counts_vst_var.rds", compress=T)
counts_vst_var <- counts_vst_var %>%
  slice(1:1000) %>%
  as.matrix()

p_4 <- pheatmap(counts_vst_var[,1:6],
                scale="row",
                color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
                annotation_col = mydata_col,
                main="Gene Expression (VST) of Top 1000 Variable Genes",
                show_rownames = F)
save_pheatmap_pdf(p_4, "QCplots/gene_heatmap.pdf", 5, 6)
save_pheatmap_png(p_4, "QCplots/gene_heatmap.png", 5, 6)

### 4. PCA plots
# make my customized function to plot PCA
plotPCA_vst <- function (object, assay, ntop = 1000) {
  rv <- rowVars(assay(object, assay))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  df <- cbind( as.data.frame(colData(object)), pca$x)
  # order points so extreme samples are more likely to get label
  ord <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
  df <- df[ord,]
  attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
  return(df)
}

pca.data <- plotPCA_vst(dds_new, assay="vst", ntop=1000)
saveRDS(pca.data, "pca.data.rds", compress = T)

for (i in c(1:4)) {
  for (j in c((i+1):5)) {
    percentVar <- round(attr(pca.data, "percentVar")$percentVar)
    p <- ggplot(pca.data, aes_string(x=paste0("PC",i), y=paste0("PC",j), color="condition")) +
      geom_point(size=3) +
      geom_text(aes(label = sample), vjust=-0.5, hjust=0.5) +
      xlab(paste0("PC",i,": ",percentVar[i],"% variance")) +
      ylab(paste0("PC",j,": ",percentVar[j],"% variance")) +
      theme_bw()
    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".pdf"), plot = p, width=5, height=4)
    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".png"), plot = p, width=5, height=4, units = "in", dpi = 300)
  }
}

### create folders
setwd("../")
if (!dir.exists("03-DE")){
  dir.create("03-DE")
} 
setwd("./03-DE")

bg_ID <- readRDS("Z:/ResearchHome/SharedResources/Immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/03-DE/bg_ID.rds")

### create folders
if (!dir.exists("DEfiles")){
  dir.create("DEfiles")
} 

if (!dir.exists("MAplots")){
  dir.create("MAplots")
} 

if (!dir.exists("Volcanoplots")){
  dir.create("Volcanoplots")
} 

if (!dir.exists("Heatmaps")){
  dir.create("Heatmaps")
} 

### generate FC and FDR via DESeq2
library(DESeq2)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
# read in dds object and metadata
dds_new <- readRDS("../02-QC/dds_new.rds")
meta_data <- readRDS("../02-QC/meta_data.rds")

### run DESeq to get result objects in loop
library(ggplot2)
library(readr)
library(tibble)
library(pheatmap)
counts_vst <- dds_new@assays@data$vst
dds_new_res <- DESeq(dds_new)
resultsNames(dds_new_res)
# "Intercept"              "condition_OLAH_vs_Ctrl"
saveRDS(dds_new_res, "dds_new_res.rds", compress = T)

### 1. DE results
res <- results(dds_new_res)
saveRDS(res, paste0("DEfiles/DE-",sampleB,"vs",sampleA,".rds"), compress = T)

# join gene_name
res_df <- as.data.frame(res)
res_symbol <- res_df %>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  left_join(bg_ID, by="gene")

### 2. MA plots
pdf(file = paste0("MAplots/MAplot-",sampleB,"vs",sampleA,".pdf"),
    width=7, height=5)
MAplot <- plotMA(object=res, alpha=0.05,
                 main=paste0(sampleB,"vs",sampleA),
                 colNonSig="black",
                 colSig="red",
                 colLine="blue")
dev.off()

### 3. volcano plots
vol <- EnhancedVolcano(res_symbol,
                       lab = NA,
                       x = 'log2FoldChange',
                       y = 'padj',
                       ylab = bquote(~-Log[10] ~ italic(FDR)),
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = paste0(sampleB,"vs",sampleA),
                       legendLabels=c('Not DE &\nabsolute log2FC < 1',
                                      'Not DE &\nabsolute log2FC > 1',
                                      'FDR < 0.05 &\nabsolute log2FC < 1',
                                      'FDR < 0.05 &\nabsolute log2FC > 1'),
                       legendPosition = 'right')
ggsave(filename = paste0("Volcanoplots/Volcanoplot-",sampleB,"vs",sampleA,".pdf"), plot=vol, width = 8, height = 7)
ggsave(filename = paste0("Volcanoplots/Volcanoplot-",sampleB,"vs",sampleA,".png"), plot=vol, width = 8, height = 7, units = "in", dpi=300)

vol_lab <- EnhancedVolcano(res_symbol,
                           lab = res_symbol$gene_name,
                           x = 'log2FoldChange',
                           y = 'padj',
                           ylab = bquote(~-Log[10] ~ italic(FDR)),
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           title = paste0(sampleB,"vs",sampleA),
                           legendLabels=c('Not DE &\nabsolute log2FC < 1',
                                          'Not DE &\nabsolute log2FC > 1',
                                          'FDR < 0.05 &\nabsolute log2FC < 1',
                                          'FDR < 0.05 &\nabsolute log2FC > 1'),
                           legendPosition = 'right')
ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-",sampleB,"vs",sampleA,".pdf"), plot=vol_lab, width = 8, height = 7)
ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-",sampleB,"vs",sampleA,".png"), plot=vol_lab, width = 8, height = 7, units = "in", dpi=300)








### 4. heatmaps
res_symbol <- res_symbol %>%
  mutate(change=case_when(padj < 0.05 ~ "DE",
                          TRUE ~ "Not")) %>%
  mutate(direction=case_when(change == "DE" & log2FoldChange > 0 ~ "Up",
                             change == "DE" & log2FoldChange < 0 ~ "Down",
                             TRUE ~ "Not"))
write.csv(res_symbol, paste0("DEfiles/DE-",sampleB,"vs",sampleA,".csv"))

DE_genes <- res_symbol[res_symbol$change=="DE",]$gene

counts_vst_DE <- counts_vst[rownames(counts_vst) %in% DE_genes,]
# 找出 DE_genes 中有但不在  counts_vst_DE的行名中的基因
setdiff(DE_genes, rownames(counts_vst_DE))
setdiff(rownames(counts_vst_DE),DE_genes)
rownames(counts_vst_DE)[duplicated(rownames(counts_vst_DE))]
DE_genes[duplicated(DE_genes)]   # "TEC"

mydata_col <- data.frame(condition=meta_data$condition, row.names = meta_data$sample)
heat_DE_1 <- pheatmap(counts_vst_DE,
                      scale="row",
                      color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
                      annotation_col = mydata_col,
                      main=paste0(sampleB,"vs",sampleA," FDR < 0.05"),
                      show_rownames = F)
save_pheatmap_pdf(heat_DE_1, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".pdf"), 5, 6)
save_pheatmap_png(heat_DE_1, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".png"), 5, 6)

# row names
tree_row <- data.frame(gene = rownames(counts_vst_DE[heat_DE_1$tree_row[["order"]],])) %>%
  left_join(bg_ID, by="gene")
saveRDS(tree_row, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".rds"), compress = T)
write.csv(tree_row, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".csv"))


### Create folders
setwd("Z:/ResearchHome/SharedResources/Immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/04-exclude_outliers/")
if (!dir.exists("04-Pathway")){
  dir.create("04-Pathway")
}
setwd("./04-Pathway")

if (!dir.exists("Pathways_tables")){
  dir.create("Pathways_tables")
}

if (!dir.exists("Pathways_figures")){
  dir.create("Pathways_figures")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(readr)

bg_ID <- readRDS("Z:/ResearchHome/SharedResources/Immunoinformatics/jqu/project/OLAH/analysis/16-Raghu_bulk_RNAseq/analysis/03-DE/bg_ID.rds")
bg <- bg_ID$ENTREZID
meta_data <- readRDS("../02-QC/meta_data.rds")
treatment <- levels(meta_data$treatment)

### Loops
res_symbol <- read.csv(paste0("../03-DE/DEfiles/DE-",sampleB,"vs",sampleA,".csv"))

## Inner loop for Up or Down
for (j in c("Up", "Down")){
  DEG <- res_symbol[res_symbol$direction == j,]$ENTREZID %>%
    as.character()
  
  # KEGG
  KEGG <- enrichKEGG(gene = DEG,
                     universe = bg,
                     organism ='hsa',
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     use_internal_data =FALSE)
  KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write_csv(as.data.frame(KEGG@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".KEGG.csv"))
  saveRDS(KEGG,paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".KEGG.rds"),compress=T)
  dot_KEGG <- dotplot(KEGG)
  ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".KEGG.pdf"), plot = dot_KEGG, width = 7, height = 7)  
  ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".KEGG.png"), plot = dot_KEGG, width = 7, height = 7, units = "in", dpi=300)  
  
  # GO sub-categories
  for (GO in c("CC", "MF", "BP")) {
    GO_enrich <- enrichGO(gene = DEG,
                          universe = bg,
                          OrgDb = "org.Hs.eg.db",
                          keyType = "ENTREZID",
                          ont = GO,
                          pvalueCutoff  = 0.1,
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.1,
                          readable=T)
    GO_enrich <- simplify(GO_enrich,cutoff=0.7,
                          by="p.adjust",
                          select_fun=min,
                          measure = "Wang")
    write_csv(as.data.frame(GO_enrich@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".",GO,".csv"))
    saveRDS(GO_enrich,paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".",GO,".rds"),compress=T)
    dot_GO <- dotplot(GO_enrich)
    ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".",GO,".pdf"), plot = dot_GO, width = 7, height = 7)
    ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".",GO,".png"), plot = dot_GO, width = 7, height = 7, units = "in", dpi=300)
    
  }  
}

### Do not separate Up and Down
DEG <- res_symbol[res_symbol$change == "DE",]$ENTREZID %>%
  as.character()

# KEGG
KEGG <- enrichKEGG(gene = DEG,
                   universe = bg,
                   organism ='hsa',
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1,
                   use_internal_data =FALSE)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write_csv(as.data.frame(KEGG@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".DE.KEGG.csv"))
saveRDS(KEGG,paste0("Pathways_tables/",sampleB,"vs",sampleA,".DE.KEGG.rds"),compress=T)
dot_KEGG <- dotplot(KEGG)
ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".DE.KEGG.pdf"), plot = dot_KEGG, width = 7, height = 7)  
ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".DE.KEGG.png"), plot = dot_KEGG, width = 7, height = 7, units = "in", dpi=300)  

# GO sub-categories
for (GO in c("CC", "MF", "BP")) {
  GO_enrich <- enrichGO(gene = DEG,
                        universe = bg,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID",
                        ont = GO,
                        pvalueCutoff  = 0.1,
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.1,
                        readable=T)
  GO_enrich <- simplify(GO_enrich,cutoff=0.7,
                        by="p.adjust",
                        select_fun=min,
                        measure = "Wang")
  write_csv(as.data.frame(GO_enrich@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".DE.",GO,".csv"))
  saveRDS(GO_enrich,paste0("Pathways_tables/",sampleB,"vs",sampleA,".DE.",GO,".rds"),compress=T)
  dot_GO <- dotplot(GO_enrich)
  ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".DE.",GO,".pdf"), plot = dot_GO, width = 7, height = 7)
  ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".DE.",GO,".png"), plot = dot_GO, width = 7, height = 7, units = "in", dpi=300)
}
gc()

```
