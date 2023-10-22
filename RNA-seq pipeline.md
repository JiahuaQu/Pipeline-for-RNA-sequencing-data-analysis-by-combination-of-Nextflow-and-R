RNA-seq pipeline

# Upstream

```bash
1.	cd ~
2.	mkdir 01-Nextflow
3.	cd ~/01-Nextflow
4.	module load Sali anaconda
5.	conda create -n NextFlow-nf-core
6.	source activate NextFlow-nf-core   
7.	python -V     # Python 3.8.10

1.	# install nextflow
2.	conda install nextflow -y
3.	conda list | grep nextflow   # version: 22.10.6
4.	 
5.	# install nf-core
6.	conda search nf-core
7.	conda install nf-core -y
8.	conda list | grep nf-core   # 2.7.2
9.	
10.	# install rnaseq pipeline
11.	nf-core download nf-core/rnaseq
12.	# select version：3.10.1 
13.	# download singularity

1.	cd ~/01-Nextflow
2.	mkdir -p reference_genome/mouse/Ensembl/GRCm38
3.	cd reference_genome/mouse/Ensembl/GRCm38
4.	nohup aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ . --exclude "*" --include "genes.gtf" &
5.	nohup aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/ ./ &

1.	cd ~/01-Nextflow
2.	mkdir data

1.	cd ~/01-Nextflow/data
2.	wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py
3.	chmod a+x fastq_dir_to_samplesheet.py
4.	./fastq_dir_to_samplesheet.py -r1 _R1_001.fastq.gz ~/data/fastq samplesheet.csv

1.	cd ~/01-Nextflow
2.	mkdir analysis
3.	mkdir output
4.	cd ~/01-Nextflow/analysis
5.	vim nextflow.config
```

```groovy
1.	// Global default params, used in configs
2.	params {
3.	 
4.	    // Input options
5.	    input = '~/data/fastq/samplesheet.csv'
6.	 
7.	    // References
8.	    fasta = '~/reference_genome/mouse/Ensembl/GRCm38/genome.fa'
9.	    gtf = '~/reference_genome/mouse/Ensembl/GRCm38/genes.gtf'
10.	    gencode = false
11.	    save_reference = true
12.	     
13.	    // QC
14.	    skip_bigwig = true
15.	    skip_stringtie = true
16.	    skip_preseq = true
17.	    skip_qualimap = true
18.	    deseq2_vst = true
19.	
20.	    // Boilerplate options
21.	    outdir = '~/output'
22.	    email = "jiahua.qu@hotmail.com"
23.	 
24.	    // Max resource options
25.	    // Defaults only, expecting to be overwritten
26.	    max_memory = '500.GB'
27.	    max_cpus = 12
28.	    max_time = '240.h'
29.	}
```

```bash
1.	cd ~/01-Nextflow/analysis
2.	vim my_script.sh
```

```bash
1.	#! /usr/bin/env bash
2.	#
3.	#$ -S /bin/bash
4.	#$ -cwd
5.	#$ -r y
6.	#$ -j y
7.	#$ -pe smp 6     
8.	#$ -l mem_free=45G
9.	#$ -l scratch=270G   
10.	#$ -l h_rt=24:00:00
11.	 
12.	module load Sali anaconda
13.	source activate NextFlow-nf-core   
14.	 
15.	# set path
16.	log=~/01-Nextflow/analysis/my.log
17.	workflow=~/01-Nextflow/nf-core-rnaseq-3.10.1/workflow
18.	profile=singularity
19.	config=~/01-Nextflow/analysis/nextflow.config
20.	 
21.	# run pipeline
22.	nextflow -log $log run $workflow -profile $profile -resume -c $config
23.	  
24.	conda deactivate
25.	 
26.	[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" 
```

```bash
1.	qsub my_script.sh
```

# Midstream

```R
1.	### create a folder
2.	setwd("~/")
3.	if (!file.exists("02-QC")){
4.	  dir.create("02-QC")
5.	} 
6.	setwd("~/02-QC") 
7.	
8.	### read in rds
9.	# use the salmon.merged.gene_counts.rds
10.	gene_counts <- readRDS("~/01-Nextflow/output/star_salmon/salmon.merged.gene_counts.rds")   
11.	
12.	# extract and process count matrix
13.	counts <- gene_counts@assays@data$counts
14.	# round non-integer to integer for loading into DESeq2
15.	counts_round <- round(counts, digits = 0)
16.	# remove rows whose values are all zero
17.	library(dplyr)
18.	counts_round <- counts_round %>%
19.	  filter(rowSums(.) > 0)
20.	 
21.	# extract meta_data
22.	library(tidyr)
23.	meta_data <- data.frame(sample=gene_counts@colData@rownames) %>%
24.	  separate(sample, into=c("treatment","genotype","num"), sep = "\\.", remove = FALSE) %>%
25.	  unite(col = "condition", treatment,genotype, remove = FALSE)
26.	meta_data$condition <- factor(meta_data$condition)
27.	meta_data$treatment <- factor(meta_data$treatment)
28.	meta_data$genotype <- factor(meta_data$genotype, levels = c("con", "mut"))
29.	saveRDS(meta_data, "meta_data.rds", compress = T)

1.	library(DESeq2)
2.	# Load in Salmon output
3.	dds_new <- DESeqDataSetFromMatrix(countData = counts_round,
4.	                                  colData = meta_data,
5.	                                  design= ~ condition)
6.	 
7.	# Normalization
8.	dds_new <- estimateSizeFactors(dds_new)
9.	dds_rlog <- rlog(dds_new)
10.	dds_vst <- varianceStabilizingTransformation(dds_new)
11.	assay(dds_new, "rlog") <- assay(dds_rlog)
12.	assay(dds_new, "vst") <- assay(dds_vst)
13.	saveRDS(dds_new, "dds_new.rds", compress = T)
14.	saveRDS(dds_rlog, "dds_rlog.rds", compress = T)
15.	saveRDS(dds_vst, "dds_vst.rds", compress = T)

1.	### create a folder
2.	setwd("~/02-QC")
3.	if (!file.exists("QCplots")){
4.	  dir.create("QCplots")
5.	}
6.	 
7.	### 1. boxplot
8.	library(DESeq2)
9.	## 1.1. counts_vst
10.	counts_vst <- dds_new@assays@data$vst
11.	 
12.	library(tibble)
13.	dat <- counts_vst %>%
14.	  as.data.frame() %>%
15.	  rownames_to_column(var="gene") %>%
16.	  gather(key = "sample", value = "counts_vst", -gene) %>%
17.	  left_join(meta_data, by="sample")
18.	 
19.	library(ggplot2)
20.	p_1 <- ggplot(data=dat,aes(x=sample,y=counts_vst,fill=condition)) +
21.	  geom_boxplot()+
22.	  labs(x="Sample",y="VST Normalized Read Counts")+
23.	  ggtitle("VST Normalized Read Count by Sample") +
24.	  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
25.	ggsave(filename = "QCplots/boxplot-counts_vst.pdf", plot = p_1, width = 10, height = 5)
26.	 
27.	## 1.2. counts_rlog
28.	counts_rlog <- dds_new@assays@data$rlog
29.	 
30.	library(tibble)
31.	dat <- counts_rlog %>%
32.	  as.data.frame() %>%
33.	  rownames_to_column(var="gene") %>%
34.	  gather(key = "sample", value = "counts_rlog", -gene) %>%
35.	  left_join(meta_data, by="sample")
36.	 
37.	library(ggplot2)
38.	p_2 <- ggplot(data=dat,aes(x=sample,y=counts_rlog,fill=condition)) +
39.	  geom_boxplot()+
40.	  labs(x="Sample",y="rlog Normalized Read Counts")+
41.	  ggtitle("rlog Normalized Read Count by Sample") +
42.	  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
43.	ggsave(filename = "QCplots/boxplot-counts_rlog.pdf", plot = p_2, width = 10, height = 5)
44.	 
45.	### 2. sample correlation heatmap
46.	sampleDists <- dist(t(counts_vst))
47.	sampleDistMatrix <- as.matrix(sampleDists)
48.	 
49.	library(RColorBrewer)
50.	colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
51.	mydata_col = data.frame(condition=meta_data$condition, row.names = meta_data$sample)
52.	
53.	library(pheatmap)
54.	p_3 <- pheatmap(
55.	  sampleDistMatrix,
56.	  clustering_distance_rows=sampleDists,
57.	  clustering_distance_cols=sampleDists,
58.	  col=colors,
59.	  annotation_col = mydata_col,
60.	  main="Euclidean distance between samples"
61.	)
62.	
63.	### make my customized function to save heatmap
64.	save_pheatmap_pdf <- function(x, filename, width, height) {
65.	  stopifnot(!missing(x))
66.	  stopifnot(!missing(filename))
67.	  pdf(filename, width=width, height=height)
68.	  grid::grid.newpage()
69.	  grid::grid.draw(x$gtable)
70.	  dev.off()
71.	}
72.	save_pheatmap_pdf(p_3, "QCplots/sample_heatmap.pdf", 7, 6.5)
73.	
74.	### 3. gene expression heatmap of top 1000 variable genes
75.	counts_vst_var <- as.data.frame(counts_vst)
76.	counts_vst_var$var <- apply(counts_vst,1,var)
77.	counts_vst_var <- counts_vst_var %>%
78.	saveRDS(counts_vst_var, "counts_vst_var.rds", compress=T)
79.	counts_vst_var <- counts_vst_var %>%
80.	  slice(1:1000) %>%
81.	  as.matrix()
82.	 
83.	p_4 <- pheatmap(counts_vst_var[,1:24],
84.	              scale="row",
85.	              color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
86.	              annotation_col = mydata_col,
87.	              main="Gene Expression (VST) of Top 1000 Variable Genes",
88.	              show_rownames = F)
89.	save_pheatmap_pdf(p_4, "QCplots/gene_heatmap.pdf", 7, 6.5)
90.	 
91.	### 4. PCA plots
92.	# make my customized function to plot PCA
93.	plotPCA_vst <- function (object, assay, ntop = 1000) {
94.	  rv <- rowVars(assay(object, assay))
95.	  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
96.	  pca <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
97.	  percentVar <- pca$sdev^2/sum(pca$sdev^2)
98.	  df <- cbind( as.data.frame(colData(object)), pca$x)
99.	  # order points so extreme samples are more likely to get label
100.	  ord <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
101.	  df <- df[ord,]
102.	  attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
103.	  return(df)
104.	}
105.	 
106.	pca.data <- plotPCA_vst(dds_new, assay="vst", ntop=1000)
107.	saveRDS(pca.data, "pca.data.rds", compress = T)
108.	 
109.	for (i in c(1:4)) {
110.	  for (j in c((i+1):5)) {
111.	    percentVar <- round(attr(pca.data, "percentVar")$percentVar)
112.	    p <- ggplot(pca.data, aes_string(x=paste0("PC",i), y=paste0("PC",j), color="condition")) +
113.	      geom_point(size=3) +
114.	      xlab(paste0("PC",i,": ",percentVar[i],"% variance")) +
115.	      ylab(paste0("PC",j,": ",percentVar[j],"% variance")) +
116.	      theme_bw()
117.	    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".pdf"), plot = p, width=6, height=5)
118.	  }
119.	}
120.	
121.	### 5. t-SNE plots
122.	library(Rtsne)
123.	# make my customized function to plot t-SNE
124.	plottsne_vst <- function (object, assay, perpl, ntop = 1000) {
125.	  rv <- rowVars(assay(object, assay))
126.	  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
127.	  set.seed(123)
128.	  tsne <- Rtsne(t(assay(object, assay)[select, ]), perplexity=perpl)
129.	  df <- as.data.frame(tsne$Y)
130.	  names(df) <- c("Dimention1","Dimention2")
131.	  df <- cbind(df, meta_data)
132.	  return(df)
133.	}
134.	 
135.	for (i in 1:6) {
136.	  tsne.data <- plottsne_vst(dds_new, assay="vst", perpl=i, ntop=1000)
137.	  saveRDS(tsne.data, paste0("tsne.data-perplexity_",i,".rds"), compress = T)
138.	  
139.	  p <- ggplot(tsne.data, aes_string(x="Dimention1", y="Dimention2", color="condition")) +
140.	    geom_point(size=3) +
141.	    xlab("Dimention 1") +
142.	    ylab("Dimention 2") +
143.	    theme_bw()
144.	  ggsave(filename = paste0("QCplots/tSNEplot-perplexity_",i,".pdf"), plot = p, width=6, height=5)
145.	}

1.	### create folders
2.	setwd("~/")
3.	if (!file.exists("03-DE")){
4.	  dir.create("03-DE")
5.	} 
6.	setwd("~/03-DE")
7.	
8.	### read in gene information obtained from Nextflow
9.	library(readr)
10.	bg_ID <- read_delim("~/01-Nextflow/output/star_salmon/salmon_tx2gene.tsv", 
11.	                    delim = "\t", escape_double = FALSE, 
12.	                    col_names = FALSE, trim_ws = TRUE)
13.	names(bg_ID) <- c("transcript", "gene", "gene_name")
14.	bg_ID <- bg_ID[,-1] %>% unique()
15.	
16.	library(clusterProfiler)
17.	library(org.Mm.eg.db)
18.	#keytypes(org.Mm.eg.db)
19.	bg_ENTREZID <- bitr(bg_ID$gene_name,
20.	                    fromType = "SYMBOL",
21.	                    toType = "ENTREZID",
22.	                    OrgDb = "org.Mm.eg.db")
23.	bg_ID <- left_join(bg_ID, bg_ENTREZID, by=c("gene_name"="SYMBOL"))
24.	saveRDS(bg_ID, "bg_ID.rds", compress = T)

1.	### create folders
2.	setwd("~/03-DE")
3.	
4.	if (!file.exists("DEfiles")){
5.	  dir.create("DEfiles")
6.	} 
7.	 
8.	if (!file.exists("MAplots")){
9.	  dir.create("MAplots")
10.	} 
11.	
12.	if (!file.exists("Volcanoplots")){
13.	  dir.create("Volcanoplots")
14.	} 
15.	
16.	if (!file.exists("Heatmaps")){
17.	  dir.create("Heatmaps")
18.	} 
19.	
20.	### generate FC and FDR via DESeq2
21.	library(DESeq2)
22.	library(EnhancedVolcano)
23.	# read in dds object and metadata
24.	dds_new <- readRDS("~/02-QC/dds_new.rds")
25.	meta_data <- readRDS("~/02-QC/meta_data.rds")
26.	
27.	### run DESeq to get result objects in loop
28.	library(ggplot2)
29.	library(readr)
30.	library(tibble)
31.	library(pheatmap)
32.	counts_vst <- dds_new@assays@data$vst
33.	dds_new_res <- DESeq(dds_new)
34.	saveRDS(dds_new_res, "dds_new_res.rds", compress = T)
35.	 
36.	for (i in levels(meta_data$treatment)) {
37.	
38.	  ### 1. DE results
39.	  sampleA <- paste(i, "con", sep = "_")
40.	  sampleB <- paste(i, "mut", sep = "_")
41.	  condition <- "condition"
42.	  contrastV <- c(condition, sampleA, sampleB)
43.	  res <- results(dds_new_res,  contrast=contrastV)
44.	  saveRDS(res, paste0("DEfiles/DE-",sampleB,"vs",sampleA,".rds"), compress = T)
45.	  
46.	  # join gene_name
47.	  res_df <- as.data.frame(res)
48.	  res_symbol <- res_df %>%
49.	    as.data.frame() %>%
50.	    rownames_to_column(var="gene") %>%
51.	    left_join(bg_ID, by="gene")
52.	   
53.	  ### 2. MA plots
54.	  pdf(file = paste0("MAplots/MAplot-",sampleB,"vs",sampleA,".pdf"),
55.	      width=7, height=5)
56.	  MAplot <- plotMA(object=res, alpha=0.05,
57.	                   main=paste0(sampleB,"vs",sampleA),
58.	                   colNonSig="black",
59.	                   colSig="red",
60.	                   colLine="blue")
61.	  dev.off()
62.	  
63.	  ### 3. volcano plots
64.	  vol <- EnhancedVolcano(res_symbol,
65.	                  lab = NA,
66.	                  x = 'log2FoldChange',
67.	                  y = 'padj',
68.	                  ylab = bquote(~-Log[10] ~ italic(FDR)),
69.	                  pCutoff = 0.05,
70.	                  FCcutoff = 2,
71.	                  title = paste0(sampleB,"vs",sampleA),
72.	                  legendLabels=c('Not DE &\nabsolute FC < 2',
73.	                                 'Not DE &\nabsolute FC > 2',
74.	                                 'FDR < 0.05 &\nabsolute FC < 2',
75.	                                 'FDR < 0.05 &\nabsolute FC > 2'),
76.	                  legendPosition = 'right')
77.	  ggsave(filename = paste0("Volcanoplots/Volcanoplot-",sampleB,"vs",sampleA,".pdf"), plot=vol, width = 7, height = 7)
78.	   
79.	  vol_lab <- EnhancedVolcano(res_symbol,
80.	                             lab = res_symbol$gene_name,
81.	                             x = 'log2FoldChange',
82.	                             y = 'padj',
83.	                             ylab = bquote(~-Log[10] ~ italic(FDR)),
84.	                             pCutoff = 0.05,
85.	                             FCcutoff = 2,
86.	                             title = paste0(sampleB,"vs",sampleA),
87.	                             legendLabels=c('Not DE &\nabsolute FC < 2',
88.	                                            'Not DE &\nabsolute FC > 2',
89.	                                            'FDR < 0.05 &\nabsolute FC < 2',
90.	                                            'FDR < 0.05 &\nabsolute FC > 2'),
91.	                             legendPosition = 'right')
92.	  ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-",sampleB,"vs",sampleA,".pdf"), plot=vol_lab, width = 7, height = 7)
93.	  
94.	  ### 4. heatmaps
95.	  res_symbol <- res_symbol %>%
96.	    mutate(change=case_when(padj < 0.05 ~ "DE",
97.	                            TRUE ~ "Not")) %>%
98.	    mutate(direction=case_when(change == "DE" & log2FoldChange > 0 ~ "Up",
99.	                               change == "DE" & log2FoldChange < 0 ~ "Down",
100.	                               TRUE ~ "Not"))
101.	  write.csv(res_symbol, paste0("DEfiles/DE-",sampleB,"vs",sampleA,".csv"))
102.	  
103.	  counts_vst_DE <- counts_vst[res_symbol$change=="DE",]
104.	  mydata_col <- data.frame(condition=meta_data$condition, row.names = meta_data$sample)
105.	  heat_DE_1 <- pheatmap(counts_vst_DE,
106.	                        scale="row",
107.	                        color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
108.	                        annotation_col = mydata_col,
109.	                        main=paste0(sampleB,"vs",sampleA," FDR < 0.05"),
110.	                        show_rownames = F)
111.	  save_pheatmap_pdf(heat_DE_1, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".pdf"), 7, 6.5)
112.	   
113.	  # row names
114.	  tree_row <- data.frame(gene = rownames(counts_vst_DE[heat_DE_1$tree_row[["order"]],])) %>%
115.	    left_join(bg_ID, by="gene")
116.	  saveRDS(tree_row, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".rds"), compress = T)
117.	  write.csv(tree_row, paste0("Heatmaps/Heatmap-",sampleB,"vs",sampleA,".csv"))
118.	   
119.	  meta_data_2 <- meta_data[meta_data$treatment==i,]
120.	  mydata_col_2 = data.frame(condition=meta_data_2$condition, row.names = meta_data_2$sample)
121.	  counts_vst_DE_2 <- counts_vst_DE[,meta_data_2[,1]]
122.	  heat_DE_2 <- pheatmap(counts_vst_DE_2,
123.	                        scale="row",
124.	                        color = colorRampPalette(c("darkblue", "white", "darkred"))(1000),
125.	                        annotation_col = mydata_col_2,
126.	                        main=paste0(sampleB,"vs",sampleA," FDR < 0.05"),
127.	                        show_rownames = F)
128.	  save_pheatmap_pdf(heat_DE_2, paste0("Heatmaps/Heatmap2-",sampleB,"vs",sampleA,".pdf"), 7, 6.5)
129.	   
130.	  # row names
131.	  tree_row_2 <- data.frame(gene = rownames(counts_vst_DE[heat_DE_2$tree_row[["order"]],])) %>%
132.	    left_join(bg_ID, by="gene")
133.	  saveRDS(tree_row_2, paste0("Heatmaps/Heatmap2-",sampleB,"vs",sampleA,".rds"), compress = T)
134.	  write.csv(tree_row_2, paste0("Heatmaps/Heatmap2-",sampleB,"vs",sampleA,".csv"))
135.	}
```

# Downstream

```R
1.	### Create folders
2.	setwd("~")
3.	if (!file.exists("04-Pathway")){
4.	  dir.create("04-Pathway")
5.	}
6.	setwd("~/04-Pathway")
7.	
8.	if (!file.exists("Pathways_tables")){
9.	  dir.create("Pathways_tables")
10.	}
11.	
12.	if (!file.exists("Pathways_figures")){
13.	  dir.create("Pathways_figures")
14.	}
15.	
16.	library(clusterProfiler)
17.	library(org.Mm.eg.db)
18.	library(ggplot2)
19.	library(readr)
20.	
21.	bg_ID <- readRDS("../03-DE/bg_ID.rds")
22.	bg <- bg_ID$ENTREZID
23.	meta_data <- readRDS("../02-QC/meta_data.rds")
24.	treatment <- levels(meta_data$treatment)
25.	
26.	### Loops
27.	for (i in treatment) {
28.	  sampleA <- paste(i, "con", sep = "_")
29.	  sampleB <- paste(i, "mut", sep = "_")
30.	  # Read in DE table
31.	  res_symbol <- read.csv(paste0("../03-DE/DEfiles/DE-",sampleB,"vs",sampleA,".csv"))
32.	  
33.	  ## Inner loop for Up or Down
34.	  for (j in c("Up", "Down")){
35.	    DEG <- res_symbol[res_symbol$direction == j,]$ENTREZID %>%
36.	      as.character()
37.	    
38.	    # KEGG
39.	    KEGG <- enrichKEGG(gene = DEG,
40.	                       universe = bg,
41.	                       organism ='mmu',
42.	                       pvalueCutoff = 0.1,
43.	                       qvalueCutoff = 0.1,
44.	                       use_internal_data =FALSE)
45.	    KEGG <- setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
46.	    write_csv(as.data.frame(KEGG@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".KEGG.csv"))
47.	    saveRDS(KEGG,paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".KEGG.rds"),compress=T)
48.	    dot_KEGG <- dotplot(KEGG)
49.	    ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".KEGG.pdf"), plot = dot_KEGG, width = 7, height = 7)  
50.	    
51.	    # GO sub-categories
52.	    for (GO in c("CC", "MF", "BP")) {
53.	      GO_enrich <- enrichGO(gene = DEG,
54.	                            universe = bg,
55.	                            OrgDb = "org.Mm.eg.db",
56.	                            keyType = "ENTREZID",
57.	                            ont = GO,
58.	                            pvalueCutoff  = 0.1,
59.	                            pAdjustMethod = "BH",
60.	                            qvalueCutoff  = 0.1,
61.	                            readable=T)
62.	      GO_enrich <- simplify(GO_enrich,cutoff=0.7,
63.	                            by="p.adjust",
64.	                            select_fun=min,
65.	                            measure = "Wang")
66.	      write_csv(as.data.frame(GO_enrich@result),paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".",GO,".csv"))
67.	      saveRDS(GO_enrich,paste0("Pathways_tables/",sampleB,"vs",sampleA,".",j,".",GO,".rds"),compress=T)
68.	      dot_GO <- dotplot(GO_enrich)
69.	      ggsave(filename = paste0("Pathways_figures/",sampleB,"vs",sampleA,".",j,".",GO,".pdf"), plot = dot_GO, width = 7, height = 7)
70.	      }  
71.	  } 
72.	}

```

