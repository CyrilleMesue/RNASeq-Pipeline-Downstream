#!/usr/bin/env Rscript
require(optparse)

#Parse arguments from command line
options <- list(
  make_option(c("--deseq2_results_path"), action = "store", type = "character", help="Path to the input data. Output from Deseq2"),
  make_option(c("--outdir"), action = "store", default = "data", type = "character", help="Output directory."),
  make_option(c("--geneNamingTYpe"), action = "store", default = "ALIAS", help="Which gene naming convention is used. e.g ALIAS, ENSEMBL, etc"),
  make_option(c("--kegg_organism_db"), action = "store", default = "org.Hs.eg.db", help="Kegg organism database."),
  make_option(c("--kegg_organism"), action = "store", default = "hsa", help="Kegg code for given organism.")
)
arguments <- parse_args(OptionParser(option_list = options))

# Check if mandatory arguments are provided
if (is.null(arguments$deseq2_results_path)) {
  stop("Error: Mandatory arguments --deseq2_results_path is required.")
}

# create variable handles for parsed parameters
deseq2_results_path <- arguments$deseq2_results_path
outdir <- arguments$outdir
geneNamingTYpe <- arguments$geneNamingTYpe
kegg_organism <- arguments$kegg_organism
kegg_organism_db <- arguments$kegg_organism_db

# load required packages
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(kegg_organism_db, character.only = TRUE)
library('cowplot')


# # set working directory
# setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

# create output directories
dir.create(paste0(outdir, "/gsea"))
dir.create(paste0(outdir, "/kegg"))

#############################
# Load and Prepare Input Data
#############################

# reading in data from deseq2
df = read.csv(deseq2_results_path, header=TRUE, sep = "\t")

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GeneName

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#############################
# Gene Set Enrichment
#############################
# perform gene set enrichment analysis
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = geneNamingTYpe,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = kegg_organism_db, 
             pAdjustMethod = "none")

# save enrichment results
write.table(gse, file = paste0(outdir,"/gsea/gene-set-enrichment-results.tsv"), row.names = F, sep="\t")

###########################################
# Gene Set Enrichment Results Visualization
###########################################
# import the dose package
require(DOSE)

# Dot plot of enriched terms
dotplot1 <- dotplot(gse, showCategory=20, split=".sign", font.size =7, color = "p.adjust",label_format = 100)
ggsave(filename = paste0(outdir, "/gsea/gseGO-dotplot-enriched-terms.png"), plot = dotplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/gsea/gseGO-dotplot-enriched-terms.pdf"), plot = dotplot1 , device = "pdf", width = 8, height = 6)

# Dot Plot of activated & Suppressed GO terms
dotplot2 <- dotplot(gse, showCategory=20, split=".sign", font.size =7, label_format = 100) + facet_grid(.~.sign)
ggsave(filename = paste0(outdir, "/gsea/gseGO-dotplot-activated-suppressed-go-terms.png"), plot = dotplot2 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/gsea/gseGO-dotplot-activated-suppressed-go-terms.pdf"), plot = dotplot2 , device = "pdf", width = 8, height = 6)

# Gene-Concept Network plot of enriched terms
cnetplot1 <- cnetplot(gse, categorySize="pvalue", foldChange=gene_list, font.size =7)
ggsave(filename = paste0(outdir, "/gsea/gseGO-cnetplot-enriched-terms.png"), plot = cnetplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/gsea/gseGO-cnetplot-enriched-terms.pdf"), plot = cnetplot1 , device = "pdf", width = 8, height = 6)

# Gene-Concept Network plot of enriched terms Circular
cnetplot2 <- cnetplot(gse, foldChange=gene_list, circular = TRUE, colorEdge = TRUE, font.size =5)
ggsave(filename = paste0(outdir, "/gsea/gseGO-cnetplot-enriched-terms-circular.png"), plot = cnetplot2 , bg = "white",device = "png", width = 16, height = 8)
ggsave(filename = paste0(outdir, "/gsea/gseGO-cnetplot-enriched-terms-circular.pdf"), plot = cnetplot2 , device = "pdf", width = 16, height = 8)

# Heatmap plot of enriched terms. default (A), foldChange=geneList (B) 
p1 <- heatplot(gse, showCategory=5, label_format = 50)
p2 <- heatplot(gse, foldChange=gene_list, showCategory=5, label_format = 50)
heatplot1 <- cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
ggsave(filename = paste0(outdir, "/gsea/gseGO-heatplot-enriched-terms.png"), plot = heatplot1, bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/gsea/gseGO-heatplot-enriched-terms.pdf"), plot = heatplot1 , device = "pdf", width = 8, height = 6)

# Enrichment Map
x2 <- pairwise_termsim(gse)
emapplot1 <- emapplot(x2, showCategory = 20)
ggsave(filename = paste0(outdir, "/gsea/gseGO-enrichment-map.png"), plot = emapplot1 , bg = "white",device = "png", width = 14, height = 8)
ggsave(filename = paste0(outdir, "/gsea/gseGO-enrichment-map.pdf"), plot = emapplot1 , device = "pdf", width = 14, height = 8)

# Ridgeplot for gene set enrichment analysis.
ridgeplot1 <- ridgeplot(gse, label_format = 100, showCategory =20) + labs(x = "enrichment distribution")
ggsave(filename = paste0(outdir, "/gsea/gseGO-ridgeplot-gsea.png"), plot = ridgeplot1 , bg = "white",device = "png", width = 14, height = 8)
ggsave(filename = paste0(outdir, "/gsea/gseGO-ridgeplot-gsea.pdf"), plot = ridgeplot1 , device = "pdf", width = 14, height = 8)

# gseaplot for GSEA result(by = "runningScore"). by = "runningScore" (A), by = "preranked" (B), default (C) 
gseaplot1 <- gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
ggsave(filename = paste0(outdir, "/gsea/gseGO-gseaplot-firstpathway.png"), plot = gseaplot1 , bg = "white",device = "png", width = 16, height = 8)
ggsave(filename = paste0(outdir, "/gsea/gseGO-gseaplot-firstpathway.pdf"), plot = gseaplot1 , device = "pdf", width = 14, height = 8)

# Pmcplot of enrichment analysis (pubmed trend of enriched terms)
terms <- gse$Description[1:5]
terms <- gsub("(.{50})", "\\1\n", terms)  
pmcplot1 <- pmcplot(terms, 2010:2023)
ggsave(filename = paste0(outdir, "/gsea/gseGO-pubmed-trend-of-enriched-terms-number_and_proportion.png"), plot = pmcplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/gsea/gseGO-pubmed-trend-of-enriched-terms-number_and_proportion.pdf"), plot = pmcplot1 , device = "pdf", width = 8, height = 6)

################################################################################
# KEGG Gene Set Enrichment Analysis
################################################################################
# Convert gene IDs for gseKEGG function. We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = geneNamingTYpe, toType = "ENTREZID", OrgDb=kegg_organism_db)
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c(geneNamingTYpe)]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$GeneName %in% dedup_ids[[geneNamingTYpe]],]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# create a gseKEGG object

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

# save kegg pathway analysis results
write.table(kk2, file = paste0(outdir,"/kegg/kegg-analysis-results.tsv"), row.names = F, sep="\t")

# Dot plot of enriched terms kegg
dotplot1 <- dotplot(kk2, showCategory=20, split=".sign", font.size =7, label_format = 100)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-dotplot-enriched-terms.png"), plot = dotplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-dotplot-enriched-terms.pdf"), plot = dotplot1 , device = "pdf", width = 8, height = 6)

# Dot Plot of activated & Suppressed GO terms kegg
dotplot2 <- dotplot(kk2, showCategory=20, split=".sign", font.size =7, label_format = 100) + facet_grid(.~.sign)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-dotplot-activated-suppressed-go-terms.png"), plot = dotplot2 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-dotplot-activated-suppressed-go-terms.pdf"), plot = dotplot2 , device = "pdf", width = 8, height = 6)

# Gene-Concept Network plot of enriched terms kegg
cnetplot1 <- cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, font.size =7)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-cnetplot-enriched-terms.png"), plot = cnetplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-cnetplot-enriched-terms.pdf"), plot = cnetplot1 , device = "pdf", width = 8, height = 6)

# Gene-Concept Network plot of enriched terms Circular kegg
cnetplot2 <- cnetplot(kk2, foldChange=gene_list, circular = TRUE, colorEdge = TRUE, font.size =5, showCategory=3)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-cnetplot-enriched-terms-circular.png"), plot = cnetplot2 , bg = "white",device = "png", width = 12, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-cnetplot-enriched-terms-circular.pdf"), plot = cnetplot2 , device = "pdf", width = 12, height = 6)

# Heatmap plot of enriched terms. default (A), foldChange=geneList (B) kegg
p1 <- heatplot(kk2, showCategory=3, label_format = 50)
p2 <- heatplot(kk2, foldChange=gene_list, showCategory=3, label_format = 50)
heatplot1 <- cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
heatplot1
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-heatplot-enriched-terms.png"), plot = heatplot1, bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-heatplot-enriched-terms.pdf"), plot = heatplot1 , device = "pdf", width = 8, height = 6)

# Enrichment Map kegg
x2 <- pairwise_termsim(kk2)
emapplot <- emapplot(x2, showCategory = 20)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-enrichment-map.png"), plot = emapplot1 , bg = "white",device = "png", width = 14, height = 8)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-enrichment-map.pdf"), plot = emapplot1 , device = "pdf", width = 14, height = 8)

# Ridgeplot for gene set enrichment analysis. kegg
ridgeplot1 <- ridgeplot(kk2, label_format = 50, showCategory =20) + labs(x = "enrichment distribution")
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-ridgeplot-gsea.png"), plot = ridgeplot1 , bg = "white",device = "png", width = 14, height = 8)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-ridgeplot-gsea.pdf"), plot = ridgeplot1 , device = "pdf", width = 14, height = 8)

# gseaplot for GSEA result(by = "runningScore"). by = "runningScore" (A), by = "preranked" (B), default (C) kegg
gseaplot1 <- gseaplot(kk2, by = "all", title = gse$Description[1], geneSetID = 1)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-gseaplot-firstpathway.png"), plot = gseaplot1 , bg = "white",device = "png", width = 16, height = 8)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-gseaplot-firstpathway.pdf"), plot = gseaplot1 , device = "pdf", width = 14, height = 8)

# Pmcplot of enrichment analysis (pubmed trend of enriched terms)
terms <- kk2$Description[1:5]
terms <- gsub("(.{25})", "\\1\n", terms)  
pmcplot1 <- pmcplot(terms, 2010:2023)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-pubmed-trend-of-enriched-terms-number_and_proportion.png"), plot = pmcplot1 , bg = "white",device = "png", width = 8, height = 6)
ggsave(filename = paste0(outdir, "/kegg/gseKEGG-pubmed-trend-of-enriched-terms-number_and_proportion.pdf"), plot = pmcplot1 , device = "pdf", width = 8, height = 6)