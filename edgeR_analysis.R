# Reading in file with NSEMBL gene IDs and corresponding gene symbols + assign column names
gene_symbol_mapping <- read.table("ENSG_ID2Name.txt", header = FALSE, stringsAsFactors = FALSE, row.names=1)
colnames(gene_symbol_mapping) <- c("gene_name")

# Reading in .cvs file with read counts 
raw_count_data = read.table("gene_read_counts_table_all_final_assignment.tsv", sep =' ', header = TRUE, stringsAsFactors = TRUE, row.names = 1)

# Trim last 5 rows containing other information
raw_count_data2 <- head(raw_count_data, -5)

# Filtering out reads where less than 25% of samples have counts >25
quant <- apply(raw_count_data2,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
filtered_raw_count_data <- raw_count_data2[keep,]

# Load edgeR 
library(edgeR)

# Make class labels
class <- factor( c( rep("Untreated",3), rep("Dex",3) ))

# Getting the gene names from gene_symbol_mapping for the genes in filtered_raw_count_data variable
genes <- rownames(filtered_raw_count_data)
gene_names <- gene_symbol_mapping[genes,1]

# Making a DGElist object with data
assignment_DGElist <- DGEList(counts = filtered_raw_count_data, genes = genes, group = class)

# Normalize counts
normalized_DGElist <- calcNormFactors(assignment_DGElist)

# Create a multi-dimensional scaling plot 
MDS_plot <- plotMDS(normalized_DGElist, main = "Multidimensional Scaling Plot for Treated (Dex) and Untreated Samples", cex.main = 1.0)

# Export MDS plot 
dev.copy(png, "MDS_plot.png", width = 1000, height = 800, res = 150)
dev.off() 

# Estimate dispersion
normalized_DGElist <- estimateCommonDisp(normalized_DGElist , verbose=TRUE)

# Perform differential expression testing 
differential_testing <- exactTest(normalized_DGElist, pair=c("Untreated", "Dex"))

# Printing out top 10 most differentially expressed genes by absolute log2 FC and getting the gene symbols for them
top10 <- topTags(differential_testing, n = 10, sort.by = "logFC")$table
top10_gene_symbols <- gene_symbol_mapping[top10$genes, 1]
write.table(top10_gene_symbols, file = "top10_gene_symbols.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

