###setting directory
setwd("E:/BIOINFO Project")

if(!require("BiocManager",quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("tidyverse")

###DESeq2 Preliminary Analysis

### Loading Libraries
library("DESeq2")
library("tidyverse")
library("ggplot2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")


raw_count_data=read.csv("GSE229611_Count_Bulk_27HC_Treat.csv")
sample_info=read.csv("GSE229611_ColData.csv")

row.names(raw_count_data) <- raw_count_data[,1]
raw_count_data <- raw_count_data[,-1]

row.names(sample_info) <- sample_info[,1]

sample_info$Design <- factor(sample_info$Design)

dds_object <- DESeqDataSetFromMatrix(countData = raw_count_data, colData=sample_info, design=~Design)

dds_object

prdds_object <- DESeq(dds_object)

prdds_object

final_wo_fdr <- results(prdds_object)
final_wo_fdr

write.csv(final_wo_fdr, "Akshith_Output.csv")

final_w_fdr <- results(prdds_object, alpha = 0.1)
final_w_fdr

write.csv(final_w_fdr, "Akshith_output_FDRcutoffincluded.csv")

#MA plot
plotMA(final_w_fdr, ylim=c(-5,5))

#Volcano plot
res_all <- data.frame(final_w_fdr) %>% mutate(threshold = padj < 0.05)              
# Create the volcano plot
ggplot(res_all) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), ymin = -3, ymax = 10) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25))) +
        scale_y_continuous(limits=c(-3,10)) +
        scale_x_continuous(limits=c(-10,10))
