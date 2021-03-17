# Monil Gandhi - Biologist

#File Description - 1) Creates a heatmap of the top 150 differentially expressed genes.
#                   2) Comparison of the DAVID functional annotation clustering data to the reference paper data

# install packages and load libraries
install.packages("pacman")
require(pacman)  
library(pacman)  
library(reshape2)
library(grid)
library(gridExtra)
library(magrittr)
library(tibble)

#run pacman for downloading the necessary libraries
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggpubr, grid, gridExtra, pheatmap, gdata, magrittr) 


# read the CSV's for the different developmental phases into a dataframe
fpkm_matrix <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", sep = "\t")
gene_table_P0_1 <-  read.csv("/projectnb/bf528/users/van-gogh/project_2/programmer/P0_1_cufflinks/genes.fpkm_tracking", sep = "\t")
gene_table_P0_2 <-  read.csv("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking",sep = "\t")
gene_table_P4_1 <-  read.csv("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking", sep = "\t")
gene_table_P4_2 <-  read.csv("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking", sep = "\t")
gene_table_P7_1 <-  read.csv("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking", sep = "\t")
gene_table_P7_2 <-  read.csv("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking", sep = "\t")
gene_table_Ad_1 <- read.csv("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking", sep = "\t")
gene_table_Ad_2 <- read.csv("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking", sep = "\t")
diff_exp <- read.csv("/projectnb/bf528/users/van-gogh/project_2/programmer/cuffdiff_out/gene_exp.diff", sep = "\t", stringsAsFactors = F)


# 7.3

#concatenate the different phases and their FPKM values
combined_df <- merge(fpkm_matrix, gene_table_P0_1, by="tracking_id")
names(combined_df)[names(combined_df) == 'FPKM'] <- 'P0_1_FPKM'
ordered_df_q_value <- diff_exp[order(diff_exp$q_value),]
diff_exp_1k <- ordered_df_q_value[1:1000,]
diff_exp_significant <- diff_exp_1k[diff_exp_1k$significant == "yes", ]
df_final_ordered_diff_exp <- diff_exp_significant[1:200,]
fpkm_value_diff_exp = combined_df[combined_df$gene_short_name %in% df_final_ordered_diff_exp$gene,]
fpkm_diff_exp_final = fpkm_value_diff_exp[complete.cases(fpkm_value_diff_exp), ]
fpkm_diff_exp_final_nd = fpkm_diff_exp_final[!duplicated(fpkm_diff_exp_final$gene_short_name),]

#create the final dataframe 
fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[,c("Ad_1_FPKM","Ad_2_FPKM","P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "gene_short_name")]
rownames(fpkm_diff_exp_final_nd)<-fpkm_diff_exp_final_nd$gene_short_name
fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[,c("Ad_1_FPKM","Ad_2_FPKM","P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM")]
fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[!(apply(fpkm_diff_exp_final_nd, 1, function(y) any(y == 0))),]

# create the matrix out of the dataframe to provide it as an input to the heatmap
matrix_df <- as.matrix(fpkm_diff_exp_final_nd)

#plot the heatmap
pheatmap(matrix_df, scale = "row", fontsize_row = 4, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_row = T,
         cluster_col = T,
         show_rownames = T,
         show_colnames = T,
         cex = 0.8,
         legend = T,
         fontsize = 10)

#save the heatmap
ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/heatMap.png")


#7.2

#upregulated 

#read the up regulated csv
up_regulated <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/up_reg_gene_output.csv", header = F)

#read the reference paper csv
up_regulated_reference <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/up_reference_research_paper.csv", header = F)

#create a column for overlapping row called similar_row_values after comparing the two csv's
similar_row_values <- match(up_regulated[,2], up_regulated_reference[,2])
similar_row_values <- ifelse(is.na(similar_row_values), "No", "Yes")
go_terms_up_regulated <- cbind(up_regulated, similar_row_values)
table(go_terms_up_regulated$similar_row_values)

#write the csv for the upregulated genes
write.csv(go_terms_up_regulated,"up_regulated_GO_clusters.csv")

#downregulated

#read the up downregulated csv
down_regulated <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/down_reg_gene_output.csv", header = F)

#read the reference paper csv
down_regulated_reference <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/down_reference_research_paper.csv", header = F)

#create a column for overlapping row called similar_row_values after comparing the two csv's
similar_row_values <- match(down_regulated[,2], down_regulated_reference[,2])
similar_row_values <- ifelse(is.na(similar_row_values), "No", "Yes")
go_terms_down_regulated <- cbind(down_regulated, similar_row_values)
table(go_terms_down_regulated$similar_row_values)

#write the csv for the downregulated genes
write.csv(go_terms_down_regulated,"down_regulated_GO_clusters.csv")












