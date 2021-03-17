install.packages("pacman")

# Then load the package by using either of the following:
require(pacman)  # Gives a confirmation message.
library(pacman)  # No message.
library(reshape2)
library(grid)
library(gridExtra)
library(magrittr)

library(tibble)




# Or, by using "pacman::p_load" you can use the p_load
# function from pacman without actually loading pacman.
# These are packages I load every time.
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggpubr, grid, gridExtra, pheatmap, gdata, magrittr) 


#tracking_id.Ad_1_FPKM.Ad_2_FPKM.P0_2_FPKM.P4_1_FPKM.P4_2_FPKM.P7_1_FPKM.P7_2_FPKM
fpkm_matrix <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", sep = "\t")

#tracking_id.class_code.nearest_ref_id.gene_id.gene_short_name.tss_id.locus.length.coverage.FPKM.FPKM_conf_lo.FPKM_conf_hi.FPKM_status
gene_table_P0_1 <-  read.csv("/projectnb/bf528/users/van-gogh/project_2/programmer/P0_1_cufflinks/genes.fpkm_tracking", sep = "\t")
gene_table_P0_2 <-  read.csv("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking",sep = "\t")
gene_table_P4_1 <-  read.csv("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking", sep = "\t")
gene_table_P4_2 <-  read.csv("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking", sep = "\t")
gene_table_P7_1 <-  read.csv("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking", sep = "\t")
gene_table_P7_2 <-  read.csv("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking", sep = "\t")

gene_table_Ad_1 <- read.csv("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking", sep = "\t")
gene_table_Ad_2 <- read.csv("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking", sep = "\t")

diff_exp <- read.csv("/projectnb/bf528/users/van-gogh/project_2/programmer/cuffdiff_out/gene_exp.diff", sep = "\t", stringsAsFactors = F)


combined_df <- merge(fpkm_matrix, gene_table_P0_1, by="tracking_id")

names(combined_df)[names(combined_df) == 'FPKM'] <- 'P0_1_FPKM'

ordered_df_q_value <- diff_exp[order(diff_exp$q_value),]

diff_exp_1k <- ordered_df_q_value[1:1000,]

diff_exp_significant <- diff_exp_1k[diff_exp_1k$significant == "yes", ]


df_final_ordered_diff_exp <- diff_exp_significant[1:200,]

fpkm_value_diff_exp = combined_df[combined_df$gene_short_name %in% df_final_ordered_diff_exp$gene,]


fpkm_diff_exp_final = fpkm_value_diff_exp[complete.cases(fpkm_value_diff_exp), ]

fpkm_diff_exp_final_nd = fpkm_diff_exp_final[!duplicated(fpkm_diff_exp_final$gene_short_name),]


fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[,c("Ad_1_FPKM","Ad_2_FPKM","P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "gene_short_name")]

rownames(fpkm_diff_exp_final_nd)<-fpkm_diff_exp_final_nd$gene_short_name

fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[,c("Ad_1_FPKM","Ad_2_FPKM","P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM")]

fpkm_diff_exp_final_nd

fpkm_diff_exp_final_nd <- fpkm_diff_exp_final_nd[!(apply(fpkm_diff_exp_final_nd, 1, function(y) any(y == 0))),]

matrix_df <- as.matrix(fpkm_diff_exp_final_nd)

matrix_df

set_rownames(matrix_df("Genes"))
set_colnames(matrix_df("Samples"))



pheatmap(matrix_df, scale = "row", fontsize_row = 4, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         cluster_row = T,
         cluster_col = T,
         show_rownames = T,
         show_colnames = T,
         cex = 0.8,
         legend = T,
         fontsize = 10)

grid.text("Genes", y=-0.07, gp=gpar(fontsize=16))
grid.text("Samples", x=-0.07, rot=90, gp=gpar(fontsize=16))

# ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/heatMap.png")


# dev.off()






up_regulated <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/up_reg_gene_output.csv", header = F)

up_regulated_reference <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/up_reference_research_paper.csv", header = F)

similar_row_values <- match(up_regulated[,2], up_regulated_reference[,2])

similar_row_values <- ifelse(is.na(similar_row_values), "No", "Yes")

go_terms_up_regulated <- cbind(up_regulated, similar_row_values)

# go_terms_up_regulated_clusters <- table_upgo[1:90,]

# write.csv(go_terms_up_regulated,"up_regulated_GO_clusters.csv")




down_regulated <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/down_reg_gene_output.csv", header = F)

down_regulated_reference <- read.csv("/projectnb/bf528/users/van-gogh/project_2/biologist/down_reference_research_paper.csv", header = F)

similar_row_values <- match(down_regulated[,2], down_regulated_reference[,2])

similar_row_values <- ifelse(is.na(similar_row_values), "No", "Yes")

go_terms_down_regulated <- cbind(down_regulated, similar_row_values)

# go_terms_up_regulated_clusters <- table_upgo[1:90,]

# write.csv(go_terms_down_regulated,"down_regulated_GO_clusters.csv")
















