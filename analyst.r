setwd("/projectnb/bf528/users/van-gogh/project_2/analyst")

df <- read.table("/projectnb/bf528/users/van-gogh/project_2/programmer/cuffdiff_out/gene_exp.diff", header=TRUE)
df_ok <- subset(df, status == 'OK')
df_checked <- subset(df_ok, value_1 > 1 | value_2 > 1)
sort_q <- df_checked[order(df_checked$q_value), ]

top_10_table <- sort_q[1:10, c("gene", "value_1", "value_2", "log2.fold_change.", "p_value", "q_value")]

write.csv(top_10_table, "top_10_genes.csv")

#Plot histogram for log2.foldchange
#For all the genes
hist(sort_q$log2.fold_change., breaks = 30, main="Log2 Fold Change for All Genes", xlab = "Log2 Fold Change Value")
#For significant genes
df_sig <- subset(df, df$significant=='yes')
hist(df_sig$log2.fold_change., breaks = 30, main="Log2 Fold Change for Significant Genes", xlab = "Log2 Fold Change Value")

up_reg_genes <- subset(df_sig, df_sig$log2.fold_change.>0)
down_reg_genes <- subset(df_sig, df_sig$log2.fold_change.<0)
dim(up_reg_genes)
dim(down_reg_genes)

write.csv(up_reg_genes$gene, row.names=FALSE, quote =FALSE, file = "up_reg_genes.csv")
write.csv(down_reg_genes$gene, row.names=FALSE, quote =FALSE, file = "down_reg_genes.csv")
