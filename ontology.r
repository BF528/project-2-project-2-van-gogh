up_c1 <- read.table("up_reg_gene_output_c1.txt", sep = '\t', header=TRUE)
up_c2 <- read.table("up_reg_gene_output_c2.txt", sep = '\t', header=TRUE)
up_c3 <- read.table("up_reg_gene_output_c3.txt", sep = '\t', header=TRUE)
dn_c1 <- read.table("down_reg_gene_output_c1.txt", sep = '\t', header=TRUE)
dn_c2 <- read.table("down_reg_gene_output_c2.txt", sep = '\t', header=TRUE)
dn_c3 <- read.table("down_reg_gene_output_c3.txt", sep = '\t', header=TRUE)
dn_c4 <- read.table("down_reg_gene_output_c4.txt", sep = '\t', header=TRUE)

#Genes that are either up or down regulated 
up_c1_gene <- unlist(strsplit(up_c1$Genes[1:5], ","))
up_c2_gene <- unlist(strsplit(up_c2$Genes[1:5], ",")) 
up_c3_gene <- unlist(strsplit(up_c3$Genes[1:5], ","))
dn_c1_gene <- unlist(strsplit(dn_c1$Genes[1:5], ","))
dn_c2_gene <- unlist(strsplit(dn_c2$Genes[1:3], ","))
dn_c3_gene <- unlist(strsplit(dn_c3$Genes[1:4], ","))
dn_c4_gene <- unlist(strsplit(dn_c4$Genes[1:4], ","))

#Extract all the 
up_c1_n <- sapply(strsplit(up_c1$Term, "~"), "[[", 2)[1:5]
up_c2_n <- sapply(strsplit(up_c2$Term, "~"), "[[", 2)[1:5]
up_c3_n <- sapply(strsplit(up_c3$Term, "~"), "[[", 2)[1:5]
dn_c1_n <- sapply(strsplit(dn_c1$Term, "~"), "[[", 2)[1:5]
dn_c2_n <- sapply(strsplit(dn_c2$Term, "~"), "[[", 2)[1:3]
dn_c3_n <- sapply(strsplit(dn_c3$Term, "~"), "[[", 2)[1:4]
dn_c4_n <- sapply(strsplit(dn_c4$Term, "~"), "[[", 2)[1:4]

# Top clustering ontology 
u1 <- data.frame(Term=up_c1_n, Count=up_c1$Count[1:5], Cluster="up_c1")
u2 <- data.frame(Term=up_c2_n, Count=up_c2$Count[1:5], Cluster="up_c2")
u3 <- data.frame(Term=up_c3_n, Count=up_c3$Count[1:5], Cluster="up_c3")
d1 <- data.frame(Term=dn_c1_n, Count=dn_c1$Count[1:5], Cluster="dn_c1")
d2 <- data.frame(Term=dn_c2_n, Count=dn_c2$Count[1:3], Cluster="dn_c2")
d3 <- data.frame(Term=dn_c3_n, Count=dn_c3$Count[1:4], Cluster="dn_c3")
d4 <- data.frame(Term=dn_c4_n, Count=dn_c4$Count[1:4], Cluster="dn_c4")
df_table <- rbind(u1, u2, u3, d1, d2, d3, d4)
write.csv(df_table, "ontology_table.csv")
