# Monil Gandhi - Biologist

#File Description - Created line plots for the different phases of Sacromere, Mitochondria and Cell Cycle


# install packages and load libraries
install.packages("pacman")
require(pacman)  
library(pacman)  
library(reshape2)
library(grid)
library(gridExtra)

#run pacman for downloading the necessary libraries
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
              ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
             stringr, tidyr, ggpubr, grid, gridExtra) 

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

#7.1

#list of genes for each Sacromere, Mitochondria and Cell Cylce
list_sarcomere_gene <- list("Pdlim5","Pygm","Myoz2","Des","Csrp3", "Tcap","Cryab")
list_mitochondria_gene <- list("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
list_cell_cyle_gene <- list("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27", "E2f1","Bora","Cdc45","Rad51","Aurkb","Cdc23")

#function to return the fpkm values of the list of genes for Sacromere, Mitochondria
fpkm_values <- function(list_gene_names, gene_phase_table) {
  list_fpkm <- c()
  i <- 1
  
  for (gene in list_gene_names){
    val_gene <- subset(gene_phase_table, gene_short_name == gene)

    val_fpkm <-  val_gene$FPKM

    list_fpkm[i] <- val_fpkm

    i <- i + 1

   
  }
  return(list_fpkm)
}

#function to return the fpkm values of the list of genes for cell cycle
fpkm_values_cycle_cycle <- function(list_gene_names, gene_phase_table) {
  list_fpkm <- c()
  i <- 1
  
  for (gene in list_gene_names){
    val_gene <- subset(gene_phase_table, gene_short_name == gene)
    
    val_fpkm <-  val_gene$FPKM
    
    
    if (!identical(val_fpkm, numeric(0))){
      list_fpkm[i] = val_fpkm}
    else{
      list_fpkm[i] = 0}
    
    i <- i + 1
  
  }
  return(list_fpkm)
}

#function to plot the line graph 
plot_function <- function(list_gene, df_to_return){
  
  df_to_return <- as.data.frame(t(df_to_return))
  
  colnames(df_to_return) <- list_gene
  
  df_to_return <- cbind(Phase = rownames(df_to_return), df_to_return)
  rownames(df_to_return) <- 1:nrow(df_to_return)
  
  df_to_return = melt(df_to_return)

  colnames(df_to_return)[2] <- "Genes"
  
  return(df_to_return)
  
}


#sarcomere_1

#get the FPKM vlaues for the sarcomere
S_P0_1 <- fpkm_values(list_sarcomere_gene, gene_table_P0_1)
S_P4_1 <- fpkm_values(list_sarcomere_gene, gene_table_P4_1)
S_P7_1 <- fpkm_values(list_sarcomere_gene, gene_table_P7_1)
S_Ad_1 <- fpkm_values(list_sarcomere_gene, gene_table_Ad_1)

#create the dataframe of different phases
df_sarcomere<- data.frame(S_P0_1, S_P4_1,  
                          S_P7_1,  S_Ad_1)

#plot the line plot
df_to_plot <- plot_function(list_sarcomere_gene,df_sarcomere)
level_order <- c("S_P0_1", "S_P4_1","S_P7_1",  "S_Ad_1")

#save the line plot
pdf(file = "/projectnb/bf528/users/van-gogh/project_2/biologist/Sarcomere_Sample_1.png",  
    width = 4, # The width of the plot in inches
    height = 4)

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

sarcomere_1 <- p + ggtitle("Sarcomere_Sample_1") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))

ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Sarcomere_Sample_1.png")


#sarcomere_2

#get the FPKM values 
S_P0_2 = fpkm_values(list_sarcomere_gene, gene_table_P0_2)
S_P4_2 = fpkm_values(list_sarcomere_gene, gene_table_P4_2)
S_P7_2 = fpkm_values(list_sarcomere_gene, gene_table_P7_2)
S_Ad_2 = fpkm_values(list_sarcomere_gene, gene_table_Ad_2)

#create the dataframe of different phases
df_sarcomere<- data.frame(S_P0_2, S_P4_2,  
                          S_P7_2,  S_Ad_2)

df_to_plot <- plot_function(list_sarcomere_gene,df_sarcomere)

level_order <- c("S_P0_2", "S_P4_2","S_P7_2",  "S_Ad_2")

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

sarcomere_2 <- p + ggtitle("Sarcomere_Sample_2") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))

ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Sarcomere_Sample_2.png")


#mitochondria_1

#get the FPKM values 
M_P0_1 = fpkm_values(list_mitochondria_gene, gene_table_P0_1)
M_P0_1[is.na(M_P0_1)] <- 0

M_P4_1 = fpkm_values(list_mitochondria_gene, gene_table_P4_1)
M_P4_1[is.na(M_P4_1)] <- 0

M_P7_1 = fpkm_values(list_mitochondria_gene, gene_table_P7_1)
M_P7_1[is.na(M_P7_1)] <- 0

M_Ad_1 = fpkm_values(list_mitochondria_gene, gene_table_Ad_1)
M_Ad_1[is.na(M_Ad_1)] <- 0

#create the dataframe of different phases
df_mitochondria<- data.frame(M_P0_1, M_P4_1,  
                          M_P7_1,  M_Ad_1)

df_to_plot <- plot_function(list_mitochondria_gene,df_mitochondria)

level_order <- c("M_P0_1", "M_P4_1","M_P7_1",  "M_Ad_1")

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

mict_1 <- p + ggtitle("Mitochondria_Sample_1") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))

ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Mitochondria_Sample_1.png")


#mitochondria_2

#get the FPKM values 
M_P0_2 = fpkm_values(list_mitochondria_gene, gene_table_P0_2)
M_P0_2[is.na(M_P0_2)] <- 0

M_P4_2 = fpkm_values(list_mitochondria_gene, gene_table_P4_2)
M_P4_2[is.na(M_P4_2)] <- 0

M_P7_2 = fpkm_values(list_mitochondria_gene, gene_table_P7_2)
M_P7_2[is.na(M_P7_2)] <- 0

M_Ad_2 = fpkm_values(list_mitochondria_gene, gene_table_Ad_2)
M_Ad_2[is.na(M_Ad_2)] <- 0

#create the dataframe of different phases
df_mitochondria<- data.frame(M_P0_2, M_P4_2,  
                             M_P7_2,  M_Ad_2)

df_to_plot <- plot_function(list_mitochondria_gene,df_mitochondria)

level_order <- c("M_P0_2", "M_P4_2","M_P7_2",  "M_Ad_2")

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

mict_2 <- p + ggtitle("Mitochondria_Sample_2") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))


ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Mitochondria_Sample_2.png")

#cycle_cycle_1

#get the FPKM values 
C_P0_1 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P0_1)
C_P4_1 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P4_1)
C_P7_1 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P7_1)
C_Ad_1 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_Ad_1)

#create the dataframe of different phases
df_cycle_cycle<- data.frame(C_P0_1, C_P4_1,  
                             C_P7_1,  C_Ad_1)

df_to_plot <- plot_function(list_cell_cyle_gene,df_cycle_cycle)

level_order <- c("C_P0_1", "C_P4_1","C_P7_1",  "C_Ad_1")

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

cell_cycle_1 <- p + ggtitle("Cell_Cycle_Sample_1") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))

ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Cell_Cycle_Sample_1.png")

#cycle_cycle_2

#get the FPKM values 
C_P0_2 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P0_2)
C_P4_2 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P4_2)
C_P7_2 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_P7_2)
C_Ad_2 = fpkm_values_cycle_cycle(list_cell_cyle_gene, gene_table_Ad_2)

#create the dataframe of different phases
df_cycle_cycle<- data.frame(C_P0_2, C_P4_2,  
                            C_P7_2,  C_Ad_2)

df_to_plot <- plot_function(list_cell_cyle_gene,df_cycle_cycle)

level_order <- c("C_P0_2", "C_P4_2","C_P7_2",  "C_Ad_2")

# plot the line plot
p <- ggplot(data=df_to_plot, aes(x=factor(Phase, level = level_order), y=value , group=Genes)) +
  geom_line(aes(color=Genes))+
  geom_point()

cell_cycle_2 <- p + ggtitle("Cell_Cycle_Sample_2") +
  xlab("Phases") + ylab("FPKM") + theme(plot.title = element_text(hjust = 0.5))

ggsave("/projectnb/bf528/users/van-gogh/project_2/biologist/Cell_Cycle_Sample_2.png")






