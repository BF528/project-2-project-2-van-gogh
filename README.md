# Project Description
This project focuses on reproducing the results from O’Meara et al. The objectives of the study were to determine if myocytes revert the transcriptional phenotype to a less differentiated state during regeneration and to systematically interrogate the transcriptional data to identify and validate potential regulators of this process.

In reproducing the study we downloaded and analyzed P0 sequencing data, aligned and normalized the reads, evaluated quality control metrics, and interpreted differential expression values using DAVID. By doing so, we were able to better understand how neonatal mice are able to regenerate their heart tissue but lose this ability later in development.

O’Meara et al. Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circ Res. Feb 2015. PMID: 25477501

# Contributors
* Data Curator: Andrew Gjelsteen (@agjelste)
* Programmer: Elysha Sameth (@esameth)
* Analyst: Lindsay Wang (@LindsayW007)
* Biologist: Monil Gandhi (@gandhimonil9823)

# Repository Contents
## Data Curator
### run_extract.qsub
* Dependencies: sratoolkit
* Execution: `qsub run_extract.qsub`
* Outputs: two fastq files
* Extracts two fastq files from sra file
* sratoolkit parameters are `fastq-dump -I --split-files <sra_file.sra> -O <output_directory_path>`

### fastqc command
* Dependencies: java, fastqc
* Execution: `fastqc <file1.fastq> <file2.fastq> -o <output_directory_path>`
* Outputs: two html files
* Performs quality assessment on the two fastq input files to output two html files containing tables and images of quality metrics.

## Programmer
### run_tophat.qsub
* Dependencies: python2, samtools-0.1.19, bowtie2, boost, tophat
* Execution: `qsub run_tophat.qsub`
* Outputs: `P0_1_tophat/accepted_hits.bam`
* Aligns the reads to the mouse genome reference called mm9
* TopHat parameters are `-r 200`, `-G`, `--segment-length=20`, `--segment-mismatches=1`, `--no-novel-juncs`

### run_rseqc.qsub
* Dependencies: python3, samtools, rseqc
* Execution: `qsub run_rseqc.qsub`
* Outputs: quality control metric in `bam_stat.txt`, plots of mRNA distribution and gene coverage
* `genebody_coverage.py`: Calculate the RNA-seq reads coverage over gene body
* `inner_distance.py`: Calculate the inner distance (insert size)  of RNA-seq fragments
* `bam_stat.py`: Summarizing mapping statistics of a BAM or SAM file 

### run_cufflinks.qsub
* Dependencies: cufflinks
* Execution: `qsub run_cufflinks.qsub` 
* Inputs: `P0_1_tophat/accepted_hits.bam`
* Outputs: `P0_1_cufflinks/genes.fpkm_tracking`
* Counts how reads map to genomic regions defined by an annotation
* Cufflinks parameters are `--compatible-hits-norm` 

### run_cuffdiff.qsub
* Dependencies: cufflinks
* Execution: `qsub rub_cuffdiff.qsub`
* Inputs: `P0_1_tophat/accepted_hits.bam`
* Identifies differentially expressed genes between P0 and Ad

## Analyst
### analyst.r
* Dependencies: r
* Inputs: `gene_exp.diff`
* Outputs: `Log2 Fold Change for All Genes`: Histogram of log2 fold-change for all genes
* `Log2 Fold Change for Significant Geness`: Histogram of log2 fold-change for significant genes
* `up_reg_genes.csv`: List of genes that are up-regulated 
* `down_reg_genes.csv`: List of genes that are down-regulated 

### ontology.r
* Dependencies: r
* Inputs: GO analysis results
* Outputs: `ontology_table.csv`: Contain information about Top 3 enrichment clusters' top 5 GO terms, # of genes cluster for the term, symbol of the cluster for up-regulated genes. And Top 4 enrichment clusters' top 5 GO terms, # of genes cluster for the term, symbol of the cluster for for down regulated genes.

### GO_barplot.r
* Dependencies: r
* Outputs: `up.png`: Barplot summarizing the GO analysis for up-regulated genes
* `down.png`: Barplot summarizing the GO analysis for down-regulated genes
