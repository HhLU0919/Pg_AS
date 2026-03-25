#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(edgeR))
suppressMessages(library(argparser))#https://github.com/cran/argparser

# Define command-line arguments
p <- arg_parser("Merge featureCounts(linux version) and calculate FPKM")

p <- add_argument(p, "--input_path", help="input: a directory containing the counts matrix named with '<sample>.count'", type="character",default = "./")
p <- add_argument(p, "--pattern", help="limit pattern of input files using regular expression in R language",type="character",default = "*count$")
p <- add_argument(p, "--output_path", help="output: an existent directory", type="character",default = "./")
p <- add_argument(p, "--prefix", help="give the file of output matrix a prefix like '<output_prefix>_genes.*'", type="character",default = "my",short = "-f")

# Parse command-line arguments
argv <- parse_args(p)

path <- argv$input_path
pattern <- argv$pattern
output_path <- argv$output_path
output_prefix <- argv$prefix

# Main workflow
file_name <- dir(path = path,pattern = pattern)
file <- paste0(path,"/",file_name)

# merge all count matrix
df <- read.table(file[1], header = T,comment.char = "#") %>% select(c(1,6,7))
colnames(df)[3] <- basename(file[1]) %>% str_remove("\\.\\w+$")
cat("1 count matrix has merged!\n")

for (i in 2:length(file)) {
  df_tmp <- read.table(file[i], header = T,comment.char = "#") %>% select(c(1,7))
  colnames(df_tmp)[2] <- basename(file[i]) %>% str_remove("\\.\\w+$")
  df <- df %>% full_join(df_tmp,by="Geneid")
  cat(i," count matrixs have merged!\n")
}

cat("Congratulations! All count matrixs have merged!\n")

count_mat <- df %>% select(-2) %>% tibble::column_to_rownames(var = "Geneid")
gene_length <- df %>% select(1,2) %>% tibble::column_to_rownames(var = "Geneid")

# calculate FPKM
cat("Calculating the FPKM…\n")
fpkm <- rpkm(count_mat,gene.length = gene_length$Length) %>% as.data.frame()

# write out 
cat("writing out raw counts matrix\n")
out_count_mat <- count_mat %>% rownames_to_column(var = "gene_id")
count_file <- paste0(output_path,output_prefix,"_","genes.counts")

write.table(out_count_mat,
            file = count_file,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat("writing out FPKM\n")

out_fpkm <- fpkm %>% rownames_to_column(var = "gene_id")
fpkm_file <- paste0(output_path,output_prefix,"_","genes.fpkm")

write.table(out_fpkm,
            file = fpkm_file,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat("Congratulations! All of the missions have been completed!\n")
