if (!("mlr" %in% installed.packages()[,"Package"])){
  install.packages("mlr")  
}
args = commandArgs(trailingOnly=TRUE)

seqtab <- readRDS(args[1], refhook = NULL)
silva_train_set <- args[2]
silva_species_train_set <- args[3]
silva_taxa_out_rds <- args[4]