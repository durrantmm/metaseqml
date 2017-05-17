if (!("mlr" %in% installed.packages()[,"Package"])){
  install.packages("mlr")  
  install.packages("randomForest")  
}

library(tidyr)
args = commandArgs(trailingOnly=TRUE)

output_file <- args[1]
write_tsv(data.frame(Done=c()), output_file)
