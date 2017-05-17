if (!("mlr" %in% installed.packages()[,"Package"])){
  install.packages("mlr", repos=c("https://cran.cnr.berkeley.edu/"))  
  install.packages("randomForest", repos=c("https://cran.cnr.berkeley.edu/"))  
}

library(tidyr)
args = commandArgs(trailingOnly=TRUE)

output_file <- args[1]
write_tsv(data.frame(Done=c()), output_file)

