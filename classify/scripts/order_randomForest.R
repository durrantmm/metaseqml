library(mlr)
library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

features_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, recursive=TRUE)

features <- list.files(features_dir, full.names=T)
features <- features[grepl("R1", features)]

max_lines <- 1000000
n_max <-as.integer(max_lines / length(features))
training_set <- NULL

counter <- 1
for (f in features){
  print(paste("Reading in file ", counter, "of", length(features)))
  data_in <- read_tsv(f, n_max=n_max) %>% select(-Seq, -Kingdom, -Phylum, -Class, -Family, -Genus, -Species)
  
  if (is.null(data_in)){
    training_set <- data_in
  }
  
  training_set <- rbind(training_set, data_in)
  counter <- counter + 1
}

classif.lrn = makeLearner("classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)
classif.task = makeClassifTask(data = data.frame(train_dataset), target = "Order")

mod = train(classif.lrn, classif.task)

out_rds <- file.path(output_dir, 'order.randomForest.rds')
saveRDS(mod, out_rds)
