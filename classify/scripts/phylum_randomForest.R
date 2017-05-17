library(mlr)
library(readr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

features_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, recursive=TRUE)

features <- list.files(features_dir, full.names=T)
features <- features[grepl("R1", features)]

max_lines <- 5000
n_max <-as.integer(max_lines / length(features))
training_set <- NULL
for (f in features){
  data_in <- read_tsv(f, n_max=n_max) %>% select(-Seq, -Kingdom, -Class, -Order, -Family, -Genus, -Species)

  if (is.null(data_in)){
    training_set <- data_in
  }
  
  training_set <- rbind(training_set, train_data)
}

classif.lrn = makeLearner("classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)
classif.task = makeClassifTask(data = data.frame(train_dataset), target = "Phylum")

mod = train(classif.lrn, classif.task)

out_rds <- file.path(output_dir, 'phylum.randomForest.rds')
saveRDS(mod, out_rds)
