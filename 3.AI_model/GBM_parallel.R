library(pROC)
library(caret)
library(foreach)
library(doParallel)

args <- commandArgs(trailing = T)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input> <output_prefix>\n")
}
input <- args[1]
prefix <- args[2]

data <- read.table(file = input, header = TRUE, row.names = 1)

set.seed(7)
data$Type <- factor(data$Type, levels = c(0, 1), labels = c("class0", "class1"))

is_prime <- function(n) {
  if (n <= 6) return(FALSE)
  if (n == 7 || n == 11 || n == 13 || n == 17 || n == 19 || n == 23) return(TRUE)
  if (n %% 2 == 0 || n %% 3 == 0) return(FALSE)
  if (n > 23 && n < 29) return(FALSE)
  for (i in seq(5, floor(sqrt(n)), by = 6)) {
    if (n %% i == 0 || n %% (i + 2) == 0) return(FALSE)
  }
  return(TRUE)
}

generate_primes <- function(n, min_value) {
  primes <- c()
  num <- min_value
  while (length(primes) < n) {
    if (is_prime(num)) {
      primes <- c(primes, num)
    }
    num <- num + 1
  }
  return(primes)
}

prime_seeds <- generate_primes(101, 7)

fixed_seeds <- vector("list", length = 101)

for (i in seq_len(100)) {
  fixed_seeds[[i]] <- rep(prime_seeds[i], 480)
}
fixed_seeds[[101]] <- prime_seeds[101]

ctrl <- trainControl(method = "repeatedcv", 
                     number = 10, 
                     repeats = 10,
                     savePredictions = TRUE, 
                     classProbs = TRUE, 
                     search = "grid",
                     #seeds = fixed_seeds,
                     summaryFunction = twoClassSummary, 
                     verboseIter = TRUE)

param_grid <- expand.grid(
  n.trees = c(50,70,90,100,150,200), 
  interaction.depth = c(1,5,7,9), 
  shrinkage = c(0.01,0.05,0.1,0.2), 
  n.minobsinnode = c(5,10,15,20,30) 
)

prime_seeds_new <- generate_primes(100, 7)

# backend
#cl <- makeCluster(detectCores() - 1) #reserve one core
#registerDoParallel(cl)
registerDoParallel(cores = detectCores() - 1)

for (rep in 1:100) { 
  set.seed(prime_seeds_new[rep])
  folds <- createFolds(data$Type, k = 10) 
  
  # foreach 10 fold
  foreach(fold = 1:10, .packages = c("caret", "pROC"), .verbose = TRUE) %dopar% {
    train_indices <- unlist(folds[-fold])
    test_indices <- unlist(folds[fold]) 
    
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    set.seed(7)
    gbm_grid <- train(Type ~ ., 
                      data = train, 
                      method = "gbm", 
                      trControl = ctrl, 
                      tuneGrid = param_grid, 
                      metric = "ROC")

    best_param <- gbm_grid$bestTune
    
    param_grid_best <- expand.grid(
      n.trees = best_param$n.trees,
      interaction.depth = best_param$interaction.depth,
      shrinkage = best_param$shrinkage,
      n.minobsinnode = best_param$n.minobsinnode
    )
    set.seed(7)
    model <- train(Type ~ ., 
                   data = train, 
                   method = "gbm", 
                   trControl = trainControl(method = "none",
                                            savePredictions = TRUE, 
                                            classProbs = TRUE, 
                                            summaryFunction = twoClassSummary, 
                                            verboseIter = TRUE), 
                   tuneGrid = param_grid_best, 
                   metric = "ROC")
    
    name=paste0(prefix,".rep_",rep,"_fold_",fold,".rds")
    saveRDS(model, name)

    train_pred <- predict(model, newdata = train, type = "prob")
    train_Type <- train$Type
    train_pred_t <- data.frame(sid = rownames(train),
                               Type = train_Type,
                               pred = train_pred$class1)

    test_pred <- predict(model, newdata = test, type = "prob")
    test_Type <- test$Type
    test_pred_t <- data.frame(sid = rownames(test),
                              Type = test_Type,
                              pred = test_pred$class1)

    write.table(train_pred_t, file = paste0(prefix,".train.rep_", rep, "_fold_", fold, ".txt"), sep = "\t", col.names = F, row.names = FALSE, quote = FALSE)
    write.table(test_pred_t, file = paste0(prefix,".test.rep_", rep, "_fold_", fold, ".txt"), sep = "\t", col.names = F, row.names = FALSE, quote = FALSE)

    # output best param
    best_params_file <- paste0(prefix,".best_params_rep_", rep, "_fold_", fold, ".txt")
    write.table(as.data.frame(t(best_param)), file = best_params_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

# stop backend
#stopCluster(cl)
stopImplicitCluster()
