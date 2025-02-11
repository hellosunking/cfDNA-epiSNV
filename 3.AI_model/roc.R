library(pROC)
library(binom)

options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output>\n")
}

input_file <- args[1]
plot_title <- args[2]
output_pdf <- args[3]

pdf(output_pdf, width = 8, height = 8)

pROC.p <- function(rr) {
  v <- var(rr)
  b <- rr$auc - 0.5
  se <- sqrt(v)
  z <- b / se
  p <- 2 * pt(-abs(z), df = Inf)
  p
}

dat <- read.table(input_file, head = T)

dat$Group <- ifelse(dat$Group == "Control", 0, 
                    ifelse(dat$Group == "BRCA", 1, 
                           ifelse(dat$Group == "COREAD", 2, 
                                  ifelse(dat$Group == "ESCA", 3, 
                                         ifelse(dat$Group == "STAD", 4, 
                                                ifelse(dat$Group == "LIHC", 5, 
                                                       ifelse(dat$Group == "NSCLC", 6, 
                                                              ifelse(dat$Group == "PACA", 7, dat$Group))))))))


ctr <- subset(dat, Group == "0")[, c("Sample", "Group", "pred")]
BRCA <- subset(dat, Group == "1")[, c("Sample", "Group", "pred")]
COREAD <- subset(dat, Group == "2")[, c("Sample", "Group", "pred")]
ESCA <- subset(dat, Group == "3")[, c("Sample", "Group", "pred")]
STAD <- subset(dat, Group == "4")[, c("Sample", "Group", "pred")]
LIHC <- subset(dat, Group == "5")[, c("Sample", "Group", "pred")]
NSCLC <- subset(dat, Group == "6")[, c("Sample", "Group", "pred")]
PACA <- subset(dat, Group == "7")[, c("Sample", "Group", "pred")]

brca <- rbind(ctr, BRCA)
coread <- rbind(ctr, COREAD)
esca <- rbind(ctr, ESCA)
stad <- rbind(ctr, STAD)
lihc <- rbind(ctr, LIHC)
nsclc <- rbind(ctr, NSCLC)
paca <- rbind(ctr, PACA)

roc.1 <- roc(brca, Group, pred)
auc.1 <- sprintf("%.4f", roc.1$auc)
p.1 <- sprintf("%.4e", pROC.p(roc.1))

roc.2 <- roc(coread, Group, pred)
auc.2 <- sprintf("%.4f", roc.2$auc)
p.2 <- sprintf("%.4e", pROC.p(roc.2))

roc.3 <- roc(esca, Group, pred)
auc.3 <- sprintf("%.4f", roc.3$auc)
p.3 <- sprintf("%.4e", pROC.p(roc.3))

roc.4 <- roc(stad, Group, pred)
auc.4 <- sprintf("%.4f", roc.4$auc)
p.4 <- sprintf("%.4e", pROC.p(roc.4))

roc.5 <- roc(lihc, Group, pred)
auc.5 <- sprintf("%.4f", roc.5$auc)
p.5 <- sprintf("%.4e", pROC.p(roc.5))

roc.6 <- roc(nsclc, Group, pred)
auc.6 <- sprintf("%.4f", roc.6$auc)
p.6 <- sprintf("%.4e", pROC.p(roc.6))

roc.7 <- roc(paca, Group, pred)
auc.7 <- sprintf("%.4f", roc.7$auc)
p.7 <- sprintf("%.4e", pROC.p(roc.7))

plot(roc.1, lwd = 3, col = "darkorange", cex.lab = 1.2, main = plot_title)
plot(roc.2, lwd = 3, col = "blue", add = T)
plot(roc.3, lwd = 3, col = "deepskyblue", add = T)
plot(roc.4, lwd = 3, col = "turquoise4", add = T)
plot(roc.5, lwd = 3, col = "#CD2626", add = T)
plot(roc.6, lwd = 3, col = "purple4", add = T)
plot(roc.7, lwd = 3, col = "darkorange4", add = T)

legend('bottomright', legend = paste0(c("BRCA", "COREAD", "ESCA", "STAD", "LIHC", "NSCLC", "PACA"), 
                                      ", AUC=", c(auc.1, auc.2, auc.3, auc.4, auc.5, auc.6, auc.7), 
                                      ", P=", c(p.1, p.2, p.3, p.4, p.5, p.6, p.7)),
       lty = 1, bty = 'n', cex = 1.1, col = c("darkorange", "blue", "deepskyblue", "turquoise4", "#CD2626", "purple4", "darkorange4"))

results <- list()

for (i in 1:7) {
  roc_obj <- get(paste0("roc.", i))
  
  coords <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity", "precision"))
  target_coord <- coords[coords[, "specificity"] >= 0.95, ][which.max(coords[coords[, "specificity"] >= 0.95, "sensitivity"]), ]

  threshold_raw <- target_coord["threshold"]
  sensitivity_raw <- target_coord["sensitivity"]
  specificity_raw <- target_coord["specificity"]

  threshold <- threshold_raw$threshold
  sensitivity <- sensitivity_raw$sensitivity
  specificity <- specificity_raw$specificity
  
  num_positive_detected <- sum(dat$pred >= threshold & dat$Group == i)

  ci <- binom.confint(num_positive_detected, length(roc_obj$case), conf.level = 0.95, methods = "exact")

  results[[i]] <- list(sensitivity = sensitivity, specificity = specificity, 
                       ci_lower = ci$lower, ci_upper = ci$upper, 
                       num_positive_detected = num_positive_detected)
  
  cat(paste0("Group ", i, ": Sensitivity at 95% specificity = ", sensitivity, 
             ", Threshold = ", threshold,
             ", 95% CI = [", ci$lower, ", ", ci$upper, 
             "], Positive cases detected = ", num_positive_detected, "\n"))

  roc_obj <- get(paste0("roc.", i))

  coords <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity", "precision"))
  target_coord <- coords[coords[, "specificity"] >= 0.99, ][which.max(coords[coords[, "specificity"] >= 0.99, "sensitivity"]), ]

  threshold_raw <- target_coord["threshold"]
  sensitivity_raw <- target_coord["sensitivity"]
  specificity_raw <- target_coord["specificity"]

  threshold <- threshold_raw$threshold
  sensitivity <- sensitivity_raw$sensitivity
  specificity <- specificity_raw$specificity

  num_positive_detected <- sum(dat$pred >= threshold & dat$Group == i)

  ci <- binom.confint(num_positive_detected, length(roc_obj$case), conf.level = 0.95, methods = "exact")

  results[[i]] <- list(sensitivity = sensitivity, specificity = specificity,
                       ci_lower = ci$lower, ci_upper = ci$upper,
                       num_positive_detected = num_positive_detected)

  cat(paste0("Group ", i, ": Sensitivity at 99% specificity = ", sensitivity,
             ", Threshold = ", threshold,
             ", 95% CI = [", ci$lower, ", ", ci$upper,
             "], Positive cases detected = ", num_positive_detected, "\n"))
}


