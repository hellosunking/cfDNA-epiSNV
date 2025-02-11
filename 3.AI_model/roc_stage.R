library( pROC );
library(binom);

options( stringsAsFactors=F );

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output>\n")
}

input_file <- args[1]
plot_title <- args[2]
output_pdf <- args[3]

pdf( output_pdf, width = 8, height = 8);

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

dat = read.table( input_file, head=T );

dat$Stage <- ifelse(dat$Stage == "Control", 0, 
                    ifelse(dat$Stage == "StageI", 1, 
                           ifelse(dat$Stage == "StageII", 2, 
                                  ifelse(dat$Stage == "StageIII", 3, 
                                         ifelse(dat$Stage == "StageIV", 4, 
                                                ifelse(dat$Stage == "StageNA", 5, dat$Stage))))))

ctr = subset(dat, Stage=="0")[, c("Sample", "Stage", "pred")]
s1 = subset(dat, Stage=="1")[, c("Sample", "Stage", "pred")]
s2 = subset(dat, Stage=="2")[, c("Sample", "Stage", "pred")]
s3 = subset(dat, Stage=="3")[, c("Sample", "Stage", "pred")]
s4 = subset(dat, Stage=="4")[, c("Sample", "Stage", "pred")]
s5 = subset(dat, Stage=="5")[, c("Sample", "Stage", "pred")]

S1 = rbind(ctr,s1);
S2 = rbind(ctr,s2);
S3 = rbind(ctr,s3);
S4 = rbind(ctr,s4);
S5 = rbind(ctr,s5);

## use pROC, which is BETTER than ROCR
roc.1 = roc( S1, Stage, pred );
auc.1 = sprintf("%.4f", roc.1$auc );
p.1   = sprintf("%.4e", pROC.p(roc.1));
roc.2 = roc( S2, Stage, pred );
auc.2 = sprintf("%.4f", roc.2$auc );
p.2   = sprintf("%.4e", pROC.p(roc.2));
roc.3 = roc( S3, Stage, pred );
auc.3 = sprintf("%.4f", roc.3$auc );
p.3   = sprintf("%.4e", pROC.p(roc.3));
roc.4 = roc( S4, Stage, pred );
auc.4 = sprintf("%.4f", roc.4$auc );
p.4   = sprintf("%.4e", pROC.p(roc.4));
roc.5 = roc( S5, Stage, pred );
auc.5 = sprintf("%.4f", roc.5$auc );
p.5   = sprintf("%.4e", pROC.p(roc.5));

plot( roc.1, lwd=3, col="deepskyblue", cex.lab=1.2, main=plot_title);
plot( roc.2, lwd=3, col="darkorange4", add=T );
plot( roc.3, lwd=3, col="darkorange", add=T );
plot( roc.4, lwd=3, col="#CD2626", add=T );
plot( roc.5, lwd=3, col="turquoise4", add=T );

legend( 'bottomright', legend=paste0(c("StageI", "StageII", "StageIII", "StageIV", "StageNA"), ", AUC=", c(auc.1, auc.2, auc.3, auc.4, auc.5), ", P=", c(p.1, p.2, p.3, p.4, p.5)),
	            lty=1, bty='n', cex=1.1, col=c("deepskyblue", "darkorange4", "darkorange", "#CD2626", "turquoise4") );

results <- list()

for (i in 1:5) {
  roc_obj <- get(paste0("roc.", i))
  
  coords <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity", "precision"))
  target_coord <- coords[coords[, "specificity"] >= 0.95, ][which.max(coords[coords[, "specificity"] >= 0.95, "sensitivity"]), ]

  threshold_raw <- target_coord["threshold"]
  sensitivity_raw <- target_coord["sensitivity"]
  specificity_raw <- target_coord["specificity"]

  threshold <- threshold_raw$threshold
  sensitivity <- sensitivity_raw$sensitivity
  specificity <- specificity_raw$specificity
  
  num_positive_detected <- sum(dat$pred >= threshold & dat$Stage == i)

  ci <- binom.confint(num_positive_detected, length(roc_obj$case), conf.level = 0.95, methods = "exact")

  results[[i]] <- list(sensitivity = sensitivity, specificity = specificity, 
                       ci_lower = ci$lower, ci_upper = ci$upper, 
                       num_positive_detected = num_positive_detected)
  
  cat(paste0("Stage ", i, ": Sensitivity at 95% specificity = ", sensitivity, 
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

  num_positive_detected <- sum(dat$pred >= threshold & dat$Stage == i)

  ci <- binom.confint(num_positive_detected, length(roc_obj$case), conf.level = 0.95, methods = "exact")

  results[[i]] <- list(sensitivity = sensitivity, specificity = specificity,
                       ci_lower = ci$lower, ci_upper = ci$upper,
                       num_positive_detected = num_positive_detected)
  
  cat(paste0("Stage ", i, ": Sensitivity at 99% specificity = ", sensitivity,
             ", Threshold = ", threshold,
             ", 95% CI = [", ci$lower, ", ", ci$upper,
             "], Positive cases detected = ", num_positive_detected, "\n"))
}

