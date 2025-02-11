library(pROC)

options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output_prefix>\n")
}

input_file <- args[1]
plot_title <- args[2]
output <- args[3]

output_pdf <- paste0(output, ".methylation.roc.pdf")
title_final <- paste0("Methylation.", plot_title)

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


ctr <- subset(dat, Group == "0")[, c("Sample", "Group", "Diff")]
BRCA <- subset(dat, Group == "1")[, c("Sample", "Group", "Diff")]
COREAD <- subset(dat, Group == "2")[, c("Sample", "Group", "Diff")]
ESCA <- subset(dat, Group == "3")[, c("Sample", "Group", "Diff")]
STAD <- subset(dat, Group == "4")[, c("Sample", "Group", "Diff")]
LIHC <- subset(dat, Group == "5")[, c("Sample", "Group", "Diff")]
NSCLC <- subset(dat, Group == "6")[, c("Sample", "Group", "Diff")]
PACA <- subset(dat, Group == "7")[, c("Sample", "Group", "Diff")]

brca <- rbind(ctr, BRCA)
coread <- rbind(ctr, COREAD)
esca <- rbind(ctr, ESCA)
stad <- rbind(ctr, STAD)
lihc <- rbind(ctr, LIHC)
nsclc <- rbind(ctr, NSCLC)
paca <- rbind(ctr, PACA)

roc.1 <- roc(brca, Group, Diff)
auc.1 <- sprintf("%.4f", roc.1$auc)
p.1 <- sprintf("%.4e", pROC.p(roc.1))

roc.2 <- roc(coread, Group, Diff)
auc.2 <- sprintf("%.4f", roc.2$auc)
p.2 <- sprintf("%.4e", pROC.p(roc.2))

roc.3 <- roc(esca, Group, Diff)
auc.3 <- sprintf("%.4f", roc.3$auc)
p.3 <- sprintf("%.4e", pROC.p(roc.3))

roc.4 <- roc(stad, Group, Diff)
auc.4 <- sprintf("%.4f", roc.4$auc)
p.4 <- sprintf("%.4e", pROC.p(roc.4))

roc.5 <- roc(lihc, Group, Diff)
auc.5 <- sprintf("%.4f", roc.5$auc)
p.5 <- sprintf("%.4e", pROC.p(roc.5))

roc.6 <- roc(nsclc, Group, Diff)
auc.6 <- sprintf("%.4f", roc.6$auc)
p.6 <- sprintf("%.4e", pROC.p(roc.6))

roc.7 <- roc(paca, Group, Diff)
auc.7 <- sprintf("%.4f", roc.7$auc)
p.7 <- sprintf("%.4e", pROC.p(roc.7))

plot(roc.1, lwd = 3, col = "darkorange", cex.lab = 1.2, main = title_final)
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

