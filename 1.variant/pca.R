library(ggplot2)
library(ggrepel)

### args
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_file> <sample_info> <output_prefix>\n")
}

input.file <- args[1]
class.file <- args[2]
prefix     <- args[3]

### read
rt1   <- read.table(file=input.file, header=T,check.names=F)
rt1   <-  as.matrix( rt1[,-1] )
rt1   <- scale(rt1)
rt2   <- read.table(file=class.file)
class <- as.character( rt2[,2] )

### stat
pca <- prcomp(rt1)
s   <- summary(pca)
write.table(pca$x, paste(prefix,".x.txt",sep=""), quote = F, sep="\t")
write.table(s$importance, paste(prefix,".eig.txt",sep=""), quote = F, sep="\t")


x  <- pca$x[,1]
y  <- pca$x[,2]

pot.mds <- cbind(x,y,rt2)
colnames(pot.mds) <- c("PC1","PC2","sample","group")
PCC1 =round(s$importance[2,1]*100,d=2)
PCC2 =round(s$importance[2,2]*100,d=2)

pdf(file= paste(prefix, ".PCA.pdf", sep="") )
ggplot(data=pot.mds, aes(PC1, PC2)) +
  geom_point(aes(color=group), size=3) +
  scale_color_manual(values=c("HCC" = "#CD2626", "Control" = "grey50")) +
  scale_fill_manual(values=c("HCC" = "#CD2626", "Control" = "grey50")) +
  labs(title="PCA plot", x=paste("PC1", PCC1, " %"), y=paste("PC2", PCC2, " %")) +
  theme_bw() +
  theme(axis.title = element_text(family = "serif", face = "bold", size = 18, colour = "black")) +
  theme(axis.text = element_text(family = "serif", face = "bold", size = 16, color="black")) +
  theme(panel.grid=element_blank()) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0)

dev.off()
