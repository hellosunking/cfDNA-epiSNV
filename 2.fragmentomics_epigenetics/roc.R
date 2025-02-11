library( pROC );

options( stringsAsFactors=F );
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output>\n")
}

input_file <- args[1]
plot_title <- args[2]
outfileName <- args[3]

pdf( outfileName );

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

dat = read.table( input_file, head=T );

dat$Group <- ifelse(dat$Group == "HCC", 1, ifelse(dat$Group == "Control", 0, dat$Group))
dat$Group <- factor(dat$Group, levels = c(0, 1))

ctr = subset(dat, Group=="0");
hcc = subset(dat, Group=="1");

roc.1 = roc( dat, Group, Diff );
auc.1 = sprintf("%.4f", roc.1$auc );
p.1   = sprintf("%.4e", pROC.p(roc.1));

plot( roc.1, lwd=2, col="black", cex.lab=1.2 );
legend( 'bottomright', legend=paste0(c("Diff_size"), ", AUC=", c(auc.1), ", P=", c(p.1)),
                    lty=rep(1, 3), bty='n', cex=1.1, col=c("black") );

text(x=0.5, y=0.95, labels=plot_title, cex=1.5, font=2, col="black", xpd=NA)
dev.off()
