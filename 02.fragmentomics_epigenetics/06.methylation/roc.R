library( ROCR );
library( pROC );

options( stringsAsFactors=F );

args <- commandArgs(trailingOnly = TRUE)
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

ctr = subset(dat, Group=="0");
COREAD = subset(dat, Group=="2");
ESCA = subset(dat, Group=="3");
STAD = subset(dat, Group=="4");
LIHC = subset(dat, Group=="5");
NSCLC = subset(dat, Group=="6");
PACA = subset(dat, Group=="7");

coread = rbind(ctr,COREAD);
esca = rbind(ctr,ESCA);
stad = rbind(ctr,STAD);
lihc = rbind(ctr,LIHC);
nsclc = rbind(ctr,NSCLC);
paca = rbind(ctr,PACA);

## use pROC, which is BETTER than ROCR
roc.2 = roc( coread, Group, pred );
auc.2 = sprintf("%.4f", roc.2$auc );
p.2   = sprintf("%.4e", pROC.p(roc.2));
roc.3 = roc( esca, Group, pred );
auc.3 = sprintf("%.4f", roc.3$auc );
p.3   = sprintf("%.4e", pROC.p(roc.3));
roc.4 = roc( stad, Group, pred );
auc.4 = sprintf("%.4f", roc.4$auc );
p.4   = sprintf("%.4e", pROC.p(roc.4));
roc.5 = roc( lihc, Group, pred );
auc.5 = sprintf("%.4f", roc.5$auc );
p.5   = sprintf("%.4e", pROC.p(roc.5));
roc.6 = roc( nsclc, Group, pred );
auc.6 = sprintf("%.4f", roc.6$auc );
p.6   = sprintf("%.4e", pROC.p(roc.6));
roc.7 = roc( paca, Group, pred );
auc.7 = sprintf("%.4f", roc.7$auc );
p.7   = sprintf("%.4e", pROC.p(roc.7));

plot( roc.2, lwd=2, col="blue", cex.lab=1.2, main=plot_title );
plot( roc.3, lwd=2, col="deepskyblue", add=T );
plot( roc.4, lwd=2, col="turquoise4", add=T );
plot( roc.5, lwd=2, col="#CD2626", add=T );
plot( roc.6, lwd=2, col="purple4", add=T );
plot( roc.7, lwd=2, col="darkorange4", add=T );

legend( 'bottomright', legend=paste0(c("COREAD", "ESCA", "STAD", "LIHC", "NSCLC", "PACA"), ", AUC=", c(auc.2, auc.3, auc.4, auc.5, auc.6, auc.7), ", P=", c(p.2, p.3, p.4, p.5, p.6, p.7)),
	            lty=1, bty='n', cex=1.1, col=c("blue", "deepskyblue", "turquoise4", "#CD2626", "purple4", "darkorange4") );
