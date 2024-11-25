library( ROCR );
library( pROC );

options( stringsAsFactors=F );

args <- commandArgs(trailingOnly = TRUE)
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

ctr = subset(dat, Group=="0");
hcc = subset(dat, Group=="1");

dat = rbind(ctr, hcc)

## use pROC, which is BETTER than ROCR
roc.1 = roc(dat, Group, sbs_contribution );	## roc(data, response, predictor)
auc.1 = sprintf("%.4f", roc.1$auc );
p.1   = sprintf("%.4e", pROC.p(roc.1));

plot( roc.1, lwd=2, col="black", cex.lab=1.2);
legend( 'bottomright', legend=paste0(c("SBS_contribution"), ", AUC=", c(auc.1), ", P=", c(p.1)),
	            lty=rep(1, 3), bty='n', cex=1.1, col=c("black") );
