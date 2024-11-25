library( ROCR );
library( pROC );

options( stringsAsFactors=F );

outfileName='roc.pdf';
pdf( outfileName );

## use ROCR
plot.ROC <- function( ctr, hcc, message="HCC vs Ctr" ) {
	if( median(ctr) > median(hcc) ) {	## I need HCC > CTR
		ctr = -ctr;
		hcc = -hcc;
	}
	roc = cbind( c(ctr, hcc), c(rep(0,length(ctr)), rep(1,length(hcc))) );
	pred= prediction(roc[,1], roc[,2]);
	perf= performance(pred, "tpr", "fpr");
	auc = performance(pred, "auc");
	auc = sprintf("%.4f", auc@y.values);
	plot( perf, lwd=2, main=paste(message, ", AUC=", auc, sep=""));
	abline( a=0, b=1, col='grey');

	auc;
}

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

dat = read.table( "diff.txt", head=T );

ctr = subset(dat, Group=="0");
hcc = subset(dat, Group=="1");

## ROC using pred
#plot.ROC( ctr$pred, hcc$pred, "pred" );

## use pROC, which is BETTER than ROCR
roc.1 = roc( dat, Group, Diff );	## roc(data, response, predictor)
auc.1 = sprintf("%.4f", roc.1$auc );
p.1   = sprintf("%.4e", pROC.p(roc.1));

#roc.2 = roc( dat, Type, Eindex );	## roc(data, response, predictor)
#auc.2 = sprintf("%.4f", roc.2$auc );
#p.2   = sprintf("%.4e", pROC.p(roc.2));

## ROC comparison
#delong = roc.test( roc.1, roc.2 );

#plot( roc.1, lwd=2, col="blue", cex.lab=1.2, main=paste0("P=", delong$p.value) );
plot( roc.1, lwd=2, col="black", cex.lab=1.2 );
#plot( roc.2, lwd=2, col='red', add=T );
#legend( 'bottomright', legend=paste0(c("pred", "E-index"), ", AUC=", c(auc.1, auc.2), ", P=", c(p.1, p.2)),
	            #lty=rep(1, 3), bty='n', cex=1.1, col=c("blue", "red") );
legend( 'bottomright', legend=paste0(c("Delta_S150"), ", AUC=", c(auc.1), ", P=", c(p.1)),
                    lty=rep(1, 3), bty='n', cex=1.1, col=c("black") );

