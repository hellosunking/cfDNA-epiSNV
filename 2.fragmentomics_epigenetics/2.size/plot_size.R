argv = commandArgs(T);
if( length(argv) < 1 ) {
	print( 'usage: R --slave --args <sid> [label1=Wt] [label2=Mut] < plot.R' );
	q();
}

label.1="Wt"
label.2="Mut"

if( length(argv) > 1 ) {
	label.1 = argv[2]
	if( length(argv) > 2 ) {
		label.2 = argv[3]
	}
}

sid = argv[1];
outfileName = paste0( sid, ".size.pdf" );
pdf( outfileName );

t = read.table( paste0(sid, ".Mut.size") );
n = read.table( paste0(sid, ".Wt.size") );

t[,2] = t[,2]/sum(t[,2])*100;
n[,2] = n[,2]/sum(n[,2])*100;

plot( n[,2] ~ n[,1], col='blue', lwd=2, xlim=c(50, 250), type='l',
		xlab="cfDNA fragment size (bp)", ylab="Frequency (%)", main=sid );
lines(t[,2] ~ t[,1], col='red', lwd=2 );
legend( 'topleft', c(label.2, label.1), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

plot( n[,2] ~ n[,1], col='blue', lwd=2, xlim=c(50, 600), type='l',
		xlab="cfDNA fragment size (bp)", ylab="Frequency (%)", main=sid );
lines(t[,2] ~ t[,1], col='red', lwd=2 );
legend( 'topleft', c(label.2, label.1), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1 );

