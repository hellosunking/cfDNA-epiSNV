argv = commandArgs(T);
if( length(argv) != 1 ) {
	print( 'usage: R --slave --args <sid> < plot.R' );
	q();
}

sid=argv[1];

outfileName = paste0(sid, ".size.vs.end.pdf");
pdf( outfileName, width=10, height=10 );

## load data
altU  = read.table(paste(sid,  ".Mut.all.U.dist", sep=""));
altD  = read.table(paste(sid,  ".Mut.all.D.dist", sep=""));
refU = read.table(paste(sid, ".Wt.all.U.dist", sep=""));
refD = read.table(paste(sid, ".Wt.all.D.dist", sep=""));

## normalize data
altU$V2  = altU$V2  / sum( altU$V2 )*100;
altD$V2  = altD$V2  / sum( altD$V2 )*100;
refU$V2 = refU$V2 / sum(refU$V2 )*100;
refD$V2 = refD$V2 / sum(refD$V2 )*100;

## plot 4: merged alt vs ref
plot( -1000, xlim=c(-100, 100), ylim=range(altU$V2, altD$V2, refU$V2, refD$V2),
	 xlab="Relative position to nucleosome center (bp)", ylab="Normalized end density",
	 main=paste( sid, ": alt vs ref fragments", sep=""));
abline( v=c(-73, -50, 50, 73), col='grey', lty=2);
lines(altD$V2~altD$V1, col='purple',       lwd=2);
lines(refD$V2~refD$V1, col='yellow3', lwd=2);
lines(altU$V2~altU$V1, col='red',  lwd=2);
lines(refU$V2~refU$V1, col='blue', lwd=2);
legend( 'top', c('Mut_DNA - D end', 'Wt_DNA - D end', 'Mut_DNA - U end', 'Wt_DNA - U end'),
	   col=c('purple','yellow3', 'red', 'blue'), lty=c(1,1,1,1), bty='n', cex=1.4, lwd=2 );

