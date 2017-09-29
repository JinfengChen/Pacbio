plot_genotype <- function(bedfile, y1, y2, title){
   gt <- read.table(bedfile)
   xleft <- gt[,2]/1000000
   xright <- gt[,3]/1000000
   ybottom <- rep(y1, length(xleft))
   ytop <- rep(y2, length(xleft))
   col <- as.vector(gt[,5])
   rect(xleft, ybottom, xright, ytop, col=col, border=NA)
   text(-5, y1+1, adj=1, labels=title, xpd=TRUE, cex=1.4)
}

pdf("FAIRCHILD.1chr.GATK.plot_genotype.pdf", width = 10, height = 8)
par(mfrow=c(2, 1))
par(mar=c(4,10,4,2))
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(1.5, 2.5))
#het SNP density
x <- read.table("FCM.vcf.tab.het.1chr.GATK.density", header=TRUE)
y <- read.table("Fairchild_v1.fasta.1chr.chrlen")
plot(x[,2]/1000000, x[,5], col='cadetblue', type='l', main="", xlab="Chromosomal Position (Mb)",  ylab="Heterozygous site/kb", xlim=c(0, 350), ylim=c(0, 12), frame.plot=FALSE, cex=1.4, cex.axis=1.4, cex.lab=1.4)
for (i in y[,2]){
abline(v=i/1000000, col='black', lty=3)
}
legend(240, 13, "Heterozygous SNP", col="cadetblue", lty=1, bty='n', cex=1.4)

#snp <- read.table("FAIRCHILD.vcf.tab.het.1chr.GATK.density", header=TRUE)
#te  <- read.table("Clementine.SlidingWin.Repeat.1chr.density", header=TRUE)
#chrs <- read.table("Cclementina_v1.0_scaffolds.chrlen.1chr")
#plot(te[,3]/1000000-0.1, te[,10], col='black', type='l', main="", xlab="",  ylab="Proportion", xlim=c(0, 290), ylim=c(0, 1), frame.plot=FALSE, cex=1.4, cex.axis=1.4, cex.lab=1.4)
#for (i in chrs[,2]){
#abline(v=i/1000000, col='black', lty=3)
#}
#legend(220, 1.1, "Transposable elements", col="black", lty=1, bty='n', cex=1.4)

#par(mar=c(20,10,1,2))
#plot(snp[,2]/1000000, snp[,5], col='cadetblue', type='l', main="", xlab="Chromosomal Position (Mb)",  ylab="Heterozygous site/kb", xlim=c(0, 290), ylim=c(0, 20), frame.plot=FALSE, cex=1.4, cex.axis=1.4, cex.lab=1.4)
#for (i in chrs[,2]){
#abline(v=i/1000000, col='black', lty=3)
#}
#legend(220, 22, "Heterozygous SNP", col="cadetblue", lty=1, bty='n', cex=1.4)
plot(c(0, 350), c(0, -50), type= "n", xlab = "", ylab = "", frame.plot=FALSE, axes=FALSE)
#GT.10x_up.1chr.txt
plot_genotype('GT.10x_up.1chr.txt', -5, -3, '10x_up')
plot_genotype('GT.10x_down.1chr.txt', -10, -8, '10x_down')

plot_genotype('GT.Pollen41.1chr.txt', -20, -18, 'Pollen41')
plot_genotype('GT.Pollen42.1chr.txt', -25, -23, 'Pollen42')
plot_genotype('GT.Pollen43.1chr.txt', -30, -28, 'Pollen43')
plot_genotype('GT.Pollen44.1chr.txt', -35, -33, 'Pollen44')
plot_genotype('GT.Pollen45.1chr.txt', -40, -38, 'Pollen45')
plot_genotype('GT.hap1.1chr.txt', -45, -43, 'Haplotype1')
plot_genotype('GT.hap2.1chr.txt', -50, -48, 'Haplotype2')
#Fairchild_v1.fasta.1chr.chrlen
chrs <- read.table("Fairchild_v1.fasta.1chr.chrlen")
for (i in chrs[,2]){
abline(v=i/1000000, col='black', lty=3)
}

axis(1, seq(0, 350, by=50), line=0, labels=c("","","","","","","",""), cex=1.4)
text(seq(0, 350, by=50), rep(-55,8), labels=seq(0, 350, by=50), xpd=TRUE, cex=1.4)
text(175, -60, labels="Chromosomal Position (Mb)", xpd=TRUE, cex=1.4)

#gt <- read.table("FAIRCHILD.ancester.uniq.1chr.bed")
#xleft <- gt[,2]/1000000
#xright <- gt[,3]/1000000
#ybottom <- rep(-30, length(xleft))
#ytop <- rep(-25, length(xleft))
#col <- as.vector(gt[,5])
#rect(xleft, ybottom, xright, ytop, col=col, border=NA, xpd=TRUE)
dev.off()

