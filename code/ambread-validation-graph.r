#At the nucleotide level

### Within patient (nonclonal) - propambi

###Plot it
pdf("../figures/SX.pdf", width = 10, height = 6)
par(mar = c(4,4,1,1))
layout(matrix(1:2, ncol = 2))
xlimvals <- (10^((-6):(-1)))
ylimvals <- (10^((-4):(-1)))
plot(clonalMeans.nu, propambig.nu, xlab = "Within-patient diversity (clonal dataset)", ylab = "Proportion D-PCR nucleotide sequence calls ambiguous", pch = 16, col = rgb(0,0,0,1), main = "Nucleotide level", cex = .5, log = "xy", axes = F, ylim = range(ylimvals), xlim = range(xlimvals))
axis(1, at = xlimvals, labels = parse(text=paste("10", "^-", 6:1, sep="")))
axis(2, at = ylimvals, labels = parse(text=paste("10", "^-", 4:1, sep="")))
text(xlimvals[1]*1.2, ylimvals[length(ylimvals)]*.8, "A", cex = 2.5)

rval <- round(cor.test(clonalMeans.nu, propambig.nu)$estimate, digits = 2)
text(10^(-4), 10^(-1.85), paste("r = ", rval, sep = ""), cex = 1.75)
box()

cor.test(clonalMeans.nu, propambig.nu)

#range(log(clonalMeans.aa[clonalMeans.aa > 0]))
#range(log(ambfreeentrop.aa[ambfreeentrop.aa > 0]))
#range(log(propambig.aa[propambig.aa > 0]))

cbind(-7:-1, log(10^(-7:-1)))
#At the amino acid level

xlimvals <- (10^((-4):(-1)))
ylimvals <- (10^((-4):(-1)))
plot(clonalMeans.aa, propambig.aa,xlab = "Within-patient diversity (clonal dataset)", ylab = "Proportion D-PCR amino acid sequence calls ambiguous", pch = 16, col = rgb(0,0,0,1), main = "Amino acid level", axes = F, cex = .5, log = "xy", ylim = range(ylimvals), xlim = range(xlimvals))
axis(1, at = xlimvals, labels = parse(text=paste("10", "^-", 4:1, sep="")))
axis(2, at = ylimvals, labels = parse(text=paste("10", "^-", 4:1, sep="")))
text(xlimvals[1]*1.2, ylimvals[length(ylimvals)]*.8, "B", cex = 2.5)
rval <- round(cor.test(clonalMeans.aa, propambig.aa)$estimate, digits = 2)
text(10^(-3.5), 10^(-1.85), paste("r = ", rval, sep = ""), cex = 1.75)

box()
dev.off()

##List the correlation coefficients.

cor.test(clonalMeans.aa, propambig.aa)

