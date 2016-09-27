#Subsample plot

pdf("../figures/F3.pdf", width =6, height =6)
#pdf("~/Desktop/elife-figs/new//F3.pdf", width =6, height =6)
par(mar = c(3,4,.5,.5))
layout(matrix(1:4, ncol = 2, byrow = T), widths=c(1,.65), heights=c(2,1))

newPal <- c(rgb(202,0,32, max = 255),rgb(244,165,130, max = 255),rgb(247,247,247, max = 255),rgb(146,197,222, max = 255),rgb(5,113,176, max = 255), rgb(230,97,1, max = 255))

refs.lab <- refs.filt 

refs.lab[refs.filt == "ATV+FTC+RTV+TDF"] <- "FTC+TDF+ATV/r"
refs.lab[refs.filt == "3TC+AZT+LPV"] <- "3TC+AZT+LPV/r"
refs.lab[refs.filt == "3TC+ABC+ATV+RTV"] <- "3TC+ABC+ATV/r"
refs.lab[refs.filt == "3TC+ABC+LPV"] <- "3TC+ABC+LPV/r"

#How many in each category?
table(coder)

#Let's precompute the offset, that will allow us to plot
offset = rep(0, length(coder))
for(i in 1:6){

    #offset.var needs to be tuned in order to make the graph look ok
    offset.var <- .12
    tooff <- which(coder == i)
    if(length(tooff)%%2 == 0){
        val <- length(tooff)/2
        offset[which(coder == i)] <- seq(-val*offset.var - offset.var/2, val*offset.var + offset.var/2, by = offset.var)
    }else{
        val <- floor(length(tooff)/2)
        offset[which(coder == i)] <- seq(-val*offset.var, val*offset.var, by = offset.var)
    }
}
#Warnings are ok here

#direction of the text labeling
dir <- rep(1, length(coder))
dir[ coder == 3 ] <- 1.3
dir[coder == 4] <- 1.65
dir[coder == 5] <- 1.3
dir[coder == 6] <- 1.4
dir[coder == 1] <- .85

#Now, let's take the values from our GLMM
resamp.dat.glmm.all.lte4 <- read.table("../tmp/GLMM.all.lte4.randeffs.txt", header = TRUE)
fe.coef.all.lte4 <- read.table("../tmp/GLMM.all.lte4.fixedeffs.txt")


#resamp.dat.glmm.all.lte4 <- read.table("../../code/dat/GLMM.all.lte4.randeffs.txt", header = TRUE)
#fe.coef.all.lte4 <- read.table("../../code/dat/GLMM.all.lte4.fixedeffs.txt")


#We now convert our old column names to headers to
plotnames <- gsub("\\.", "\\+", colnames(resamp.dat.glmm.all.lte4))
plotnames[which(substr(plotnames, 1, 1) == "X")] <- substr(plotnames[which(substr(plotnames, 1, 1) == "X")], 2, max(nchar(plotnames)))

glmm.means <- apply(resamp.dat.glmm.all.lte4, 2, mean, na.rm = TRUE) + apply(fe.coef.all.lte4, 2, mean)[2]
names(glmm.means) <- plotnames

glmm.intervals <- apply(resamp.dat.glmm.all.lte4, 2, quantile, c(.025, .975), na.rm = TRUE) + apply(fe.coef.all.lte4, 2, mean)[2]
colnames(glmm.intervals) <- plotnames

#Let's reorganize glmm.means to be in the same order as coder/refs.filt...
code.order.convert <- c()
for(i in 1:length(refs.filt)){
    code.order.convert[i] <- which(plotnames == refs.filt[i])
}

glmm.1995.ordered.means <- glmm.means[code.order.convert]
glmm.1995.ordered.intervals <- glmm.intervals[,code.order.convert]


#Plot the RTI portion of this graph
###########################################################################
#PART 3: RTI lines
###########################################################################
#We want slightly different par here. 

#layout(matrix(1:2, nrow = 1))
par(mar = c(4,4,.5,.5))

lowerlim <- -2.2
upperlim = 2.5

adj <- rep(1, length(refs.filt))
adj[coder == 1] <- -1.1
adj[coder == 2] <- .4
adj[coder == 3] <-  .5
adj[coder == 4] <- .7
adj[coder == 5] <- -1.5
adj[coder == 6] <- .6

internallyordered <- list()
for(i in 1:6){
    order(which(coder == i), decreasing = TRUE)
}


plot(coder[coder > 0 & coder < 5], glmm.1995.ordered.means[coder > 0 & coder < 5], axes = F, xlab = "",ylab =  expression(paste("Change in diversity accompanying each DRM (", Delta, "DRM)", sep = "")), type = "n", ylim = c(lowerlim,upperlim), xlim = c(.5, 4.5), col =  newPal[coder[coder > 0 & coder < 5]+2])
abline(h = 0, lty = "dashed", col = "grey")
coltmp <- newPal[c(2,2,2,1)]
for(k in 1:4){
    i <- which(coder ==k)
    new.ord <- order(glmm.1995.ordered.means[i], decreasing = TRUE)
    arrows(coder[i] + offset[i], glmm.1995.ordered.intervals[1,i][new.ord] ,coder[i] + offset[i], glmm.1995.ordered.intervals[2,i][new.ord]  , length = 0, col = coltmp[coder[i]], lwd = 2.5)
    points(coder[i] + offset[i], glmm.1995.ordered.means[i][new.ord], col = coltmp[coder[i]], pch = 16, cex = 1)
    text(coder[i] + offset[i], glmm.1995.ordered.intervals[2,i][new.ord] + adj[i], refs.lab[i][new.ord], srt = 90, cex = .75, col = "black")
}
axis(1, at = 1:4, labels = c("1NRTI", "2NRTI", "3NRTI", " \n2NRTI+\nNNRTI"))
axis(2:4)
#text(.55, .0092, "A", cex = 2.5)
text(.55, 1.95, "A", cex = 2)
box()


###########################################################################
#PART 5: PI lines
###########################################################################
#We want slightly different par here. 
#par(mar = c(4,4, 2, 1))
par(mar = c(4,.5,.5,.5))
plot(coder[coder > 4], glmm.1995.ordered.means[coder > 4], axes = F, xlab = "", ylab =  "", type = "n", ylim = c(lowerlim, upperlim), xlim = c(4.5, 6.25))
abline(h = 0, lty = "dashed", col = "grey")
cols <- newPal[c(4,5)]

for(k in 5:6){
                                        #category accomodation
    i <- which(coder ==k)
    new.ord <- order(glmm.1995.ordered.means[i], decreasing = TRUE)
    arrows(coder[i] + offset[i], glmm.1995.ordered.intervals[1,i][new.ord] ,coder[i] + offset[i],glmm.1995.ordered.intervals[2,i][new.ord] , length = 0, col = cols[coder[i]-4], lwd = 2.5)
    points(coder[i] + offset[i],  glmm.1995.ordered.means[i][new.ord], col = cols[coder[i]-4], pch = 16, cex = 1)
    text(coder[i] + offset[i], glmm.1995.ordered.intervals[2,i][new.ord] + adj[i], refs.lab[i][new.ord], srt = 90, cex = .75, col = "black")

}

axis(1, at = 5:6, labels = c("2NRTI+PI", "2NRTI+PI/r"))
#axis(3:4)

#text(4.55, .0092, "B", cex = 2.5)
text(4.55, 1.95, "B", cex = 2)
box()

sum(dat$IsolateYear > 1994)

add.alpha <- function(col, alpha=1){
apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

#GLM plot (1995+)
#layout(matrix(1:2, nrow = 1))

load("../tmp/NNRTIs.all.lte4")
load("../tmp/NRTIs.all.lte4")

load("../tmp/PIrs.all.lte4")
load("../tmp/PIs.all.lte4")


par(mar = c(4,6.5,0,3.5))
lwd.val <- 2
intervals <- apply(NNRTIs.all.lte4, 2, quantile, c(.025, .975))
means <- apply(NNRTIs.all.lte4, 2, mean)
plot(0, type = "n", xlim = c(0, 4), ylim = c(0, 4), ylab = c("# ambiguous reads"), xlab = c("#DRMs") )
polygon(c(rev(seq(0, 4, length.out = 100)),seq(0, 4, length.out = 100)), c(rev(intervals[1,]), intervals[2,]), col = add.alpha(newPal[1], .25), border = FALSE)
lines(seq(0, 4, length.out = 100), means, col = newPal[1] , lwd = lwd.val )
intervals <- apply(NRTIs.all.lte4, 2, quantile, c(.025, .975))
means <- apply(NRTIs.all.lte4, 2, mean)
polygon(c(rev(seq(0, 4, length.out = 100)),seq(0, 4, length.out = 100)), c(rev(intervals[1,]), intervals[2,]), col = add.alpha(newPal[2], .25), border = FALSE)
lines(seq(0, 4, length.out = 100), means, col = newPal[2], lwd = lwd.val)
intervals <- apply(PIrs.all.lte4, 2, quantile, c(.025, .975))
means <- apply(PIrs.all.lte4, 2, mean)
text(0.14, 3.65, "C", cex = 2)



par(mar = c(4,.5,0,.5))
plot(0, type = "n", xlim = c(0, 4), ylim = c(0, 4), ylab = c("# ambiguous reads"), xlab = c("#DRMs") )
polygon(c(rev(seq(0, 4, length.out = 100)),seq(0, 4, length.out = 100)), c(rev(intervals[1,]), intervals[2,]), col = add.alpha(newPal[4], .25), border = FALSE)
lines(seq(0, 4, length.out = 100),means, col = newPal[4] , lwd = lwd.val)
intervals <- apply(PIs.all.lte4, 2, quantile, c(.025, .975))
means <- apply(PIs.all.lte4, 2, mean)
polygon(c(rev(seq(0, 4, length.out = 100)),seq(0, 4, length.out = 100)), c(rev(intervals[1,]), intervals[2,]), col = add.alpha(newPal[5], .4), border = FALSE)
lines(seq(0, 4, length.out = 100), means, col = newPal[5], lwd = lwd.val)
text(0.14, 3.65, "D", cex = 2)

dev.off()




#Non-parametric, -SH
#rand.effs <- apply(resamp.dat.glmm.all.lte4, 2, mean, na.rm = TRUE) #+ fe.coef.all.lte4
rand.effs <- glmm.means
rand.effs.names <- plotnames#[ordered.regs]
matched.drugnames <- c()
matched.effects <- c()
ordered.regs <- c()
for(i in 1:length(drugnames)){
    toAdd <- which(rand.effs.names == drugnames[i])
    if(length(toAdd) > 0){matched.effects[i] <- toAdd
                          matched.drugnames[i] <- drugnames[i]
                          ordered.regs[i] <- toAdd
                      }else{
        matched.effects[i] <- NA
        matched.drugnames[i] <- NA
        ordered.regs[i] <- NA
    }
}


sig.val <- c()
sig.neg <- c()
allfits <- matrix(data = NA, nrow = nrow(resamp.dat.glmm.all.lte4), ncol = 101)
for(i in 1:nrow(resamp.dat.glmm.all.lte4)){
    fit.i <- resamp.dat.glmm.all.lte4[i,ordered.regs[!is.na(ordered.regs)]] + fe.coef.all.lte4[i, 2]
#cbind(drugnames[-which(is.na(ordered.regs ))], colnames(fit.i))
    fit.i.lm <- lm(unlist(fit.i) ~ percentfail[-which(is.na(ordered.regs))])
    allfits[i,] <- coef(fit.i.lm)[1]+ 0:100 * coef(fit.i.lm)[2]
}



rel.all.inds <- which(dat$DRMnum <= 4)
treatsize <- c()
rel.col <- c()
for(i in matched.drugnames){
    treatsize <- c(treatsize, length(which(dat[rel.all.inds, 'Regimen'] == i)))
    cod.val <- which(refs.filt == i)
    if(length(cod.val) == 0){rel.col <- c(rel.col, "white")
    }else{
        rel.col <- c(rel.col, coder[cod.val])
    }
}
rel.col[rel.col <4] <- newPal[2]
rel.col[rel.col == 4] <- newPal[1]
rel.col[rel.col == 5] <- newPal[4]
rel.col[rel.col == 6] <- newPal[5]





pdf("../figures/F4.pdf", width =6, height =5)
#pdf("~/Desktop/elife-figs/new//F4.pdf", width =6, height =5)

par(mar = c(4,4,1, 1))
plot(percentfail, rand.effs[matched.effects], xlab = "Percentage of patients with virologic suppression after 48 weeks" , ylab =  expression(paste("Change in diversity accompanying each DRM (", Delta, "DRM)", sep = "")), type = "n", xlim = c(0, 105), ylim = c(-2, upperlim))
#abline(lm(rand.effs[matched.effects] ~ percentfail))
rands.perc <- rand.effs[matched.effects]
interval <- apply(allfits, 2, quantile, c(.025, .975), na.rm = TRUE)
polygon( c(0:100, 100:0), c(interval[1,], rev(interval[2,])), col = rgb(0,0,0,.15), border = FALSE)
abline(h = 0, lty = "dashed", col = 'grey')
for(ind in order(treatsize, decreasing = TRUE)){
    points(percentfail[ind], rand.effs[matched.effects][ind], cex = sqrt(treatsize[ind]/10), bg = rel.col[ind], col = "black", pch = 21)
}
offset.4 <- rep(.2, length(matched.drugnames))
offset.4[matched.drugnames == "AZT"] <- .24
offset.4[matched.drugnames == "3TC+AZT"] <- .24
offset.4[matched.drugnames == "3TC+AZT+NVP"] <- .3
offset.4[matched.drugnames == "3TC+AZT+NFV"] <- .26
offset.4[matched.drugnames == "3TC+AZT+EFV"] <- -.32
offset.4[matched.drugnames == "3TC+D4T+NVP"] <- .46
offset.4[matched.drugnames == "AZT+DDI"] <- -.28
offset.4[matched.drugnames == "3TC+ABC+EFV"] <- .16
offset.4[matched.drugnames == "3TC+ABC+AZT"] <- .26
offset.4[matched.drugnames == "3TC+D4T+IDV"] <- .14
offset.4[matched.drugnames == "D4T+DDI+EFV"] <- -.16
offset.4[matched.drugnames == "3TC+D4T+EFV"] <- .08
offset.4[matched.drugnames == "3TC+EFV+TDF"] <- .22
text(percentfail, rand.effs[matched.effects]-offset.4, matched.drugnames, cex = .5)
legend("topright", c('1, 2 or 3 NRTI', '2NRTI+NNRTI', '2NRTI+PI', "2NRTI+PI/r"), pch = c(21, 21, 21, 21), ncol = 1, pt.cex = c( 1, 1, 1, 1), col = "black", pt.bg = c( newPal[2], newPal[1], newPal[4], newPal[5]), box.lwd = 0, cex = .75)
points(c(-2, 9, 27), rep(-1.8,3), pch = c(21, 21, 21), cex = sqrt(c(10, 100, 500)/10), col = "black", bg = "black")
text(c(1.5, 15, 36.5), rep(-1.8, 3), c("10", "100", "500"), cex = .75)
polygon(c(42, -20, -20, 42), c(-3, -3, -1.4, -1.4))

dev.off()


coef(lm(rand.effs[matched.effects] ~ percentfail))
.98/.016
#A ten percent increase in treatment efficacy = .2 amb reads lost with each additional DRM
#these will need to be updated slightly once we run the rest of our random effects 

summary(lm(ambnum ~ IsolateYear, data = dat[dat$DRMnum == 0, ]))


