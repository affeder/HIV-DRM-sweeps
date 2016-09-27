#This script produces figure 2 (F2, F2-S1, F2-S2) and part of figure 5 (F5-S1, F5-S2)

###########################################################################
## How does having a DRM affect diversity?
###########################################################################

weight <- 1.5

drminf <- read.table("~/Dropbox/HIV_softsweeps/DRM_file_WHO.txt", header = T, stringsAsFactors = F)

#What are the most common reverse transcriptase mutations?
mutlist <- names(sort(table(unlist(hardDRMlist)[grep("RT", unlist(hardDRMlist))]), decreasing = T)[1:25])
mutpos <- as.numeric(substr(mutlist, start = 3, stop = nchar(mutlist) - 1))

#What are the ancestral states at these positions? Use the b look up chart
ances <- c()
for(i in mutpos){
    ances <- c(ances, paste(unique(drminf[drminf$Position == i & drminf$Location == "RT", ]$Pre)))
}


#How many derived versus ancestral patients are there?
deranc <- matrix(data = NA, ncol = 2, nrow = length(mutlist))
for(i in 1:length(mutpos)){

    #How many have the derived state
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist[i], x)) > 0})) == TRUE & dat$DRMnum == 1)

    dat$DRMnum[15]
    deranc[i, 1] <- length(der)

    #How many have the ancestral state (and no DRMs)
    anc <- which(dat$DRMnum == 0 & ncaaalign[,mutpos[i]] == ances[i])
    deranc[i, 2] <- length(anc)
}


   #Let's say that we need at least 5 derived mutations in order to be included?
enough <- as.numeric(deranc[,1] >= 5)

mutlist <- mutlist[which(enough == 1)]
mutpos <- mutpos[which(enough == 1)]
ances <- ances[which(enough == 1)]
meanmat <- matrix(data = NA, ncol = 5, nrow = length(mutpos))

cbind(mutlist, ances)

for(i in 1:length(mutpos)){

#Patients with derived state at the DRM in question AND no other drms - ambiguous or not
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist[i], x)) > 0})) == TRUE & dat$DRMnum == 1)
i
    #What patients do not have ANY fixed DRMs #and have nonambiguous ancestral at the site in question
    anc <- which(dat$DRMnum == 0 & ncaaalign[,mutpos[i]] == ances[i])

    #Only allow well defined indices
    sto <- (list(dat[anc,]$ambnum, dat[der,]$ambnum))
    tt <- t.test(sto[[1]], sto[[2]])
    meanmat[i,] <- c(tt$estimate, tt$p.value, tt$conf.int[1:2])

}
meanmat

mutlist.new <- names(sort(table(unlist(hardDRMlist)[grep("PI", unlist(hardDRMlist))]), decreasing = T)[1:25])
mutpos.new <- as.numeric(substr(mutlist.new, start = 3, stop = nchar(mutlist.new) - 1))

mutlist.new
mutpos.new

#What are the ancestral states at these positions? Use the b look up chart
ances <- c()
for(i in mutpos.new){
    ances <- c(ances, paste(unique(drminf[drminf$Position == i & drminf$Location == "PI", ]$Pre)))
}

deranc <- matrix(data = NA, ncol = 2, nrow = length(mutlist.new))

for(i in 1:length(mutpos.new)){

    #How many have the derived state
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist.new[i], x)) > 0})) == TRUE & dat$DRMnum == 1)

    deranc[i, 1] <- length(der)

    #How many have the ancestral state (and no DRMs)
    anc <- which(dat$DRMnum == 0 & ncpraaalign[,mutpos.new[i]] == ances[i])
    deranc[i, 2] <- length(anc)
}


#Let's say that we need at least 5 derived mutations in order to be included?
enough <- as.numeric(deranc[,1] >= 5)

mutlist.new <- mutlist.new[which(enough == 1)]
mutpos.new <- mutpos.new[which(enough == 1)]
ances <- ances[which(enough == 1)]
#meanmat <- matrix(data = NA, ncol = 5, nrow = length(mutpos))
#stos <- matrix(data = NA, ncol = 2, nrow = length(mutpos))


stos.new <- matrix(data = NA, ncol = 2, nrow = length(mutpos.new))
for(i in 1:length(mutpos.new)){
print(i)
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist.new[i], x)) > 0})) == TRUE & dat$DRMnum == 1)
    der.p <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist.new[i], x)) > 0})) == TRUE)

    stos.new[i,1] <- length(der)
    stos.new[i,2] <- length(der.p)
}


cbind(mutlist.new, stos.new)

for(i in 1:length(mutpos.new)){

#Patients with derived state at the DRM in question AND no other drms - ambiguous or not
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutlist.new[i], x)) > 0})) == TRUE & dat$DRMnum == 1)

    #What patients do not have ANY fixed DRMs #and have nonambiguous ancestral at the site in question
    anc <- which(dat$DRMnum == 0 & ncpraaalign[,mutpos.new[i]] == ances[i])

    sto <- (list(dat[anc,]$ambnum, dat[der,]$ambnum))
    tt <- t.test(sto[[1]], sto[[2]])
    meanmat <- rbind(meanmat, c(tt$estimate, tt$p.value, tt$conf.int[1:2]))

}


pdf("../figures/F2.pdf", width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
par(mar = c(5,4,1,1))
means <- meanmat[,2] - meanmat[,1]  
plot(0, type = "n", ylim = c(-10,10), xlim = c(1, nrow(meanmat)), axes = F, ylab =  "Change in diversity (Ambiguous calls)", xlab = "", main = "")
abline(h = 0, lty = "dashed", col = "grey")
text(1.5, 9.5, "A", cex = 2.5)
points(means, cex = 1, bg = "grey", pch = 21, lwd = weight)
arrows(1:nrow(meanmat), -meanmat[,4], 1:nrow(meanmat), -meanmat[,5], length = 0, lwd = weight)
axis(2)
axis(1,at = 1:nrow(meanmat), labels = c(mutlist, mutlist.new), las = 2)
box()

###########################################################################
## How does having multiple DRMs affect diversity?
###########################################################################

gotPI <- grep(PI, dat$Regimen)
gotNoPI <- grep(PI, dat$Regimen, invert = T)
liststo <- list()
noPIsto <- append(liststo, rep(0, 11))
PIsto <- append(liststo, rep(0, 11))
for(i in 0:10){
    inds <- intersect(which(dat$DRMnum ==  i), gotPI)
    PIsto[[i+1]] <- dat[inds,]$ambnum
    inds <- intersect(which(dat$DRMnum ==  i), gotNoPI)
    noPIsto[[i+1]] <- dat[inds,]$ambnum
}
liststo <- list()
liststo <- append(liststo, rep(0, 11))
for(i in 0:10){
    inds <- which(dat$DRMnum ==  i)
    liststo[[i+1]] <- dat[inds,]$ambnum
}
means <- (unlist(lapply(liststo, mean)))
ses <- unlist(lapply(liststo, function(x){c(sd(x)/sqrt(length(x)))} ))
plot(0:10, means, ylim = c(min(means - ses), max(means + ses)), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", pch = 21, bg = "grey", lwd = weight)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
text(9.85, 7.2, "B", cex = 2.5)
dev.off()


####################################
# Does this decrease in diversity look different for different subtypes?
####################################

pdf("../figures/F2-S2.pdf", width =10, height = 6)

subtypenames <- names(sort(table(dat$Subtype), decreasing = T))[1:8]
subs <- list()

layout(matrix(1:8, ncol = 4, byrow = T))
for(j in 1:length(subtypenames)){

    liststo <- list()
    liststo <- append(liststo, rep(0, 11))
    subspec <- which(dat$Subtype == subtypenames[j])

    for(i in 0:10){
        inds <- intersect(which(dat$DRMnum ==  i), subspec)
        liststo[[i+1]] <- dat[inds,]$ambnum
    }

    unlist(lapply(liststo, mean))
    means <- (unlist(lapply(liststo, mean)))
    ses <- unlist(lapply(liststo, function(x){c(sd(x)/sqrt(length(x)))} ))
    toPlot <- which(!is.na(means))

    means <- means[toPlot]
    ses <- ses[toPlot]
    #gee
    plot(toPlot, means, ylim = c(0, 20), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = paste("Subtype ", subtypenames[j], sep = ""))
    arrows(toPlot, means - ses, toPlot, means + ses, length = 0)
}

dev.off()


##########################################
#What is the distribution of DRMs for different subtypes
##########################################

pdf("../figures/F5-S1.pdf", width =12, height = 8)

subs <- list()
for(i in 1:length(drugnames)){
length(drugnames)
    subs[[i]] <- dat[dat$Regimen == drugnames[i], ]$DRMnum

}

par(mar = c(2,2,4,2))
layout(matrix(1:length(drugnames), ncol = 7))
for(i in 1:length(drugnames)){
    hist(subs[[i]], main = paste(drugnames[i]), breaks = seq(-.5, 10.5, by = 1))
}

dev.off()




####################################
# Does this decrease in diversity look different for different types of drugs?
####################################

pdf("../figures/F2-S1.pdf", width =8, height = 8)
#must be run after F5
nrtis.l <- refs.filt[coder > 0 & coder < 4]
nnrtis.l <- refs.filt[coder==  4]
pi.l <- refs.filt[coder == 5]
pir.l <- refs.filt[coder == 6]
nrti.b <- c()
for(i in nrtis.l){
    nrti.b <- c(nrti.b, which(dat$Regimen == i))
}
nnrti.b <- c()
for(i in nnrtis.l){
    nnrti.b <- c(nnrti.b, which(dat$Regimen == i))
}
pi.b <- c()
for(i in pi.l){
    pi.b <- c(pi.b, which(dat$Regimen == i))
}
pir.b <- c()
for(i in pir.l){
    pir.b <- c(pir.b, which(dat$Regimen == i))
}
liststo <- list()
nrtisto <- append(liststo, rep(0, 11))
nnrtisto <- append(liststo, rep(0, 11))
pisto <- append(liststo, rep(0, 11))
pirsto <- append(liststo, rep(0, 11))
for(i in 0:10){
    inds <- intersect(which(dat$DRMnum ==  i), nrti.b)
    nrtisto[[i+1]] <- dat[inds,]$ambnum
    inds <- intersect(which(dat$DRMnum ==  i), nnrti.b)
    nnrtisto[[i+1]] <- dat[inds,]$ambnum
    inds <- intersect(which(dat$DRMnum ==  i), pi.b)
    pisto[[i+1]] <- dat[inds,]$ambnum
    inds <- intersect(which(dat$DRMnum ==  i), pir.b)
    pirsto[[i+1]] <- dat[inds,]$ambnum  
}
uppery <- 0
lowery <- 10
                                        #c(min(means - ses), max(means + ses))
yval.letter <- 9.5
layout(matrix(1:4, ncol = 2))
par(mar = c(4,4,2,1))
means <- (unlist(lapply(nrtisto, mean)))
ses <- unlist(lapply(nrtisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "NRTI (1, 2 or 3)", pch = 21, bg = "grey", lwd = weight)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
text(.4, yval.letter, "A", cex = 2.5)
means <- (unlist(lapply(nnrtisto, mean)))
ses <- unlist(lapply(nnrtisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI + NNRTI", pch = 21, bg = "grey", lwd = weight)
text(.4, yval.letter, "C", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
means <- (unlist(lapply(pisto, mean)))
ses <- unlist(lapply(pisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI + PI", pch = 21, bg = "grey", lwd = weight)
text(.4, yval.letter, "B", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
means <- (unlist(lapply(pirsto, mean)))
ses <- unlist(lapply(pirsto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI+PI/r", pch = 21, bg = "grey", lwd = weight)
text(9.6, yval.letter, "D", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)

dev.off()



###########################################################################
#Information used in the introduction
#categorical data analysis
###########################################################################
#Diversity calculations
#Compute the relative diversity levels
#Generate Figure 5, S2

years <- sort(unique(dat$IsolateYear))

onedrug <- c()
twodrug <- c()
threedrug <- c()
moredrug <- c()
samplesize <- c()

for(i in 1:length(years)){
    tmp1 <-  dat[dat$IsolateYear == years[i],]
    tmp2 <-  dat[dat$IsolateYear == years[i],]$numMeds
#    PIsbyyear[i] <- length(grep(PI, tmp1$Regimen))
#    PIsbyyear[i] <- length(grep(PI, tmp1$Regimen))
    
    samplesize[i] <- length(tmp2)
    onedrug[i] <- sum(tmp2 == 1)/length(tmp2)
    twodrug[i] <- sum(tmp2 == 2)/length(tmp2)
    threedrug[i] <- sum(tmp2 >= 3)/length(tmp2)

    moredrug[i] <- sum(tmp2 > 3)/length(tmp2)

}
newPal <- c(rgb(202,0,32, max = 255),rgb(244,165,130, max = 255),rgb(247,247,247, max = 255),rgb(146,197,222, max = 255),rgb(5,113,176, max = 255), rgb(230,97,1, max = 255))
cols <- newPal[c(1, 6, 5)]


pdf("../figures/F5-S2.pdf", width = 9, height = 3.5)

layout(matrix(1:3, ncol = 3))
par(mar = c(4,4,3,1))
cols <- rainbow(3)
cols <- c(rgb(202,0,32, max = 255), rgb(230,97,1, max = 255),rgb(5,113,176, max = 255))
barplot(t(cbind(onedrug, twodrug, threedrug)), names = years, las = 2, col = cols, ylab = "Proportion of patients treated with X drugs")
text(1.4, .95, "A", col = "white", cex = 2.5)


#now, how about a cummulative plot of all the different treatments
uniqreg <- c()
for(i in 1:length(years)){
    uniqreg[i] <-  length(unique(dat[dat$IsolateYear <= years[i],]$Regimen))
}

plot(years, uniqreg, type = "l", xlab = "Year", ylab = "Cummulative Unique Regimens", lwd = 2)
text(1989.5, 200, "B", col = "black", cex = 2.5)

###Part A

#How many ambiguous reads are in each year
ambperyear <- list()
fullyears <- seq(min(years), max(years), by = 1)
for(i in 1:length(fullyears)){
#    ambperyear[[i]] <- dat[dat$IsolateYear == fullyears[i],]$RTnumAmbig/dat[dat$IsolateYear == fullyears[i],]$RTlen
    ambperyear[[i]] <- dat[dat$IsolateYear == fullyears[i],]$ambig
}

#How does year affect the number of ambiguous reads per year
pre95 <- c()
pos94 <- c()
for(i in 1:length(years)){
    if(years[i] < 1995){
        pre95 <- c(pre95, ambperyear[[i]])
    }
    if(years[i] >=1995){
        pos94 <- c(pos94, ambperyear[[i]])
    }
}
t.test(pre95, pos94)

#Relationship between year and number of ambiguous reads?
summary(lm(dat$ambig ~ dat$IsolateYear))
#Ok, so, actually, it's really barely signicant
#What happens if we remove the first few years?
summary(lm(dat[dat$IsolateYear >= 1995,]$ambig ~ dat[dat$IsolateYear >= 1995,]$IsolateYear))
#Ok, no signal. Ok, great.

boxplot( ambperyear, xaxt = "n",  xlab = "Year", ylab = "% D-PCR amino acid calls ambiguous")
axis(1, at = (2 + seq(0,20, by = 5)), labels = seq(1990, 2010, by = 5))
text(1.2, .345, "C", col = "black", cex = 2.5)

dev.off()





