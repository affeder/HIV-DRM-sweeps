#This script produces figure 2 (F2, F2-S1, F2-S2) and part of figure 5 (F5-S1, F5-S2)

###########################################################################
## Fig 2. How does having a DRM affect diversity?
###########################################################################

##Don't forget to set the working directory! 

weight <- 1.5 # This is for the linewidth for the figures
source("read-in-data.r") 
drminf <- read.table("../dat/DRM_file_WHO.txt", header = T, stringsAsFactors = F)

#######################################
## First, look at RT
#######################################

#What are the most common reverse transcriptase mutations in our dataset?
mutlist <- names(sort(table(unlist(hardDRMlist)[grep("RT", unlist(hardDRMlist))]), decreasing = T)[1:25])
mutpos <- as.numeric(substr(mutlist, start = 3, stop = nchar(mutlist) - 1))

mutinfoRT<-data.frame(mutlist=mutlist, mutpos=mutpos) 

#What are the ancestral states at these positions? 

for(i in 1:length(mutpos)){
    mutinfoRT$ances[i] = unique(drminf[drminf$Position == mutpos[i] & drminf$Location == "RT", ]$Pre)
}

#How many derived versus ancestral patients are there?
for(i in 1:length(mutpos)){
    mutinfoRT$der[i]<-length(which(unlist(lapply(hardDRMlist, function(x){length(grep(mutinfoRT$mutlist[i], x)) > 0}))== TRUE & dat$DRMnum == 1))
    mutinfoRT$anc[i]<- length(which(dat$DRMnum == 0 & ncaaalign[,mutinfoRT$mutpos[i]] == mutinfoRT$ances[i]))
}

#We want to test whether patients with a specific DRM (but no other DRMs) have fewer amb reads than patients with no DRMs. 
#We consider only mutations that occur at least 5 times.

mutinfoRT$enough<-as.numeric(mutinfoRT$der >= 5)

mutinfoRT$tt.estimate1<-NA
mutinfoRT$tt.estimate2<-NA
mutinfoRT$p.value<-NA
mutinfoRT$conf.intL<-NA
mutinfoRT$conf.intR<-NA

for(i in 1:length(mutinfoRT[,1])){
    if (mutinfoRT$enough[i]==1){
    #Patients with derived state at the DRM in question AND no other drms - ambiguous or not
    der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutinfoRT$mutlist[i], x)) > 0})) == TRUE & dat$DRMnum == 1)
    #What patients do not have ANY fixed DRMs #and have nonambiguous ancestral at the site in question
    anc <- which(dat$DRMnum == 0 & ncaaalign[,mutinfoRT$mutpos[i]] == mutinfoRT$ances[i])
    #Only allow well defined indices
    sto <- (list(dat[anc,]$ambnum, dat[der,]$ambnum))
    tt <- t.test(sto[[1]], sto[[2]])
    mutinfoRT$tt.estimate1[i]<-tt$estimate[1]
    mutinfoRT$tt.estimate2[i]<-tt$estimate[2]
    mutinfoRT$p.value[i]<-tt$p.value
    mutinfoRT$conf.intL[i]<-tt$conf.int[1]
    mutinfoRT$conf.intR[i]<-tt$conf.int[2]
    print(mutinfoRT$mutlist[i])
    print(tt)
}}

#######################################
## Second, look at Protease
#######################################

mutlist.new <- names(sort(table(unlist(hardDRMlist)[grep("PI", unlist(hardDRMlist))]), decreasing = T)[1:25])
mutpos.new <- as.numeric(substr(mutlist.new, start = 3, stop = nchar(mutlist.new) - 1))

mutinfoPI<-data.frame(mutlist=mutlist.new, mutpos=mutpos.new)

#What are the ancestral states at these positions? Use the b look up chart
for(i in 1:length(mutinfoPI$mutpos)){
    mutinfoPI$ances[i] = unique(drminf[drminf$Position == mutinfoPI$mutpos[i] & drminf$Location == "PI", ]$Pre)
}

for(i in 1:length(mutinfoPI$mutpos)){
    mutinfoPI$der[i]<-length(which(unlist(lapply(hardDRMlist, function(x){length(grep(mutinfoPI$mutlist[i], x)) > 0}))== TRUE & dat$DRMnum == 1))
    mutinfoPI$anc[i]<- length(which(dat$DRMnum == 0 & ncpraaalign[,mutinfoPI$mutpos[i]] == mutinfoPI$ances[i]))
}

# As before, 5 mutations (without other mutations) must be necessary to proceed)

mutinfoPI$enough<-as.numeric(mutinfoPI$der >= 5)

mutinfoPI$tt.estimate1<-NA
mutinfoPI$tt.estimate2<-NA
mutinfoPI$p.value<-NA
mutinfoPI$conf.intL<-NA
mutinfoPI$conf.intR<-NA

for(i in 1:length(mutinfoPI[,1])){
    if (mutinfoPI$enough[i]==1){
        print (i)
        #Patients with derived state at the DRM in question AND no other drms - ambiguous or not
        der <- which(unlist(lapply(hardDRMlist, function(x){length(grep(mutinfoPI$mutlist[i], x)) > 0})) == TRUE & dat$DRMnum == 1)
        #What patients do not have ANY fixed DRMs #and have nonambiguous ancestral at the site in question
        anc <- which(dat$DRMnum == 0 & ncpraaalign[,mutinfoPI$mutpos[i]] == mutinfoPI$ances[i])
        #Only allow well defined indices
        sto <- (list(dat[anc,]$ambnum, dat[der,]$ambnum))
        tt <- t.test(sto[[1]], sto[[2]])
        mutinfoPI$tt.estimate1[i]<-tt$estimate[1]
        mutinfoPI$tt.estimate2[i]<-tt$estimate[2]
        mutinfoPI$p.value[i]<-tt$p.value
        mutinfoPI$conf.intL[i]<-tt$conf.int[1]
        mutinfoPI$conf.intR[i]<-tt$conf.int[2]
        print(mutinfoPI$mutlist[i])
        print(tt)
    }
}

###########################################################################
## Figure 2A 
###########################################################################

mutinfo<-rbind(mutinfoRT,mutinfoPI)
write.csv(mutinfo, file = "../tmp/mutinfo.csv")

mutinfo<-mutinfo[mutinfo$enough==1,]
pdf("../figures/F2.pdf", width = 10, height = 5)

layout(matrix(1:2, nrow = 1))
par(mar = c(5,4,1,1))
means <- mutinfo$tt.estimate2 - mutinfo$tt.estimate1
plot(0, type = "n", ylim = c(-10,10), xlim = c(1, length(means)), axes = F, ylab =  "Change in diversity (Ambiguous calls)", xlab = "", main = "")
abline(h = 0, lty = "dashed", col = "grey")
text(1.5, 9.5, "A", cex = 2.5)
points(means, cex = 1, bg = "grey", pch = 21, lwd = weight)
arrows(1:length(means), -mutinfo$conf.intL, 1:length(means), -mutinfo$conf.intR, length = 0, lwd = weight)
axis(2)
axis(1,at = 1:length(means), labels = mutinfo$mutlist, las = 2)
box()

###########################################################################
## Figure 2B:  How does having multiple DRMs affect diversity? 
###########################################################################

liststo <- list()
liststo <- append(liststo, rep(0, 11))
for(i in 0:10){
    liststo[[i+1]] <- dat[which(dat$DRMnum ==  i),]$ambnum
}
means <- (unlist(lapply(liststo, mean))) # get mean # ambig sites given # DRMs
ses <- unlist(lapply(liststo, function(x){c(sd(x)/sqrt(length(x)))} )) # get standard error given # ambig sites given # DRMs
plot(0:10, means, ylim = c(min(means - ses), max(means + ses)), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", pch = 21, bg = "grey", lwd = weight)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
text(9.85, 7.2, "B", cex = 2.5)

dev.off()

EffectNumDRMs<-data.frame(NumDRMs=0:10, meanNumAmb=means, stErr= ses)
write.csv(EffectNumDRMs, file = "../tmp/EffectNumDRMs.csv")

########################################################
# Stats for #DRMs vs #amb reads
########################################################

for( i in 0:5){
    g1 <- dat[which(dat$DRMnum == i),'ambnum']
    g2 <- dat[which(dat$DRMnum > i),'ambnum']
    print(paste(i, " versus >", i, ", p value:", t.test(g1, g2)$p.value))
}

#[1] "0  versus > 0 , p value: 3.21768723697155e-28"
#[1] "1  versus > 1 , p value: 8.20597492921639e-15"
#[1] "2  versus > 2 , p value: 5.51367095155348e-05"
#[1] "3  versus > 3 , p value: 0.0503823725393725"
#[1] "4  versus > 4 , p value: 0.0206617461111769"
#[1] "5  versus > 5 , p value: 0.464019222571229"

for( i in 0:5){
    g1 <- dat[which(dat$DRMnum == i),'ambnum']
    g2 <- dat[which(dat$DRMnum == i+1),'ambnum']
    print(paste(i, " versus ", i+1, ", p value:", t.test(g1, g2)$p.value))
}

#[1] "0  versus  1 , p value: 0.000719156237522984"
#[1] "1  versus  2 , p value: 1.6246679696995e-05"
#[1] "2  versus  3 , p value: 0.0107565652917854"
#[1] "3  versus  4 , p value: 0.00245048403240502"
#[1] "4  versus  5 , p value: 0.162466953226026"
#[1] "5  versus  6 , p value: 0.29987903523199"

####################################
# Fig F2-S2 Does this decrease in diversity look different for different subtypes? 
####################################

pdf("../figures/F2-S2.pdf", width =10, height = 6)

subtypenames <- names(sort(table(dat$Subtype), decreasing = T))[1:8]

layout(matrix(1:8, ncol = 4, byrow = T))
for(j in 1:length(subtypenames)){
    print(subtypenames[j])
    liststo <- list()
    liststo <- append(liststo, rep(0, 11))
    subspec <- which(dat$Subtype == subtypenames[j])

    for(i in 0:10){
        inds <- intersect(which(dat$DRMnum ==  i), subspec)
        liststo[[i+1]] <- dat[inds,]$ambnum
    }

    means <- (unlist(lapply(liststo, mean)))
    ses <- unlist(lapply(liststo, function(x){c(sd(x)/sqrt(length(x)))} ))
    toPlot <- which(!is.na(means))

    print(toPlot)
    
    means <- means[toPlot]
    ses <- ses[toPlot]
    
    #write these summaries to files
    EffectNumDRMs_subtype<-data.frame(NumDRMs=toPlot-1, meanNumAmb=means, stErr= ses)
    write.csv(EffectNumDRMs_subtype, file = paste("../tmp/EffectNumDRMs_subtype",subtypenames[j],".csv",sep=""))
    
    plot(toPlot-1, means, ylim = c(0, 20), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = paste("Subtype ", subtypenames[j], sep = ""))
    arrows(toPlot-1, means - ses, toPlot-1, means + ses, length = 0)  
}

dev.off()


##########################################
#What is the distribution of DRMs for different subtypes
##########################################

pdf("../figures/F5-S1.pdf", width =12, height = 8)

par(oma = c(3, 3, 0, 0))
par(mar = c(2,2,4,2))
layout(matrix(1:length(drugnames), ncol = 7))
for(i in 1:length(drugnames)){
    hist(dat[dat$Regimen == drugnames[i], ]$DRMnum, main = paste(drugnames[i]), breaks = seq(-.5, 10.5, by = 1))
}

mtext("Number of DRMs", side = 1, line = 1, outer = TRUE)
mtext("Count", side = 2, line = 1, las = 0, outer = TRUE)
dev.off()

####################################
# F2-S1 Does this decrease in diversity look different for different types of drugs? 
####################################

pdf("../figures/F2-S1.pdf", width =8, height = 8)
#must be run after F5
nrtis.l <- refs.filt[coder > 0 & coder < 4] # list of NRTI treatments
nnrtis.l <- refs.filt[coder==  4] #list of NNRTI treatments
pi.l <- refs.filt[coder == 5] #list of PI treatments
pir.l <- refs.filt[coder == 6] #list of PI/r treatments
nrti.b <- c() #create list of sequences from patients on NRTI treatments
for(i in nrtis.l){
    nrti.b <- c(nrti.b, which(dat$Regimen == i))
}
nnrti.b <- c()#create list of sequences from patients on NNRTI treatments
for(i in nnrtis.l){
    nnrti.b <- c(nnrti.b, which(dat$Regimen == i))
}
pi.b <- c()#create list of sequences from patients on PI treatments
for(i in pi.l){
    pi.b <- c(pi.b, which(dat$Regimen == i))
}
pir.b <- c()#create list of sequences from patients on PI/r treatments
for(i in pir.l){
    pir.b <- c(pir.b, which(dat$Regimen == i))
}
liststo <- list()
nrtisto <- append(liststo, rep(0, 11))
nnrtisto <- append(liststo, rep(0, 11))
pisto <- append(liststo, rep(0, 11))
pirsto <- append(liststo, rep(0, 11))
for(i in 0:10){ #to make plot from 0 to 10 DRMs 
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
yval.letter <- 9.5
layout(matrix(1:4, ncol = 2))
par(mar = c(4,4,2,1))

##NRTI treatments
means <- (unlist(lapply(nrtisto, mean)))
ses <- unlist(lapply(nrtisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "NRTI (1, 2 or 3)", pch = 21, bg = "grey", lwd = weight)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)
text(.4, yval.letter, "A", cex = 2.5)

EffectNumDRMs_treatment<-data.frame(NumDRMs=0:(length(means)-1), meanNumAmb=means, stErr= ses)
write.csv(EffectNumDRMs_treatment, file = paste("../tmp/EffectNumDRMs_treatment_NRTI.csv",sep=""))

##NNRTI treatments
means <- (unlist(lapply(nnrtisto, mean)))
ses <- unlist(lapply(nnrtisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI + NNRTI", pch = 21, bg = "grey", lwd = weight)
text(.4, yval.letter, "C", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)

EffectNumDRMs_treatment<-data.frame(NumDRMs=0:(length(means)-1), meanNumAmb=means, stErr= ses)
write.csv(EffectNumDRMs_treatment, file = paste("../tmp/EffectNumDRMs_treatment_NNRTI.csv",sep=""))

##PI treatments
means <- (unlist(lapply(pisto, mean)))
ses <- unlist(lapply(pisto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI + PI", pch = 21, bg = "grey", lwd = weight)
text(.4, yval.letter, "B", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)

EffectNumDRMs_treatment<-data.frame(NumDRMs=0:(length(means)-1), meanNumAmb=means, stErr= ses)
write.csv(EffectNumDRMs_treatment, file = paste("../tmp/EffectNumDRMs_treatment_PI.csv",sep=""))

##PI/r treatments
means <- (unlist(lapply(pirsto, mean)))
ses <- unlist(lapply(pirsto, function(x){c(sd(x)/sqrt(length(x)))} ))
means[is.na(means)] <- 0
ses[is.na(ses)] <- 0
plot(0:10, means, ylim = c(uppery, lowery), xlab = "Number of DRMs", ylab = "Number of ambiguous calls", main = "2NRTI+PI/r", pch = 21, bg = "grey", lwd = weight)
text(9.6, yval.letter, "D", cex = 2.5)
arrows(0:length(means), means - ses, 0:length(means), means + ses, length = 0, lwd = weight)

EffectNumDRMs_treatment<-data.frame(NumDRMs=0:(length(means)-1), meanNumAmb=means, stErr= ses)
write.csv(EffectNumDRMs_treatment, file = paste("../tmp/EffectNumDRMs_treatment_PIr.csv",sep=""))

dev.off()

###########################################################################
#F5-S2 Information used in the introduction
###########################################################################

years <- sort(unique(dat$IsolateYear))
onedrug <- c()
twodrug <- c()
threedrug <- c()
moredrug <- c()
samplesize <- c()

for(i in 1:length(years)){
    tmp2 <-  dat[dat$IsolateYear == years[i],]$numMeds
    
    samplesize[i] <- length(tmp2)
    onedrug[i] <- sum(tmp2 == 1)/length(tmp2)
    twodrug[i] <- sum(tmp2 == 2)/length(tmp2)
    threedrug[i] <- sum(tmp2 >= 3)/length(tmp2)
    moredrug[i] <- sum(tmp2 > 3)/length(tmp2)
}

pdf("../figures/F5-S2.pdf", width = 9, height = 3.5)
layout(matrix(1:3, ncol = 3))
par(mar = c(4,4,3,1))

#######
## A ## 
#######
cols <- c(rgb(202,0,32, max = 255), rgb(230,97,1, max = 255),rgb(5,113,176, max = 255))
barplot(t(cbind(onedrug, twodrug, threedrug)), names = years, las = 2, col = cols, ylab = "Proportion of patients treated with X drugs")
text(1.4, .95, "A", col = "white", cex = 2.5)

#now a cummulative plot of all the different treatments
uniqreg <- c()
for(i in 1:length(years)){
    uniqreg[i] <-  length(unique(dat[dat$IsolateYear <= years[i],]$Regimen))
}

#######
## B ##
#######
plot(years, uniqreg, type = "l", xlab = "Year", ylab = "Cummulative Unique Regimens", lwd = 2)
text(1989.5, 200, "B", col = "black", cex = 2.5)

#How many ambiguous reads are in each year #looks good. 
ambperyear <- list()
fullyears <- seq(min(years), max(years), by = 1)
for(i in 1:length(fullyears)){
    ambperyear[[i]] <- dat[dat$IsolateYear == fullyears[i],]$ambig
}

t.test(dat[dat$IsolateYear <1995,]$ambig,dat[dat$IsolateYear >=1995,]$ambig)
print(t.test(dat[dat$IsolateYear <1995,]$ambig,dat[dat$IsolateYear >=1995,]$ambig)$p.value)
#[1] 1.081953e-17

#Relationship between year and number of ambiguous reads?
summary(lm(dat$ambig ~ dat$IsolateYear))

#Call:
#    lm(formula = dat$ambig ~ dat$IsolateYear)
#
#Residuals:
#    Min       1Q   Median       3Q      Max 
#-0.01955 -0.01752 -0.00679  0.00955  0.33624 
#
#Coefficients:
#    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)      1.950e-01  1.157e-01   1.685    0.092 .
#dat$IsolateYear -8.819e-05  5.771e-05  -1.528    0.127  
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.02109 on 6715 degrees of freedom
#Multiple R-squared:  0.0003476,    Adjusted R-squared:  0.0001988 
#F-statistic: 2.335 on 1 and 6715 DF,  p-value: 0.1265


#What if we look among sequences with no DRMs?
summary(lm(dat[dat$DRMnum == 0,]$ambig ~ dat[dat$DRMnum == 0,]$IsolateYear))


## Call:
## lm(formula = dat[dat$DRMnum == 0, ]$ambig ~ dat[dat$DRMnum == 
##     0, ]$IsolateYear)

## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.031315 -0.021243 -0.005618  0.013635  0.131508 

## Coefficients:
##                                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                        -1.4168609  0.2629257  -5.389 8.14e-08 ***
## dat[dat$DRMnum == 0, ]$IsolateYear  0.0007194  0.0001312   5.483 4.85e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 0.02486 on 1608 degrees of freedom
## Multiple R-squared:  0.01835,	Adjusted R-squared:  0.01774 
## F-statistic: 30.06 on 1 and 1608 DF,  p-value: 4.849e-08

#How many more ambiguous reads is an increase in year associated with,
# on average, per sequence?
0.0007 * 800

#######
## C ##
#######
boxplot( ambperyear, xaxt = "n",  xlab = "Year", ylab = "% D-PCR amino acid calls ambiguous")
axis(1, at = (2 + seq(0,20, by = 5)), labels = seq(1990, 2010, by = 5))
text(1.2, .345, "C", col = "black", cex = 2.5)

dev.off()


