#Resamp file
require(vcd)
require(glmmADMB)
require(parallel)
require(doParallel)


#####################
# Read in the data
#####################

drugnames <- drugnames[order(percentfail)]
percentfail <- percentfail[order(percentfail)]


###################
# Compute subsampling amounts
###################

years <- sort(unique(dat$IsolateYear))
zeroclass <- which(dat$DRMnum == 0)

year.zero.mean <- c()
for(i in 1:length(years)){
    year.zero.mean[i] <- mean(dat[intersect(zeroclass,which(dat$IsolateYear == years[i])),'ambnum'])
}


#Full (i.e., all data, not just 1995+)
fit2b <- lm(dat[zeroclass, 'ambnum'] ~ poly(dat[zeroclass, 'IsolateYear'], 1, raw=TRUE))$coefficients
poly2b <- function(year){c(fit2b[1] + fit2b[2]*year )}
down.sample <- function(year){
    toRet <- poly2b(1989)/poly2b(year)
    names(toRet) <- ""
    c(toRet)
}
fit2b

#1995+
fit2b.1995 <- lm(dat[intersect(zeroclass, which(dat$IsolateYear >=1995)), 'ambnum'] ~ poly(dat[intersect(zeroclass, which(dat$IsolateYear >=1995)), 'IsolateYear'], 1, raw=TRUE))$coefficients
poly2b.1995 <- function(year){c(fit2b.1995[1] + fit2b.1995[2]*year )}
down.sample.1995 <- function(year){
    toRet <- poly2b.1995(1995)/poly2b.1995(year)
    names(toRet) <- ""
    c(toRet)
}
fit2b.1995

#Setup plot F5-S3
pdf("../figures/F5-S3.pdf", width =9, height =5)
par(mar = c(4,4,1,1))
layout(matrix(1:2, nrow = 1))
plot(years, year.zero.mean, ylab = "Number of Ambiguous Calls", xlab = "Year", pch = 16, col = "grey", ylim = c(0, 10))
lines(years, poly2b(years))
points(years, poly2b(years), col = "black", pch = 16)
lines(seq(1995, 2013, by = 1), poly2b.1995(seq(1995, 2013, by = 1)), col = "red")
points(seq(1995, 2013, by = 1), poly2b.1995(seq(1995, 2013, by = 1)), col = "red", pch =16)
legend("bottomright", c("Mean of 0-DRM seqs.", "1989+ model fits", "1995+ model fits"), col = c("grey", "black", "red"), pch = c(16, 16, 16),  lty = c(0, 1, 1))
plot(years, 1/(poly2b(years)/min(poly2b(years))), xlab = "Year", ylab = "Subsample effect", pch = 3, cex = .5)
points(seq(1995, 2013, by = 1), 1/(poly2b(seq(1995, 2013, by = 1))/min(poly2b(seq(1995, 2013, by = 1)))), ylab = "Subsample effect", pch = 3, cex = .5, col = "red")
legend("topright", c("1989+ model", "1995+ model"), col = c("black", "red"),  pch = c( 3, 3))
dev.off()





###################
# Assess model fit
###################


#These are the specific draws that will be plotted
downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
nbfit <- goodfit(downsamp.amb, "nbinomial")
poisfit <- goodfit(downsamp.amb, "poisson")
downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample(dat$IsolateYear))
nbfit.all <- goodfit(downsamp.amb, "nbinomial")
poisfit.all <- goodfit(downsamp.amb, "poisson")

fit.plot <- function(x, maintitle, AorB, AorB.x, AorB.y){

    xlim.vals <- range(x$count)

    xlim.vals[1] <- xlim.vals[1] - 0.5
    xlim.vals[2] <- xlim.vals[2] + 0.5
    plot(0, type = "n", xlim = xlim.vals, ylim = c(min(sqrt(x$fitted) - sqrt(x$observed)), max(sqrt(x$fitted))), xlab = "Number of Ambiguous Calls", ylab = "sqrt(Frequency)", main = maintitle)
    upperlims <- sqrt(x$fitted)
    lowerlims <- sqrt(x$fitted) - sqrt(x$observed)
    offsetval <- .4
    for(i in x$count){
        polygon( c( i - offsetval, i + offsetval, i + offsetval, i - offsetval ), c(rep(lowerlims[i+1],2), rep(upperlims[i+1], 2) ), col = "grey")
    }
    points(x$count, sqrt(x$fitted), pch = 16, col = "red")
    lines(x$count, sqrt(x$fitted), pch = 16, col = "red")
    abline(h = 0)
    text(AorB.x, AorB.y, AorB, cex = 2.5)
}

pdf("../figures/F5-S4.pdf", width =8, height =4.5)
layout(matrix(1:2, nrow = 1))
par(mar = c(4,4,3,1))
fit.plot(poisfit.all, "Poisson Fits", "A", 48, 39)
fit.plot(nbfit.all, "Negative Binomial Fits", "B", 48, 48)
dev.off()

pdf("../figures/F5-S5.pdf", width =8, height =4.5)
layout(matrix(1:2, nrow = 1))
par(mar = c(4,4,3,1))
fit.plot(poisfit, "Poisson Fits", "A", 68, 34.5)
fit.plot(nbfit, "Negative Binomial Fits", "B", 68, 44.5)
dev.off()



##############################
#1995+
#Truncated to 4 DRMs
##############################


#1995+
#Truncated to 4 DRMs
#GLMM, note: GLMM code is quite slow to run
##############################

#Note, in the paper, we run 1000 iterations, but it's quite slow, so we've updated it so that only 20 iterations (parallelized to four cores)
iters<- 20

#Set up code that can be run on 4 cores
cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()
com.treat.inds <- c()

#Only run this for treatments that are abundant
for(i in refs.filt){
    com.treat.inds <- c(com.treat.inds, which(dat$Regimen == i))
}
tmpdat <- dat[com.treat.inds, ]
#The regimens need to be converted to factors
tmpdat$Regimen <- factor(tmpdat$Regimen)


#to be run 'iters' number of times
resamp.glmm.1995.lte4 <-foreach(icount(iters), .packages = 'glmmADMB') %dopar% {
    #somesample the number of ambiguous reads
    tmpdat$downsampamb <- rpois(length(tmpdat$ambnum), tmpdat$ambnum* down.sample.1995(tmpdat$IsolateYear))
                                        #Fit the GLMM
    sto.withzero <- glmmadmb(downsampamb ~ 1 + DRMnum + lens + (1 + DRMnum|Regimen), data = tmpdat[intersect(which(tmpdat$IsolateYear >= 1995), which(tmpdat$DRMnum <= 4)),], family = "nbinom", zeroInflation = FALSE)
    #return the GLMM information (more processing will be necessary)
    c(coef(sto.withzero), sto.withzero[[23]])
}

#Stop using four cores
print(Sys.time()-strt)
stopCluster(cl)

#Record the number of treatments
numTreats <- nrow(resamp.glmm.1995.lte4[[1]]$Regimen)

#Set up a matrix that will keep track of all the random effects 
resamp.dat.glmm.1995.lte4 <- matrix(data = NA, nrow = length(resamp.glmm.1995.lte4), ncol =numTreats)
for(i in 1:length(resamp.glmm.1995.lte4)){
    resamp.dat.glmm.1995.lte4[i, ] <-  (resamp.glmm.1995.lte4[[i]]$Regimen)[,2]
}
colnames(resamp.dat.glmm.1995.lte4) <- rownames(resamp.glmm.1995.lte4[[1]]$Regimen)

resamp.dat.glmm.1995.lte4
#Store the random effects coefficients (so it can be graphed without having to rerun everything)
write.table(resamp.dat.glmm.1995.lte4, "../tmp/GLMM.1995.lte4.randeffs.txt", row.names = FALSE, col.names  = TRUE , quote = FALSE)

#Set up a matrix that will keep track of all the fixed effects
fe.dat.glmm.1995.lte4 <- matrix(data = NA, nrow = length(resamp.glmm.1995.lte4), ncol = 3)
for(i in 1:length(resamp.glmm.1995.lte4)){
    fe.dat.glmm.1995.lte4[i, ] <- c(resamp.glmm.1995.lte4[[i]]$'(Intercept)', resamp.glmm.1995.lte4[[i]]$DRMnum, resamp.glmm.1995.lte4[[i]]$lens)
}


#Write the fixed effects
write.table(fe.dat.glmm.1995.lte4, "../tmp/GLMM.1995.lte4.fixedeffs.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)



#1995+
#Truncated to 4 DRMs
#negative binomial GLM
#Parallelized to run  on 4 clusters
##############################

#Set the number of iterations
iters<- 1000
#Set up the cluster
cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()

#newdata3 will encode the data.frame that has all the values we want to predict on
#It is set up so that all sequences have the same length and different numbers of DRMs and makes the predictions for NNRTIs (rows 1-100), NRTIs (rows 101-200), PIs (rows 201-300) and PI/rs (rows 301-400)
newdata3 <- data.frame(
    DRMnum = rep(seq(from = 0, to =4, length.out = 100), 4),
    #lens and adjYear will be fixed
    lens = rep(800, 400),
    adjYear = rep(2005, 400),
    NNRTI.effect = c(seq(from = 0, to =4, length.out = 100),  rep(0, 300)),
    NRTI.effect = c(rep(0, 100), seq(from = 0, to =4, length.out = 100),  rep(0, 200)),
    PI.effect = c(rep(0, 200), seq(from = 0, to =4, length.out = 100),  rep(0, 100)),
    PIr.effect = c(rep(0, 300), seq(from = 0, to =4, length.out = 100))
)

predicts <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    #downsample the number of ambiguous reads
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
    #store them
    dat$downsamp.amb <- downsamp.amb
    #Fit the GLM
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[intersect(which(dat$IsolateYear >= 1995), which(dat$DRMnum <= 4)),])
    #Use the GLM fit to predict on the matrix we made above
    predict(m.allfacts, newdata3, type = "response")
}

#Exit 4-core mode
print(Sys.time()-strt)
stopCluster(cl)


#Store this into NNRTIs
NNRTIs <- matrix(data = NA, ncol = 100, nrow = iters)
NRTIs <- matrix(data = NA, ncol = 100, nrow = iters)
PIrs <- matrix(data = NA, ncol = 100, nrow = iters)
PIs <- matrix(data = NA, ncol = 100, nrow = iters)
for(i in 1:iters){
    NNRTIs[i,] <- predicts[[i]][1:100]
    NRTIs[i,] <- predicts[[i]][101:200]
    PIrs[i,] <- predicts[[i]][201:300]
    PIs[i,] <- predicts[[i]][301:400]
}

save(NNRTIs, file = "../tmp/NNRTIs.1995.lte4")
save(NRTIs, file = "../tmp/NRTIs.1995.lte4")
save(PIrs, file = "../tmp/PIrs.1995.lte4")
save(PIs, file = "../tmp/PIs.1995.lte4")


#we would like to do this again in order to just record the slopes (as opposed to the fits) so we can test for signficance. This is almost entirely the same code as above.
#GLM, slope differences

#parallelization set up
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#Run the model iters number of times
coefs.glm.1995.lte4 <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[intersect(which(dat$IsolateYear >= 1995), which(dat$DRMnum <= 4)),])
    coef(m.allfacts)
}

#Take the coefficients from the GLM and format them so that we can perform a t-test
coefs.glm.1995.lte4.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.1995.lte4.for.ttest[i,] <- coefs.glm.1995.lte4[[i]]
}
colnames(coefs.glm.1995.lte4.for.ttest) <- names(coefs.glm.1995.lte4[[1]])

#parallelization end
print(Sys.time()-strt)
stopCluster(cl)

write.table(coefs.glm.1995.lte4.for.ttest, "../tmp/GLM.1995.lte4.fixedeffs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

colnames(coefs.glm.1995.lte4.for.ttest)
coefs.glm.1995.lte4.for.ttest[,'NNRTI.effect']
coefs.glm.1995.lte4.for.ttest[,'NRTI.effect']


#CORR2

#0 DRMs:
NoDRMs.div <- coef[,'(Intercept)'] + coef[,'lens'] * 800 

#1NNRTI
OneDRM.NNRTI.div <- coef[,'(Intercept)'] + coef[,'lens'] * 800  + coef[,'NNRTI.effect']
#OneDRM.NNRTI.div <- coef[,'(Intercept)']  + coef[,'NNRTI.effect']

mean((NoDRMs.div - OneDRM.NNRTI.div)/NoDRMs.div)
quantile((NoDRMs.div - OneDRM.NNRTI.div)/NoDRMs.div, c(.025, .975))

#1NRTI
OneDRM.NRTI.div <- coef[,'(Intercept)'] + coef[,'lens'] * 800  + coef[,'NRTI.effect']

mean((NoDRMs.div - OneDRM.NRTI.div)/NoDRMs.div)
quantile((NoDRMs.div - OneDRM.NRTI.div)/NoDRMs.div, c(.025, .975))

#1PI
OneDRM.PI.div <- coef[,'(Intercept)'] + coef[,'lens'] * 800  + coef[,'PI.effect']

mean((NoDRMs.div - OneDRM.PI.div)/NoDRMs.div)
quantile((NoDRMs.div - OneDRM.PI.div)/NoDRMs.div, c(.025, .975))


#1PIr
OneDRM.PIr.div <- coef[,'(Intercept)'] + coef[,'lens'] * 800  + coef[,'PIr.effect']

mean((NoDRMs.div - OneDRM.PIr.div)/NoDRMs.div)
quantile((NoDRMs.div - OneDRM.PIr.div)/NoDRMs.div, c(.025, .975))




##############################
#1995+
#No truncation
##############################


#1995+
#No truncation
#GLMM, again, the GLMM is very slow to run. 
##############################

#Number of iterations
iters<- 20

#Set up parallelization
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#Isolate treatments that are abundant
com.treat.inds <- c()
for(i in refs.filt){
    com.treat.inds <- c(com.treat.inds, which(dat$Regimen == i))
}
tmpdat <- dat[com.treat.inds, ]
tmpdat$Regimen <- factor(tmpdat$Regimen)

#Fit the GLMMs and record
resamp.glmm.1995.notrunc <-foreach(icount(iters), .packages = 'glmmADMB') %dopar% {
    tmpdat$downsampamb <- rpois(length(tmpdat$ambnum), tmpdat$ambnum* down.sample.1995(tmpdat$IsolateYear))
    sto.withzero <- glmmadmb(downsampamb ~ 1 + DRMnum + lens + (1 + DRMnum|Regimen), data = tmpdat[(tmpdat$IsolateYear >= 1995),], family = "nbinom", zeroInflation = FALSE)
    c(coef(sto.withzero), sto.withzero[[23]])
}
print(Sys.time()-strt)
stopCluster(cl)


numTreats <- nrow(resamp.glmm.1995.notrunc[[1]]$Regimen)

#Reformat the data
resamp.dat.glmm.1995.notrunc <- matrix(data = NA, nrow = length(resamp.glmm.1995.notrunc), ncol = numTreats)
for(i in 1:length(resamp.glmm.1995.notrunc)){
    resamp.dat.glmm.1995.notrunc[i, ] <-  (resamp.glmm.1995.notrunc[[i]]$Regimen)[,2]
}

plotnames <- rownames(resamp.glmm.1995.notrunc[[1]]$Regimen)
colnames(resamp.dat.glmm.1995.notrunc) <- plotnames


#Store the random effects coefficients (so it can be graphed without having to rerun everything)
write.table(resamp.dat.glmm.1995.notrunc, "../tmp/GLMM.1995.notrunc.randeffs.txt", row.names = FALSE, col.names  = TRUE , quote = FALSE)

#Set up a matrix that will keep track of all the fixed effects
fe.dat.glmm.1995.notrunc <- matrix(data = NA, nrow = length(resamp.glmm.1995.notrunc), ncol = 3)
for(i in 1:length(resamp.glmm.1995.notrunc)){
    fe.dat.glmm.1995.notrunc[i, ] <- c(resamp.glmm.1995.notrunc[[i]]$'(Intercept)', resamp.glmm.1995.notrunc[[i]]$DRMnum, resamp.glmm.1995.notrunc[[i]]$lens)
}

write.table(fe.dat.glmm.1995.notrunc, "../tmp/GLMM.1995.notrunc.fixedeffs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


#1995+
#No truncation
#GLM
#####################################                                        

#Number of iterations
iters<- 1000

#Set up parallelization
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#newdata3 will encode the data.frame that has all the values we want to predict on
#It is set up so that all sequences have the same length and different numbers of DRMs and makes the predictions for NNRTIs (rows 1-100), NRTIs (rows 101-200), PIs (rows 201-300) and PI/rs (rows 301-400)
newdata3 <- data.frame(
    DRMnum = rep(seq(from = 0, to =4, length.out = 100), 4),
    #lens and adjYear will be fixed
    lens = rep(800, 400),
    adjYear = rep(2005, 400),
    NNRTI.effect = c(seq(from = 0, to =4, length.out = 100),  rep(0, 300)),
    NRTI.effect = c(rep(0, 100), seq(from = 0, to =4, length.out = 100),  rep(0, 200)),
    PI.effect = c(rep(0, 200), seq(from = 0, to =4, length.out = 100),  rep(0, 100)),
    PIr.effect = c(rep(0, 300), seq(from = 0, to =4, length.out = 100))
)

#Predictions (same as above)
predicts <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[which(dat$DRMnum <= 4),])
    predict(m.allfacts, newdata3, type = "response")
}

#end parallelization
print(Sys.time()-strt)
stopCluster(cl)

NNRTIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
NRTIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
PIrs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
PIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
for(i in 1:iters){
    NNRTIs.1995.notrunc[i,] <- predicts[[i]][1:100]
    NRTIs.1995.notrunc[i,] <- predicts[[i]][101:200]
    PIrs.1995.notrunc[i,] <- predicts[[i]][201:300]
    PIs.1995.notrunc[i,] <- predicts[[i]][301:400]
}

save(NNRTIs.1995.notrunc, file = "../tmp/NNRTIs.1995.notrunc")
save(NRTIs.1995.notrunc, file = "../tmp/NRTIs.1995.notrunc")
save(PIrs.1995.notrunc, file = "../tmp/PIrs.1995.notrunc")
save(PIs.1995.notrunc, file = "../tmp/PIs.1995.notrunc")



#Repeat this analysis to have t-tests
iters<-1000
cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()

coefs.glm.1995.notrunc <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[which(dat$IsolateYear >= 1995),])
    coef(m.allfacts)
}

print(Sys.time()-strt)
stopCluster(cl)

coefs.glm.1995.notrunc.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.1995.notrunc.for.ttest[i,] <- coefs.glm.1995.notrunc[[i]]
}
colnames(coefs.glm.1995.notrunc.for.ttest) <- names(coefs.glm.1995.notrunc[[1]])


write.table(coefs.glm.1995.notrunc.for.ttest, "../tmp/GLM.1995.notrunc.fixedeffs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)


######################################
#ALL DATA
#TRUNCATED TO 4 DRMs
#GLMM, very slow
######################################

#Number of iterations
iters<- 20 

#Set up parallelization
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#As above, find all indices of abundant treatments
com.treat.inds <- c()
for(i in refs.filt){
    com.treat.inds <- c(com.treat.inds, which(dat$Regimen == i))
}
tmpdat <- dat[com.treat.inds, ]
tmpdat$Regimen <- factor(tmpdat$Regimen)

#Fit GLMMs with resampling
resamp.glmm.all.lte4 <-foreach(icount(iters), .packages = 'glmmADMB') %dopar% {
    tmpdat$downsampamb <- rpois(length(tmpdat$ambnum), tmpdat$ambnum* down.sample(tmpdat$IsolateYear))
    sto.withzero <- glmmadmb(downsampamb ~ 1 + DRMnum + lens + (1 + DRMnum|Regimen), data = tmpdat[which(tmpdat$DRMnum <= 4),], family = "nbinom", zeroInflation = FALSE)
    c(coef(sto.withzero), sto.withzero[[23]])
}


#End parallelization
print(Sys.time()-strt)
stopCluster(cl)

numTreats <- nrow(resamp.glmm.all.lte4[[1]]$Regimen)

#Reformat data
resamp.dat.glmm.all.lte4 <- matrix(data = NA, nrow = length(resamp.glmm.all.lte4), ncol = numTreats)
for(i in 1:nrow(resamp.dat.glmm.all.lte4)){
    resamp.dat.glmm.all.lte4[i, ] <-  (resamp.glmm.all.lte4[[i]]$Regimen)[,2]
}


#Figure out names
plotnames <- rownames(resamp.glmm.all.lte4[[1]]$Regimen)
colnames(resamp.dat.glmm.all.lte4) <- plotnames

                                        #Store the random effects coefficients (so it can be graphed without having to rerun everything)
write.table(resamp.dat.glmm.all.lte4, "../tmp/GLMM.all.lte4.randeffs.txt", row.names = FALSE, col.names  = TRUE , quote = FALSE)

#Set up a matrix that will keep track of all the fixed effects
fe.dat.glmm.all.lte4 <- matrix(data = NA, nrow = length(resamp.glmm.all.lte4), ncol = 3)
for(i in 1:length(resamp.glmm.all.lte4)){
    fe.dat.glmm.all.lte4[i, ] <- c(resamp.glmm.all.lte4[[i]]$'(Intercept)', resamp.glmm.all.lte4[[i]]$DRMnum, resamp.glmm.all.lte4[[i]]$lens)
}

#Write the fixed effects
write.table(fe.dat.glmm.all.lte4, "../tmp/GLMM.all.lte4.fixedeffs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


#ALL DATA
#TRUNCATED TO 4 DRMs
#GLM

#Number of iterations
iters<- 1000

#Set up parallelization
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#newdata3 will encode the data.frame that has all the values we want to predict on
#It is set up so that all sequences have the same length and different numbers of DRMs and makes the predictions for NNRTIs (rows 1-100), NRTIs (rows 101-200), PIs (rows 201-300) and PI/rs (rows 301-400)
newdata3 <- data.frame(
    DRMnum = rep(seq(from = 0, to =4, length.out = 100), 4),
    #lens and adjYear will be fixed
    lens = rep(800, 400),
    adjYear = rep(2005, 400),
    NNRTI.effect = c(seq(from = 0, to =4, length.out = 100),  rep(0, 300)),
    NRTI.effect = c(rep(0, 100), seq(from = 0, to =4, length.out = 100),  rep(0, 200)),
    PI.effect = c(rep(0, 200), seq(from = 0, to =4, length.out = 100),  rep(0, 100)),
    PIr.effect = c(rep(0, 300), seq(from = 0, to =4, length.out = 100))
)

#Fits of the glm to newdata3
predicts.all.lte4 <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[which(dat$DRMnum <= 4),])
    predict(m.allfacts, newdata3, type = "response")
}

#End parallelization
print(Sys.time()-strt)
stopCluster(cl)

#NNRTI predictions
NNRTIs.all.lte4 <- matrix(data = NA, ncol = 100, nrow = iters)
NRTIs.all.lte4 <- matrix(data = NA, ncol = 100, nrow = iters)
PIrs.all.lte4 <- matrix(data = NA, ncol = 100, nrow = iters)
PIs.all.lte4 <- matrix(data = NA, ncol = 100, nrow = iters)
for(i in 1:iters){
    NNRTIs.all.lte4[i,] <- predicts.all.lte4[[i]][1:100]
    NRTIs.all.lte4[i,] <- predicts.all.lte4[[i]][101:200]
    PIrs.all.lte4[i,] <- predicts.all.lte4[[i]][201:300]
    PIs.all.lte4[i,] <- predicts.all.lte4[[i]][301:400]
}

save(NNRTIs.all.lte4, file = "../tmp/NNRTIs.all.lte4")
save(NRTIs.all.lte4, file = "../tmp/NRTIs.all.lte4")
save(PIrs.all.lte4, file = "../tmp/PIrs.all.lte4")
save(PIs.all.lte4, file = "../tmp/PIs.all.lte4")




#GLM, slope differences
#Models refit so that we can do statistical tests
iters<-1000

cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()

coefs.glm.all.lte4 <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[ which(dat$DRMnum <= 4),])
    coef(m.allfacts)
}

coefs.glm.all.lte4.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.all.lte4.for.ttest[i,] <- coefs.glm.all.lte4[[i]]
}
colnames(coefs.glm.all.lte4.for.ttest) <- names(coefs.glm.all.lte4[[1]])

print(Sys.time()-strt)
stopCluster(cl)

write.table(coefs.glm.all.lte4.for.ttest, "../tmp/GLM.all.lte4.fixedeffs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)




#Tests


wilcox.test(coefs.glm.1995.lte4.for.ttest[,3], coefs.glm.1995.lte4.for.ttest[,4], paired = TRUE)
wilcox.test(coefs.glm.1995.lte4.for.ttest[,5], coefs.glm.1995.lte4.for.ttest[,6], paired = TRUE)

wilcox.test(coefs.glm.1995.notrunc.for.ttest[,3], coefs.glm.1995.notrunc.for.ttest[,4], paired = TRUE)
wilcox.test(coefs.glm.1995.notrunc.for.ttest[,5], coefs.glm.1995.notrunc.for.ttest[,6], paired = TRUE)

wilcox.test(coefs.glm.all.lte4.for.ttest[,3], coefs.glm.all.lte4.for.ttest[,4], paired = TRUE)
wilcox.test(coefs.glm.all.lte4.for.ttest[,5], coefs.glm.all.lte4.for.ttest[,6], paired = TRUE)

sum(coefs.glm.1995.lte4.for.ttest[,3] < coefs.glm.1995.lte4.for.ttest[,4])
sum(coefs.glm.1995.lte4.for.ttest[,5]< coefs.glm.1995.lte4.for.ttest[,6])

sum(coefs.glm.1995.notrunc.for.ttest[,3]< coefs.glm.1995.notrunc.for.ttest[,4])
sum(coefs.glm.1995.notrunc.for.ttest[,5]< coefs.glm.1995.notrunc.for.ttest[,6])

sum(coefs.glm.all.lte4.for.ttest[,3]< coefs.glm.all.lte4.for.ttest[,4])
sum(coefs.glm.all.lte4.for.ttest[,5]< coefs.glm.all.lte4.for.ttest[,6])

#Table 1
fe.1995.lte4 <- read.table("../tmp/GLMM.1995.lte4.fixedeffs.txt")
signif(apply(fe.1995.lte4, 2, mean), 2)
paste("(", apply(signif(apply(fe.1995.lte4, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")

fe.1995.notrunc <- read.table("../tmp/GLMM.1995.notrunc.fixedeffs.txt")
signif(apply(fe.1995.notrunc, 2, mean), 2)
paste("(", apply(signif(apply(fe.1995.notrunc, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")

fe.all.lte4 <- read.table("../tmp/GLMM.all.lte4.fixedeffs.txt")
signif(apply(fe.all.lte4, 2, mean), 2)
paste("(", apply(signif(apply(fe.all.lte4, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")


#Table 2
#Note, the PIr and PI order is switched from the table in the paper, but the labels are correct
signif(apply(coefs.glm.1995.lte4.for.ttest, 2, mean), 2)
paste("(", apply(signif(apply(coefs.glm.1995.lte4.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")

signif(apply(coefs.glm.1995.notrunc.for.ttest, 2, mean), 2)
paste("(", apply(signif(apply(coefs.glm.1995.notrunc.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")

signif(apply(coefs.glm.all.lte4.for.ttest, 2, mean), 2)
paste("(", apply(signif(apply(coefs.glm.all.lte4.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")
