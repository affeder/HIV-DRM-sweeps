#Resamp file
require(vcd)
require(glmmADMB)
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")
require(parallel)
require(doParallel)


#####################
# Read in the data
#####################

source("read-in-data.r") 

drugnames <- drugnames[order(percentfail)]
percentfail <- percentfail[order(percentfail)]

###################
# Compute subsampling amounts
###################

years <- sort(unique(dat$IsolateYear))
zeroclass <- which(dat$DRMnum == 0)

year.zero.mean <- c() # to get the mean number of amb reads in each year in the zero DRM class
for(i in 1:length(years)){
    year.zero.mean[i] <- mean(dat[intersect(zeroclass,which(dat$IsolateYear == years[i])),'ambnum'])
}

Zeromean<-data.frame(years=years, zeromean = 0) # this creates a dataframe with the mean number of amb reads for each year in the zero class. 
for(i in 1:length(years)){
    Zeromean$zeromean[i]<-mean(dat$ambnum[which(dat$DRMnum == 0&dat$IsolateYear == years[i])])
}


#lm(dat[zeroclass, 'ambnum'] ~ dat[zeroclass, 'IsolateYear'])$coefficients

#Full (i.e., all data, not just 1995+)
fit2b <- lm(dat[zeroclass, 'ambnum'] ~ dat[zeroclass, 'IsolateYear'])$coefficients 
poly2b <- function(year){c(fit2b[1] + fit2b[2]*year )}
down.sample <- function(year){
    #toRet stores the values to return (toRet) 
    #this value is the ratio between the average number of amb reads in the 0 
    #category in 1989 and the inputted year. This value is used to determine
    #how much p-thinning must be done in the inputted year to achive 1989 levels of
    #ambiguous reads.
    toRet <- poly2b(1989)/poly2b(year) 
    names(toRet) <- ""
    c(toRet)
}
fit2b

#1995+
fit2b.1995 <- lm(dat[intersect(zeroclass, which(dat$IsolateYear >=1995)), 'ambnum'] ~ dat[intersect(zeroclass, which(dat$IsolateYear >=1995)), 'IsolateYear'])$coefficients
poly2b.1995 <- function(year){c(fit2b.1995[1] + fit2b.1995[2]*year )}
down.sample.1995 <- function(year){
    #This function works as above, but computes p-thinning with respect to 1995
    toRet <- poly2b.1995(1995)/poly2b.1995(year)
    names(toRet) <- ""
    c(toRet)
}
fit2b.1995

#print fit2b and fit2b.1995, as they're recorded in the text
coefsToPrint <- rbind(fit2b, fit2b.1995)
colnames(coefsToPrint) <- c("Intercept", "Slope (by year)")
rownames(coefsToPrint) <- c("Relative to 1989", "Relative to 1995")
write.table(coefsToPrint, "../figures/pthinning-coefs.txt", row.names = TRUE, col.names = TRUE)

#Setup plot F5-S3
pdf("../figures/F5-S3.pdf", width =9, height =5)
par(mar = c(4,4,1,1))
layout(matrix(1:2, nrow = 1))
# F5-S3A This plot shows the mean number of amb calls for each year in the zero class
plot(years, year.zero.mean, ylab = "Number of Ambiguous Calls", xlab = "Year", pch = 16, col = "grey", ylim = c(0, 10))
# add line for poly2b model
lines(years, poly2b(years))
points(years, poly2b(years), col = "black", pch = 16)
# add line for poly2b.1995 model
lines(seq(1995, 2013, by = 1), poly2b.1995(seq(1995, 2013, by = 1)), col = "red")
points(seq(1995, 2013, by = 1), poly2b.1995(seq(1995, 2013, by = 1)), col = "red", pch =16)
legend("bottomright", c("Mean of 0-DRM seqs.", "1989+ model fits", "1995+ model fits"), col = c("grey", "black", "red"), pch = c(16, 16, 16),  lty = c(0, 1, 1))
# F5-S3B This plot shows how much the data are downsampled depending on whether we start at 1989 or 1995
plot(years, 1/(poly2b(years)/min(poly2b(years))), xlab = "Year", ylab = "Subsample effect", pch = 3, cex = .5)
points(seq(1995, 2013, by = 1), 1/(poly2b(seq(1995, 2013, by = 1))/min(poly2b(seq(1995, 2013, by = 1)))), ylab = "Subsample effect", pch = 3, cex = .5, col = "red")
legend("topright", c("1989+ model", "1995+ model"), col = c("black", "red"),  pch = c( 3, 3))
dev.off()

###################
# Assess model fit
###################

# We found that the subsampled data visually fit a negative binomial distribution much more closely than a Poisson distribution, 
#which is often used for count data (See Figure 5—figure supplement 4 and Figure 5—figure supplement 5).

#These are the specific draws that will be plotted
downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum*down.sample.1995(dat$IsolateYear)) # creates a downsampled dataset, based on the 1995+ model. 
nbfit <- goodfit(downsamp.amb, "nbinomial") # creates a neg binomial dist to fit the dataset that was just made
poisfit <- goodfit(downsamp.amb, "poisson") # creates a poisson dist to fit the dataset that was just made

downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample(dat$IsolateYear))# creates a downsampled dataset, based on all data
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
# This plot F5-S4 shows how the downsampled data using all data fit a poisson or neg binom model. 
layout(matrix(1:2, nrow = 1))
par(mar = c(4,4,3,1))
fit.plot(poisfit.all, "Poisson Fits", "A", 48, 39)
fit.plot(nbfit.all, "Negative Binomial Fits", "B", 48, 48)
dev.off()

pdf("../figures/F5-S5.pdf", width =8, height =4.5)
# This plot F5-S5 shows how the downsampled data using data from 1995 fit a poisson or neg binom model. 
layout(matrix(1:2, nrow = 1))
par(mar = c(4,4,3,1))
fit.plot(poisfit, "Poisson Fits", "A", 68, 34.5)
fit.plot(nbfit, "Negative Binomial Fits", "B", 68, 44.5)
dev.off()

##############################
# We estimate the relationship between the number of DRMs and genetic diversity by fitting a 
#genelalized linear mixed model (GLMM) with a negative binomial error distribution for our 31 abundant treatments. 
#In this model, length, the number of DRMs and an intercept term are fit as fixed effects, 
#and the number of DRMs by treatment is fit as a random effect. 
#This allows us to assess the relationship between diversity and the 
#number of ambiguous reads separately for each treatment. 
#The models were fit using the glmmADMB package in R (Fournier et al., 2012).
#Subsampled number of ambiguous reads ~ DeltaDRM,all(numDRM) + Alpha_all + Gamma(Sequence Length) + 
#                                           (Alpha_t + DeltaDRM,t(numDRM)|Regimen)
#This model was fit to 1000 datasets created by p-thinning the number of ambiguous calls. 
#The overall effect of a DRM on diversity was fit by the DeltaDRM,all term, 
#but the effect of a DRM on diversity by treatment t is fit by the random effect term, DeltaDRM,t.



#The full effect of a DRM on diversity for a given treatment was called DeltaDRM and was computed by 
#combining the treatment-specific random effect and the overall fixed effect of the model (DeltaDRM = DeltaDRM,t + DeltaDRM,all). 
#Confidence intervals were generated by excluding the highest and lowest 2.5% of estimates of DDRM among the subsamples.

#We performed the above procedure three times: 
#    1. our main analysis used the p1995;i procedure for p-thinning sequences from year i and included sequences with 4 or fewer DRMs. 
#    2. We performed this same analysis using p1995;i-thinning and including all numbers of DRMs and 
#    3. using p1989;i-thinning and including only sequences with 4 or fewer DRMs.

##############################
###
#PREPARATION
####

com.treat.inds <- c() # create list of inds from common treatments
#Only run this for treatments that are abundant
for(i in refs.filt){ #refs.filt is list of common treatments
    com.treat.inds <- c(com.treat.inds, which(dat$Regimen == i))
}
tmpdat <- dat[com.treat.inds, ]
tmpdat$Regimen <- factor(tmpdat$Regimen)#The regimens need to be converted to factors


##############################
#    1. our main analysis used the p1995;i procedure for p-thinning sequences from year i and included sequences with 4 or fewer DRMs. 
#1995+
#Truncated to 4 DRMs
##############################
#GLMM, note: GLMM code is quite slow to run
##############################

#Note, in the paper, we run 1000 iterations, but it's quite slow, so we've updated it so that only 20 iterations (parallelized to four cores)
iters<- 20

#Set up code that can be run on 4 cores
cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()

#to be run 'iters' number of times # using foreach function bc it allows us to use 4 cores. See https://beckmw.wordpress.com/tag/foreach/
resamp.glmm.1995.lte4 <-foreach(icount(iters), .packages = 'glmmADMB') %dopar% {
    #subsample the number of ambiguous reads
    tmpdat$downsampamb <- rpois(length(tmpdat$ambnum), tmpdat$ambnum* down.sample.1995(tmpdat$IsolateYear))
                                        #Fit the GLMM (lens = length seq)
    GLMM <- glmmadmb(downsampamb ~ 1 + DRMnum + lens + (1 + DRMnum|Regimen), data = tmpdat[intersect(which(tmpdat$IsolateYear >= 1995), which(tmpdat$DRMnum <= 4)),], family = "nbinom", zeroInflation = FALSE)
    #return the GLMM information (more processing will be necessary)
    c(coef(GLMM), GLMM[[23]]) #things to return (and save into resamp.glmm.1995.lte4)
    # coef(sto.withzero) returns the fixed effects coefficients 
    #sto.withzero[[23]] returns the the random effects coefficients  (i.e., the drug effects)
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

#Store the random effects coefficients (so it can be graphed without having to rerun everything)
write.table(resamp.dat.glmm.1995.lte4, "../tmp/GLMM.1995.lte4.randeffs.txt", row.names = FALSE, col.names  = TRUE , quote = FALSE)

#Set up a matrix that will keep track of all the fixed effects
fe.dat.glmm.1995.lte4 <- matrix(data = NA, nrow = length(resamp.glmm.1995.lte4), ncol = 3)
for(i in 1:length(resamp.glmm.1995.lte4)){
    fe.dat.glmm.1995.lte4[i, ] <- c(resamp.glmm.1995.lte4[[i]]$'(Intercept)', resamp.glmm.1995.lte4[[i]]$DRMnum, resamp.glmm.1995.lte4[[i]]$lens)
}
#Write the fixed effects
write.table(fe.dat.glmm.1995.lte4, "../tmp/GLMM.1995.lte4.fixedeffs.txt",  row.names = FALSE, col.names = FALSE, quote = FALSE)

#######################
#Parametric approach
#1995+
#Truncated to 4 DRMs
#negative binomial GLM 

#To compare how the effect of DRMs on genetic diversity varied between two groups of treatments, 
#we fit generalized linear models (GLMs) with a negative binomial error distribution 
#using the package pcsl (Jackman, 2015) including all sequences belonging to the 
#31 treatments that passed our threshold criteria (see above).
#These models were parametrized to fit separate slopes for the 
#four different types of treatments (1, 2 or 3 NRTIs, 2NRTIs + NNRTI, 2NRTIs + PI, 2NRTIs + PI/r). 

#Parallelized to run  on 4 clusters
##############################

#Set the number of iterations
iters<- 1000
#Set up code that can be run on 4 cores 
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
    c(predict(m.allfacts, newdata3, type = "response"), coef(m.allfacts)) # save the predictions to "predicts" along with the coefficients that will be used in a t-test
}

#Exit 4-core mode
print(Sys.time()-strt)
stopCluster(cl)


#  Store data and write to files
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

#Take the coefficients from the GLM and format them so that we can perform a t-test
coefs.glm.1995.lte4.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.1995.lte4.for.ttest[i,] <- predicts[[i]][401:406]
}
colnames(coefs.glm.1995.lte4.for.ttest) <- names(predicts[[1]][401:406])

write.table(coefs.glm.1995.lte4.for.ttest, "../tmp/GLM.1995.lte4.fixedeffs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)



#this code compute some numbers presented in the text as the effect of having an NRTI/NNRTI/PI/PI/r
coef <- coefs.glm.1995.lte4.for.ttest
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
#    2. We performed this same analysis using p1995;i-thinning and including all numbers of DRMs and 
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

#Predictions (same as above)
predicts.1995.notrunc <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample.1995(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[which(dat$DRMnum <= 4),])
    c(predict(m.allfacts, newdata3, type = "response"), coef(m.allfacts))
}

#end parallelization
print(Sys.time()-strt)
stopCluster(cl)


NNRTIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
NRTIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
PIrs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
PIs.1995.notrunc <- matrix(data = NA, ncol = 100, nrow = iters)
for(i in 1:iters){
    NNRTIs.1995.notrunc[i,] <- predicts.1995.notrunc[[i]][1:100]
    NRTIs.1995.notrunc[i,] <- predicts.1995.notrunc[[i]][101:200]
    PIrs.1995.notrunc[i,] <- predicts.1995.notrunc[[i]][201:300]
    PIs.1995.notrunc[i,] <- predicts.1995.notrunc[[i]][301:400]
}

save(NNRTIs.1995.notrunc, file = "../tmp/NNRTIs.1995.notrunc")
save(NRTIs.1995.notrunc, file = "../tmp/NRTIs.1995.notrunc")
save(PIrs.1995.notrunc, file = "../tmp/PIrs.1995.notrunc")
save(PIs.1995.notrunc, file = "../tmp/PIs.1995.notrunc")

coefs.glm.1995.notrunc.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.1995.notrunc.for.ttest[i,] <- predicts.1995.notrunc[[i]][401:406]
}
colnames(coefs.glm.1995.notrunc.for.ttest) <- names(predicts.1995.notrunc[[1]][401:406])

write.table(coefs.glm.1995.notrunc.for.ttest, "../tmp/GLM.1995.notrunc.fixedeffs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)



######################################
#    3. using p1989;i-thinning and including only sequences with 4 or fewer DRMs.
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

####################
#ALL DATA
#TRUNCATED TO 4 DRMs
#GLM
####################

#Number of iterations
iters<- 1000

#Set up parallelization
cl<-makeCluster(4)
registerDoParallel(cl)
strt<-Sys.time()

#Fits of the glm to newdata3
predicts.all.lte4 <-foreach(icount(iters), .packages = 'MASS') %dopar% {
    downsamp.amb <- rpois(length(dat$ambnum), dat$ambnum* down.sample(dat$IsolateYear))
    dat$downsamp.amb <- downsamp.amb
    m.allfacts <- glm.nb(downsamp.amb ~  lens + NNRTI.effect + NRTI.effect + PIr.effect + PI.effect , data = dat[which(dat$DRMnum <= 4),])
    c(predict(m.allfacts, newdata3, type = "response"), coef(m.allfacts))
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

coefs.glm.all.lte4.for.ttest <- matrix(data = NA, ncol = 6, nrow = iters)
for(i in 1:iters){
    coefs.glm.all.lte4.for.ttest[i,] <- predicts.all.lte4[[i]][401:406]
}
colnames(coefs.glm.all.lte4.for.ttest) <- names(predicts.all.lte4[[1]][401:406])


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
fe.1995.notrunc <- read.table("../tmp/GLMM.1995.notrunc.fixedeffs.txt")
fe.all.lte4 <- read.table("../tmp/GLMM.all.lte4.fixedeffs.txt")

tab1 <- rbind( signif(apply(fe.1995.lte4, 2, mean), 2), 
      paste("(", apply(signif(apply(fe.1995.lte4, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = ""),
      signif(apply(fe.1995.notrunc, 2, mean), 2), 
      paste("(", apply(signif(apply(fe.1995.notrunc, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = ""),
      signif(apply(fe.all.lte4, 2, mean), 2), 
paste("(", apply(signif(apply(fe.all.lte4, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = "")
      )

require(xtable)

tab1 <- cbind(c("1995+, $leq$ 4 DRMs", "", "1995+, all DRMs", "", "1989+, $leq$ 4 DRMs", ""), tab1)
colnames(tab1) <- c("  ", "$alpha_{all}$  (Intercept)", "$Delta$ (Number of DRMs)", "$gamma$ (Length)")
write(print(xtable(tab1), include.rownames = FALSE), file = "../figures/tab1.txt")

#re-read in data for table 2
coefs.glm.1995.lte4.for.ttest <- read.table("../tmp/GLM.1995.lte4.fixedeffs.txt",header = TRUE)
coefs.glm.1995.notrunc.for.ttest <- read.table("../tmp/GLM.1995.notrunc.fixedeffs.txt",header = TRUE)
coefs.glm.all.lte4.for.ttest <- read.table("../tmp/GLM.all.lte4.fixedeffs.txt",header = TRUE)

#Table 2
#Note, the PIr and PI order is switched from the table in the paper, but the labels are correct
tab2 <- rbind( signif(apply(coefs.glm.1995.lte4.for.ttest, 2, mean), 2), 
      paste("(", apply(signif(apply(coefs.glm.1995.lte4.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = ""), 
      signif(apply(coefs.glm.1995.notrunc.for.ttest, 2, mean), 2), 
      paste("(", apply(signif(apply(coefs.glm.1995.notrunc.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = ""), 
      signif(apply(coefs.glm.all.lte4.for.ttest, 2, mean), 2), 
paste("(", apply(signif(apply(coefs.glm.all.lte4.for.ttest, 2, quantile, c(.025, .975)), 2), 2, paste, collapse = ","), ")", sep = ""))

tab2 <- cbind(c("1995+, $leq$ 4 DRMs", "", "1995+, all DRMs", "", "1989+, $leq$ 4 DRMs", ""), tab2)
colnames(tab2) <- c("  ", "$alpha_{all}$  (Intercept)",  "$gamma$ (Length)", "$Delta_{2NRTI+NNRTI}$", "$Delta_{1,2,3NRTI}$", "$Delta_{2NRTI+PI/r}$", "$Delta_{2NRTI+PI}$")

write(print(xtable(tab2), include.rownames = FALSE), file = "../figures/tab2.txt")


