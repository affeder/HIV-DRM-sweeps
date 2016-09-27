#read.in.data.r
#Read in the majority of the data

PRseqcaps <- read.table("../tmp/PR_nuc.txt", stringsAsFactors = FALSE)
RTseqcaps <- read.table("../tmp/RT_nuc.txt", stringsAsFactors = FALSE)
praaalign <- read.table("../tmp/PR_aa.txt", stringsAsFactors = FALSE)
aaalign <- read.table("../tmp/RT_aa.txt", stringsAsFactors = FALSE)

load("../tmp/hardDRMlist.txt")
load("../tmp/softDRMlist.txt")
load("../tmp/DRMlist.txt")

dat <- read.table("../tmp/nclonal.dat.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
oldnames <- names(dat)
oldnames[oldnames == "newRegimen"] <- "Regimen"
names(dat) <- oldnames
names(dat)


ncaaalign <- read.table("../tmp/ncaaalign.txt", stringsAsFactors = FALSE, header = TRUE)
ncpraaalign <- read.table("../tmp/ncpraaalign.txt", stringsAsFactors = FALSE, header = TRUE)
ncRTseqcaps <- read.table("../tmp/ncRTseqcaps.txt", stringsAsFactors = FALSE, header = TRUE)
ncPRseqcaps <- read.table("../tmp/ncPRseqcaps.txt", stringsAsFactors = FALSE, header = TRUE)


NNRTI <- "EFV|RPV|ETR|DLV|NVP|NNRTI"
NRTI <- "AZT|3TC|FTC|DDI|ABC|D4T|TDF|ZDV|DDC|^NRTI"
PI <- "TPV|IDV|SQV|LPV|FPV|DRV|ATV|NFV|APV|PI"


lenfunc <- function(x){
    c(sum(regexpr('\\.|N', x) < 1))
}

rtlen <-  apply(ncRTseqcaps, 1, lenfunc)
lens <- rtlen 

ambfunc <- function(x){
    c(sum(regexpr('\\.|N|A|T|C|G', x) < 1))
}
ambnum <- apply(ncRTseqcaps, 1, ambfunc)

dat$ambnum <- ambnum
dat$lens <- lens

#table(dat[dat$Regimen == "ATV+FTC+RTV+TDF", 'DRMnum'])
#table(dat[dat$Regimen == "3TC+ABC+ATV+RTV", 'DRMnum'])
table(dat[dat$Regimen == "3TC+ABC+ATV+RTV", ]$DRMnum)
table(dat[dat$Regimen == "ATV+FTC+RTV+TDF", ]$DRMnum)

dat[4964,]$DRMnum
DRMlist[[4964]]
sort(table(unlist(DRMlist)))

#Start with all treatments
alltreats <- unique(dat$Regimen)

#For each treatment, how many sequences do we have in each DRM category? (numDRM.bt)
numDRM.bt <- matrix(data = NA, ncol = 13, nrow = length(alltreats))
for(i in 1:length(alltreats)){
    subs <- dat[dat$Regimen == alltreats[i], ]
    for(j in 0:12){
        numDRM.bt[i, j+1] <- nrow(subs[subs$DRMnum == j,])
    }
}




#Our list of filtered treatments are all those in which we have at least three categories with at least 5 sequences
refs.filt <- alltreats[which(apply(numDRM.bt, 1, function(x){ c(sum(x >= 3) >=3 )}))]

numDRM.bt[which(alltreats == "3TC+ABC+ATV"),]
alltreats == "3TC+ABC+ATV"

#Determine NRTI/NNRTI/etc and store it in reg.class
reg.class <- c()
for(i in refs.filt){
    components <- strsplit(i, split = "\\+")[[1]]
    boosted = 0
    nrti = 0
    nnrti = 0
    pi = 0
    other = 0
    for(j in components){
        if(regexpr(PI, j) > 0){pi = pi+1}
        if(regexpr(NRTI, j) > 0){nrti = nrti+1}
        if(regexpr(NNRTI, j) > 0){nnrti = nnrti+1}
        if(j == "RTV"){ boosted = 1 }
        if(j == "LPV"){ boosted = 1 }
        if((regexpr(PI, j) < 0) & (regexpr(NRTI, j) < 0) & (regexpr(NNRTI, j) < 0) & (j != "RTV")){ other = 1 }
    }
    reg.class <- c(reg.class, paste(nrti,",",nnrti,",", pi, ",", boosted, ",", other, sep = ""))
}

#We want to just recode this into a number that's easy to filter by
coder <- rep(0, length(reg.class))
coder[reg.class == "1,0,0,0,0"] <- 1 #1 NRTI
coder[reg.class == "2,0,0,0,0"] <- 2 #2 NRTI
coder[reg.class == "3,0,0,0,0"] <- 3 #3 NRTI
coder[reg.class == "2,1,0,0,0"] <- 4 #2 NRTI, 1 NNRTI
coder[reg.class == "2,0,1,0,0"] <- 5 #2 NRTI, 1 PI
coder[reg.class == "2,0,1,1,0"] <- 6 #2 NRTI, 1 PI/r

toSort <- matrix(data = NA, ncol = 2, nrow = 31)
for( i in refs.filt[coder > 0 ]){
    toSort[min(which(is.na(toSort[,1]))),] <-  c(i, as.numeric(length(which(dat$Regimen == i))))
}


cbind(refs.filt, coder)
#Ok, I also want to know how many sequences from each treatment category there are for each of the drugs in refs.filt
abundantinf <- matrix(data = NA, nrow = length(refs.filt), ncol <- 5)
for(i in 1:length(refs.filt)){

   giventreat <- dat[dat$Regimen == refs.filt[i], ]
   givenapi <- (regexpr(PI, refs.filt[i]) > 0)
   if(givenapi){
       giventreat <- giventreat[giventreat$PRlen > 0, ]
   }
   abundantinf[i,] <- (c(refs.filt[i], givenapi, nrow(giventreat),sum(giventreat$DRMnum <= 4), sum(giventreat$DRMnum > 4)))
   
}


treatments <- refs.filt[coder > 0 & coder < 5]
NNRTI.treats <- refs.filt[coder == 4]
NRTI.treats <- setdiff(treatments, NNRTI.treats)
PI.treats <- refs.filt[coder == 5]
PIr.treats <- refs.filt[coder == 6]

NNRTItreatinds <- c()
for(i in NNRTI.treats){
    NNRTItreatinds <- c(NNRTItreatinds, which(dat$Regimen == i))
}
NRTItreatinds <- c()
for(i in NRTI.treats){
    NRTItreatinds <- c(NRTItreatinds, which(dat$Regimen == i))
}
PItreatinds <- c()
for(i in PI.treats){
    PItreatinds <- c(PItreatinds, which(dat$Regimen == i))
}
PIrtreatinds <- c()
for(i in PIr.treats){
    PIrtreatinds <- c(PIrtreatinds, which(dat$Regimen == i))
}

PI.effect <- rep(0, nrow(dat))
PI.effect[PItreatinds] <- dat$DRMnum[PItreatinds]

PIr.effect <- rep(0, nrow(dat))
PIr.effect[PIrtreatinds] <- dat$DRMnum[PIrtreatinds]

NNRTI.effect <- rep(0, nrow(dat))
NNRTI.effect[NNRTItreatinds] <- dat$DRMnum[NNRTItreatinds]

NRTI.effect <- rep(0, nrow(dat))
NRTI.effect[NRTItreatinds] <- dat$DRMnum[NRTItreatinds]

dat$NNRTI.effect <- NNRTI.effect
dat$NRTI.effect <- NRTI.effect
dat$PI.effect <- PI.effect
dat$PIr.effect <- PIr.effect




##############################################
# Read in the efficacy information
##############################################

#Concert the efficacy information to a usable format
tmp <- read.table("../dat/effic3.csv", sep = ",", fill = T, header = T, stringsAsFactors = F)
names(tmp) <- c("Treatment.Name", "Study.name", "succprop", "N", "Weeks", "VL.limit")

allnames <- unique(tmp$Treatment.Name)
topnames <- names(sort(table(dat$Regimen), decreasing = T)[1:30])
drugnames <- intersect(topnames, allnames)

tmp2 <- tmp[tmp$Weeks > 47 & tmp$Weeks < 53, ]
treat.exclude <- c()
for(i in which(tmp2$VL.limit > 50)){
    if(nrow(tmp2[tmp2$Treatment.Name == tmp2[i,]$Treatment.Name,]) > 1){
        treat.exclude <- c(treat.exclude,i)
    }
}
tmp <- tmp2[-treat.exclude,]

percentfail <- c()
for(i in 1:length(drugnames)){
    drugrel <- tmp[tmp$Treatment.Name == drugnames[i],]
    print(drugnames[i])
    print(nrow(drugrel))
    percentfail[i] <- sum(as.numeric(drugrel$succprop), na.rm = T)/sum(as.numeric(drugrel$N), na.rm = T)
}
percentfail <- percentfail*100

#Ok, let's just take the non-NAs
drugnames <- drugnames[!is.na(percentfail)]
percentfail <- percentfail[!is.na(percentfail)]


relInds <- c()
for(i in 1:length(refs.filt)){
    relInds <- c(relInds, which(dat$Regimen == refs.filt[i]))
}


