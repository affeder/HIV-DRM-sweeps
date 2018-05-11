library(seqinr)
library(lmtest)

# The goal of this script is to read in the initial data file, 
# Align it in a sensible way
# Print the newly formatted data out

##############################################
#Table of contents:
#t1.1 Functions to initially align the data
#t1.2 Align the data and print out aligned copies to /tmp
#t1.3 Split the data into the clonal and non-clonal dataset
#t1.4 Call drug resistance mutations in each sequence
#Note: you can SEARCH for the above titles in this document, but you will likely have dependency problems if you don't run it in order
#This is mostly just organizational
##############################################


##############################################
# t1.1 Functions to initially align the data
##############################################

#Split all characters 
conv <- function(x){
    return(strsplit(as.vector(x), split = "")[[1]])
}

#Standardize all regimen names so that they appear as sorted alphabetically separated by +s.
#This would mean, for example, both "AZT+3TC" and "3TC+AZT" appear as "3TC+AZT"
standardizeRegimenName <- function(x){
    sort(strsplit(as.vector(x), split = '\\+')[[1]])
}

#Remove whitespace
#http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#Use the extra information about first and last position in the nucleotide sequence to align the sequences in a way that makes them
#easier to handle.
alignRT<- function(x){

#For reasons that I do not entirely understand relating to apply, I need to index these as numeric indices instead of names 
    #index 11 == RTFirstAA
    #index 12 == RTLastAA 
    first <- (as.numeric(x[11]) - 1) *3 + 1
    last <- (as.numeric(x[12])) * 3 

    #13 = RTseq
    aligned <- paste(x[13])
    longest <- max(dat$RTLastAA)*3

    while(first > 1){
        aligned <- paste(".", aligned, sep = "")
        first <- first - 1
    }

    while(last < longest){
        aligned <- paste(aligned, ".", sep = "")
        last <- last + 1
    }

    return(aligned)
}

#Align PR does the same thing as above, but it further accounts for the fact that not all sequences had nucleotide information
#If it can't align, returns NA (note, all actual alignment taken from HIV DR DB
alignPR<- function(x){

    #Index 14 = PRFirstAA
    #Index 15 = PRLastAA
    #Index 16 = PRSeq

    if(x[14] == "NULL" | x[15] == "NULL"){
        return(NA)
    }
    first <- (as.numeric(x[14]) - 1) *3 + 1
    last <- as.numeric(x[15]) * 3 

    aligned <- paste(x[16])
    longest <- max(as.numeric(dat$PRLastAA), na.rm = T)*3

    while(first > 1){
        aligned <- paste(".", aligned, sep = "")
        first <- first - 1
    }

    while(last < longest){
        aligned <- paste(aligned, ".", sep = "")
        last <- last + 1
    }

    return(aligned)
}



#Options for nucleotides
nucOpts <- c("A","C","T","G","B","D","H","K","M","N","R","S","V","W","Y","a","c","g","k","m","r","s","t","w","y","~")

#The next three functions compute metric entropy - for the final analysis, I believe we ended up using pi instead, so I'm not sure that these are used.
#nucleotide level metric entropy
ent <- function(y){
    x <- table(subset(y, y != "." & y != "N"))
    c(entropy(x)/length(x))
}
#amino acid level metric entropy (excluding Xs)
entAA <- function(y){
    x <- table(subset(y, y != "X"))
    c(entropy(x)/length(x))
}
#Metric entropy  
enttab <- function(x){
    c(entropy(x)/length(x))
}

#This function translates nucleotides into 'ambiguous' amino acids. 
#i.e., AAW would be AAA or AAT which would be K or N (return "K/N")
translateAmbigNA<- function(x){

    helper <- function(y){
        if(is.na(y)){return(c(NA))}
        if(y == "A"){return(c("A"))}
        if(y == "T"){return(c("T"))}
        if(y == "C"){return(c("C"))}
        if(y == "G"){return(c("G"))}
        if(y == "W"){return(c("A", "T"))}
        if(y == "S"){return(c("C", "G"))}
        if(y == "M"){return(c("A", "C"))}
        if(y == "K"){return(c("G", "T"))}
        if(y == "R"){return(c("A", "G"))}
        if(y == "Y"){return(c("C", "T"))}
        if(y == "B"){return(c("C", "G", "T"))}
        if(y == "D"){return(c("A", "G", "T"))}
        if(y == "H"){return(c("A", "C", "T"))}
        if(y == "V"){return(c("A", "C", "G"))}
        if(y == "."){return(c("."))}
        if(y == "~"){return(c("~"))}
        if(y == "N"){return(c("N"))}

    }

    #Compute the unambiguous versions of each nucleotide (i.e., W into A/T)
    a <- helper(x[1])
    b <- helper(x[2])
    c <- helper(x[3])

    #If none of them are defined, toss it out
    if(is.na(a[1]) | is.na(b[1]) | is.na(c[1])){ return(NA) }
    #PSP this gives a warning when a, b or c are ambiguous. I think it is not a big problem,
    # but it is a bit annoying that the warnings pop up. 

    #Otherwise, translate all combinations
    toRet <- c()
    for(i in a){
        for(j in b){
            for(k in c){
                toRet <- c(toRet, translate(c(i,j,k)))
            }
        }
    }
    return(paste(sort(toRet), collapse = "/"))
}


#For an amino acid number, return the nucleotide positions encoding that amino acid
pos <- function(x){ c(x*3 -2, x*3 -1, x*3 ) }

polycheck <- function(x){
    if(length(unique(x[x != "."]) > 1)){ return(1) };
    return(false);
}

#How many ambiguous reads are there (not Ns or .s)
numAmbig <- function(x){
    nondot <- (x[x != "."])
    return(length(grep('[atcgATCGN]', nondot, perl = T, invert = T)))#/size)
}

#Same thing, except only works on capital letters
numAmbigCaps <- function(x){
    nondot <- (x[x != "."])
    return(length(grep('[ATCGN]', nondot, perl = T, invert = T)))#/size)
}


#Compute the number of objects in an array, except if the vector just has 0
#i.e., listlen(c("K103N", "M184V", "L75V")) should be 3
# listlen(c(0)) should be 0
#c("K103N")) should be 1 
listlen <- function(x){
    if(length(x) == 0){
        return(0)
    }
    if(length(x) == 1 && x[1] == 0){ # PSP I was getting warnings about "the condition has length > 1 and only the first element will be used"
        #so I changed "x == 0" into "length(x) == 1 && x[1] == 0"
        return(0)
    }
    return(length(x))
}

#as numeric function that plays nicely with factors
asNumeric <- function(x) as.numeric(as.character(x))




#####################################################
# t1.2 Align the data and print out aligned copies to /tmp
#
# Read in the data file, align it, and print out temporary files. 
#####################################################

#Read in the initial dataset to process
#This still needs to be somewhat parsed and aligned
dat <- read.table("../dat/dataset.01.24.txt", skip = 13, sep= "\t", header = T, stringsAsFactors = FALSE)

#PSP to make all this run fast to check
#dat<-dat[sample(1:length(dat[,1]), 2000),]

#We'll need to create alignments. We're given the start and end positions of each of these sequences, but let's arrange them in matrices with missing data
RTalign <- unlist(apply(dat, 1, alignRT))
PRalign <- unlist(apply(dat, 1, alignPR))
#Note: warnings for PRalign are bad form, but expected.

#old way
#a <- dat
#PRalign <- c()
#for(i in 1:nrow(a)){
#    if(i %% 1000 == 0){ print(i) }
#    PRalign[i] <- alignPR(a[i,])
#}


#Add the alignments to the a data frame and save the data
dat <- cbind(dat, RTalign, PRalign)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)
write.table(dat, file = "../tmp/data_processed.txt", row.names = F, col.names = T)

#Convert that string of ATCG into a matrix of the first 300 amino acids
RTseqs <- matrix(data = NA, ncol = 900, nrow = nrow(dat))
for(i in 1:nrow(RTseqs)){
    if(i %% 1000 == 0){ print(i) }
    RTseqs[i,] <- (strsplit(dat$RTalign[i], split = "")[[1]])[1:900]
}
#Convert them to uppercase
RTseqcaps <- apply(RTseqs, 2, toupper)

#Num amino acids
RTaas <- 300
sites <- nrow(dat)
aaalign <- matrix(data = NA, nrow = sites, ncol = RTaas)

#Convert these into amino acids in a way that captures the diversity
for(i in 1:RTaas){
    if(i %% 10 == 0){ print(i) }
    aaalign[,i] <-(apply(RTseqcaps[,pos(i)], 1, translateAmbigNA))
}
#this is just slow

#Now, analagous processing for protease
#This handling is a little bit more complicated, because it needs to handle the case of missing protease sequence data

PRaas <- 297
PRseqs <- matrix(data = NA, ncol = PRaas, nrow = nrow(dat))
for(i in 1:nrow(PRseqs)){
    if(i %% 1000 == 0){ print(i) }
    PRseqs[i,] <- (strsplit(PRalign[i], split = "")[[1]])[1:297]
}
PRseqcaps <- apply(PRseqs, 2, toupper)

sitesPR <- 99
praaalign <- matrix(data = NA, ncol = sitesPR, nrow = nrow(PRseqcaps))
for(i in 1:sitesPR){
    if(i %% 10 == 0){ print(i) }
    praaalign[,i] <- apply(PRseqcaps[,pos(i)], 1, translateAmbigNA )    
}

#praaalign.new <- praaalign
#indval <- 30
#firstaa <- 2
#table(praaalign[which(dat$PRFirstAA == firstaa), indval] == praaalign.new[which(dat$PRFirstAA == firstaa), indval])



dim(PRseqcaps)
dim(praaalign)
dim(RTseqcaps)
dim(aaalign)


#Write to a data file, the PR and RT nucleotide and amino acid alignments
write.table(PRseqcaps, "../tmp/PR_nuc.txt", row.names= FALSE, col.names = FALSE, quote = FALSE)
write.table(praaalign, "../tmp/PR_aa.txt", row.names= FALSE, col.names = FALSE, quote = FALSE)
write.table(RTseqcaps, "../tmp/RT_nuc.txt", row.names= FALSE, col.names = FALSE, quote = FALSE)
write.table(aaalign, "../tmp/RT_aa.txt", row.names= FALSE, col.names = FALSE, quote = FALSE)


#####################################################################
# t1.3 Split the data into the clonal and non-clonal dataset
#
# This section splits the data into two datasets -
# The *nonclonal dataset* (called D-PCR) through the dataset features one seq/patient 
# The *clonal dataset* features multiple seqs/patient
#####################################################################


#Determine what is clonal and what is not

#Which Patient IDs appear more than once?
PtIDClones <- names(which((table(dat$PtID) > 1) == TRUE))
PtIDNonClones <- names(which((table(dat$PtID) == 1) == TRUE))

#Record the indices that have patientIDs that appear more than once
cloneIdent <- rep(FALSE, nrow(dat))
for(i in PtIDClones){
    cloneIdent[(which(dat$PtID == i))] = TRUE
}
#The rest of these are non-clonal
ncloneIdent <- (xor(rep(TRUE,nrow(dat)), cloneIdent))

#Shorten the final names that we use to differentiate clones and non-clones
clonal <- cloneIdent 
nclone <- ncloneIdent
nullnames <- which(dat$CloneName=="NULL")

clone <- setdiff(which(clonal == TRUE), nullnames)
#How many clonal sequences?
length(clone)

#Information about clonal sequences
summary(c(table(dat[clone,]$PtID)))
length(c(table(dat[clone,]$PtID)))


##################################################################
#t1.4 Call drug resistance mutations in each sequence
#
# This section calls DRMs and does some other housekeeping.
# For example, it also standardizes regimen names
##################################################################

#First, standardize all the regimen names so that they appear in alphabetical order, separated by a plus sign
tmp <- strsplit(as.vector(dat$RegimenName), split = '\\+')
newRegimen <- c()
numMeds <- c()
for(i in 1:length(tmp)){
    newRegimen[i] <- paste(sort(trim(tmp[[i]])), collapse = "+")
    numMeds[i] <- length(tmp[[i]])
}
head(newRegimen, n = 10)
is.factor(newRegimen)

dat$newRegimen <- newRegimen
dat$numMeds <- numMeds
#Add the new standardized regimen names, in addition to the number of meds that are present

#Handle the fixed DRM information
b <- read.table("../dat/DRM_file_WHO.txt", header = T, stringsAsFactors = F)

NNRTI <- "EFV|RPV|ETR|DLV|NVP|NNRTI"
NRTI <- "AZT|3TC|FTC|DDI|ABC|D4T|TDF|ZDV|DDC|^NRTI"
PI <- "TPV|IDV|SQV|LPV|FPV|DRV|ATV|NFV|APV|PI"

#Ok, we're now going to need to count DRMs
observedDRM <- list()

#How many nonclonal sequences do we have?
nclen <- length(which(nclone == TRUE))

#Set up lists in which we'll store the number of DRMs
hardDRMlist <- append(observedDRM, rep(0, nclen))
softDRMlist <- append(observedDRM, rep(0, nclen))

#We're only going to call DRMs in nonclonal sequences
prncaa <- (praaalign[nclone,]) # Protease
ncaa <- (aaalign[nclone,]) # RT

#For each sequence, we're going to try to decide what DRMs are present
for(i in 1:nclen){
#    i <- 4964
#    i <- which(dat[nclone,]$PtID == 118856)

    #Print occasionally 
    if(i%%500 == 0){ print(i)}

    #Lists to keep track of possible DRMs that we will ultimately add to the list
    toAddhard <- c()
    toAddsoft <- c()

    #Figure out what the treatment is for that patient
    treat <- dat[nclone,]$newRegimen[i]
#    treat <- dat[nclone,]$RegimenName[i]

    #Split apart the treatment into its component parts
    treatsplit <- strsplit(treat, split = "\\+")[[1]]
    
    #does it have PIs? NRTIs? NNRTIs?
    PI.all <- c()
    NRTI.all <- c()
    NNRTI.all <- c()
    for(k in 1:length(treatsplit)){
        PI.all[k] <- regexpr(PI, treatsplit[k])[1] > 0
        NRTI.all[k] <- regexpr(NRTI, treatsplit[k])[1] > 0
        NNRTI.all[k] <- regexpr(NNRTI, treatsplit[k])[1] > 0
    }
    PI.bool <- sum(PI.all) > 0
    NRTI.bool <- sum(NRTI.all) > 0
    NNRTI.bool <- sum(NNRTI.all) > 0

    #If it has NRTIs, test it for NRTI DRMs
    if(NRTI.bool){
        position <- b[b$Type == "NRTI", ]$Position
        postval <- b[b$Type == "NRTI", ]$Post

        for(j in 1:length(position)){
            observed.state <- ncaa[i,position[j]]
            post.list <- strsplit(postval[j], split = "")[[1]]

            #hard sweep
            hard.test <- grep(paste("^[(", postval[j], ")]$", sep = ""), observed.state, value = T)
            if(length(hard.test)>0){
                toAddhard <- c(toAddhard, paste("RT", position[j], hard.test, sep = ""))
            }

            #soft sweep
            soft.test <- grep(paste("^[", postval[j], "](\\/[",postval[j],"])+$", sep = ""), observed.state, value = T)
            if(length(soft.test)>0){
                toAddsoft <- c(toAddsoft, paste("RT", position[j],soft.test, sep = ""))
            }
        }
    }
    #If it has NNRTIs, test it for NNRTI DRMs
    if(NNRTI.bool){
       position <- b[b$Type == "NNRTI", ]$Position
       postval <- b[b$Type == "NNRTI", ]$Post

       for(j in 1:length(position)){

           observed.state <- ncaa[i,position[j]]
           post.list <- strsplit(postval[j], split = "")[[1]]

            #hard sweep
           hard.test <- grep(paste("^[(", postval[j], ")]$", sep = ""), observed.state, value = T)
           if(length(hard.test)>0){
               toAddhard <- c(toAddhard, paste("RT", position[j], hard.test, sep = ""))
           }
            
            #soft sweep
           soft.test <- grep(paste("^[", postval[j], "](\\/[",postval[j],"])+$", sep = ""), observed.state, value = T)
           if(length(soft.test)>0){
               toAddsoft <- c(toAddsoft, paste("RT", position[j],soft.test, sep = ""))
           }
        }
   }
   #If it has PIs, test it for PI DRMs
   if(PI.bool){

       position <- b[b$Type == "PI", ]$Position
       postval <- b[b$Type == "PI", ]$Post

       for(j in 1:length(position)){

#           dat[nclone,][which(prncaa[,73] == "T"),]$newRegimen
#           table(prncaa[,position[11]])
           
           observed.state <- prncaa[i,position[j]]
           post.list <- strsplit(postval[j], split = "")[[1]]

           #hard sweep
           hard.test <- grep(paste("^[(", postval[j], ")]$", sep = ""), observed.state, value = T)
           if(length(hard.test)>0){
               toAddhard <- c(toAddhard, paste("PI", position[j], hard.test, sep = ""))
           }
            
            #soft sweep
           soft.test <- grep(paste("^[", postval[j], "](\\/[",postval[j],"])+$", sep = ""), observed.state, value = T)
           if(length(soft.test)>0){
               toAddsoft <- c(toAddsoft, paste("PI", position[j],soft.test, sep = ""))
           }
       }

   }

   #Add hard and soft DRMs to the two lists 
   hardDRMlist[[i]] <- toAddhard
   softDRMlist[[i]] <- toAddsoft

    #If nothing is found, add 0
   if(length(toAddhard) == 0){hardDRMlist[[i]] <- 0}
   if(length(toAddsoft) == 0){softDRMlist[[i]] <- 0}

}


#this is a check. Let's look at an index where we think we've found a soft sweep
softlocs <- which(unlist(lapply(softDRMlist, function(x){x != 0})) == TRUE) # PSP these are the individuals (sequences) that have a soft sweep, right?
#PSP softlocs sounds like sites to me. 
ind <- softlocs[3]
softDRMlist[[ind]]
hardDRMlist[[ind]]

#Create a list that has both hard and soft sweeps
DRMlist <- list()
DRMlist <- append(DRMlist, rep(0, length(hardDRMlist)))

for(i in 1:length(hardDRMlist)){
    toAdd <-  c(hardDRMlist[[i]], softDRMlist[[i]])
    if(length(grep("^0", toAdd)) > 0){
        toAdd <- toAdd[-grep("^0", toAdd)]
    }
    DRMlist[[i]] <- toAdd
    if(length(toAdd) == 0){DRMlist[[i]] <- "0"}
}

#Numbers of hard and soft sweeps
hardDRMnum <- as.numeric(unlist(lapply(hardDRMlist, listlen)))
softDRMnum <- as.numeric(unlist(lapply(softDRMlist, listlen)))
DRMnum <- as.numeric(unlist(lapply(DRMlist, listlen)))


#how many RT ambiguous sites?
RTnumAmbig <- apply(ncaa, 1, function(x){ c(length(grep("\\/", x)))})
#How long is the RT sequence? (Not Ns or Xs)
RTlen <- apply(prncaa, 1, function(x){ c(length(grep("[XN]", x, invert = T)))})

#How many PR ambiguous sites?
PRnumAmbig <- apply(prncaa, 1, function(x){ c(length(grep("\\/", x)))})
#How many PR seqs are NA at the first position (these whole seqs will have len 0)
NAS <- which(is.na(prncaa[,1]))
#What is the length of the PR sequence?
PRlen <- apply(prncaa, 1, function(x){ c(length(grep("[XN]", x, invert = T)))})
#Store length 0 in the PR sequences with no reads
PRlen[NAS] = 0

#Revisit this
#What if remove all the drm sites?
drmsites.rt <- b[b$Location == "RT",]$Position
drmsites.pr <- b[b$Location == "PI",]$Position
RTnumAmbig <- apply(ncaa[, -drmsites.rt], 1, function(x){ c(length(grep("\\/", x)))})
RTlen <- apply(ncaa[,-drmsites.rt], 1, function(x){ c(length(grep("[XN]", x, invert = T)))})

#What if we remove all PR drm sites?
PRnumAmbig <- apply(prncaa[,-drmsites.pr ], 1, function(x){ c(length(grep("\\/", x)))})
nrpr.tmp <- prncaa
for(i in length(prncaa)){
    if(length(intersect(i, NAS)) < 0){
        nrpr.tmp[i,] <- prncaa[i, -drmsites.pr]
    }
}
PRlen <- apply(nrpr.tmp, 1, function(x){ c(length(grep("[XN]", x, invert = T)))})
PRlen[NAS] = 0


ambig <- (RTnumAmbig + PRnumAmbig)/(RTlen + PRlen)


#datold <- dat
dat <- (dat[nclone,])

#Add all of this extra information to the data list
dat$DRMnum <- DRMnum
dat$hardDRMnum <- hardDRMnum
dat$softDRMnum <- softDRMnum
dat$RTnumAmbig <- RTnumAmbig
dat$RTlen <- RTlen
dat$PRnumAmbig <- PRnumAmbig
dat$PRlen <- PRlen
dat$ambig <- as.numeric(ambig)


##################################################################
#t1.5 Final filters
#
# This section calls DRMs and does some other housekeeping.
# For example, it also standardizes regimen names
##################################################################



#Additional filters
#Here, we remove certain sequences for various different reasons

#This section is making a list of indices to be excluded based on being from non-ambiguous studies
#Studies that call no ambiguous reads in ANY of their sequences are excluded
#Studies that have any sequences with any reads that are ambiguous are included
percnoamb <- c()
studlen <- c()
#Create a list of unique studies
allstud <- unique(dat$MedlineID)

length(allstud)
#For each study, 
for( i in 1:length(allstud)){
    #look only at the data from a particular study
    subs <- dat[dat$MedlineID == allstud[i],]
    #Record the proportion of sequences that have any ambiguous reads
    percnoamb[i] <- sum(subs$RTnumAmbig ==  0)/nrow(subs)
    #And now many sequences there are from that study
    studlen[i] <- nrow(subs)
}

#noambigs are all studies with no sequences with any ambiguous reads
noambigs <- allstud[which(percnoamb == 1)]

toExclude <- grep(paste(noambigs, collapse = "|"), dat$MedlineID, invert = F)
#toExclude is all patients coming from studies that call no ambiguous reads


#We will also exclude patients receiving PIs with no sequenced protease
#If we include these patients, they will also have fewer DRMs called than are possibly actually there, since we won't be able to observe any DRMs in the unsequenced region
gotaPI <- which(regexpr(PI, dat$Regimen) > -1)
PRlen0 <- which(dat$PRlen == 0)
#intersection of those who got a PI and have PR length of 0
noPRexcludes <- intersect(gotaPI, PRlen0)

#Add these guys to the "toExclude" list
toExclude <- c(noPRexcludes , toExclude)

#Let's also exclude any patient who didn't get a well defined treatment

#What are all the treatments
alltreats <- unique(dat$Regimen)

#Filter out things we should definitely exclude:
#If any treatment is encoded just as "RTI"/"NRTI"/"NNRTI" or "PI" - we don't want to group things since we don't what they are
treats.filt <- alltreats[(regexpr('RTI', alltreats) < 0) & (regexpr('PI', alltreats) < 0)]
#Exclude single dose nevarapine - these patients are generally only treated once (not our study population)
#PSP I added if statements here to make sure that we don't remove everything when there is nothing to remove. 
if (length(which(treats.filt == "NVP"))>0)treats.filt <- treats.filt[-which(treats.filt == "NVP")]
#Exclude anything with the treatment "unknown"
if (length(which(regexpr("Unknown", treats.filt) >= 0))>0)treats.filt <- treats.filt[-which(regexpr("Unknown", treats.filt) >= 0)]
#Remove anything treated with aAPA (experimental treatment)
if (length(which(regexpr("aAPA", treats.filt) >= 0))>0)treats.filt <- treats.filt[-which(regexpr("aAPA", treats.filt) >= 0)]
#Many drugs
more.than.5.drugs <- which(unlist(lapply(strsplit(treats.filt, split = "\\+"), length) > 4))
if (length(more.than.5.drugs)>0) treats.filt <- treats.filt[-more.than.5.drugs]

#These are all the regimens we want to consider valid firstline therapies
treats.filt

#How much data do we have left, if we only include patients that were given one of these things?

#for each treatment on the list
allowableindices <- c()
for(i in 1:length(treats.filt)){
    #Record which indices map to this treatment
    allowableindices <- c(allowableindices, which(dat$Regimen == treats.filt[i]))
}
length(allowableindices)



#Ok, so allowable indices are our well-defined patients
allowableindices 

#toExclude are our patients coming from studies with no ambiguous reads called
toExclude

#remove toExclude from allowable
allow <- setdiff(allowableindices, toExclude)

#Ok, allow only these patients:
dat <- dat[allow,]

#We'll exclude these things from our hard/soft DRM lists too
hardDRMlist <- hardDRMlist[allow]
softDRMlist <- softDRMlist[allow]
DRMlist <- DRMlist[allow]

#And our alignments
ncaaalign <- ncaa[allow, ]
ncpraaalign <- prncaa[allow, ]
ncRTseqcaps <- (RTseqcaps[nclone,])[allow,]
ncPRseqcaps <- (PRseqcaps[nclone,])[allow,]

write.table(dat, "../tmp/nclonal.dat.txt", row.names= FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

save(hardDRMlist, file = "../tmp/hardDRMlist.txt") #PSP is it correct that these are things that are not human readable? 
save(softDRMlist, file = "../tmp/softDRMlist.txt")
save(DRMlist, file = "../tmp/DRMlist.txt")

write.table(ncaaalign, "../tmp/ncaaalign.txt", row.names= FALSE, col.names = TRUE, quote = FALSE)
write.table(ncpraaalign, "../tmp/ncpraaalign.txt", row.names= FALSE, col.names = TRUE, quote = FALSE)
write.table(ncRTseqcaps, "../tmp/ncRTseqcaps.txt", row.names= FALSE, col.names = TRUE, quote = FALSE)
write.table(ncPRseqcaps, "../tmp/ncPRseqcaps.txt", row.names= FALSE, col.names = TRUE, quote = FALSE)





