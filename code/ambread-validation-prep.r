#Prep for the validation

#We need to re-read in the whole dataset, because we use nonclonal reads
all <- read.table("~/Dropbox/HIV_softsweeps/dataset.01.24.txt", skip = 13, sep= "\t", header = T)

#mults is all ids that are present more than once in the db
ids <- table(all$PtID)
mults <- names(ids[which(ids > 1)])

#At the nucleotide level

### Within patient (nonclonal) - propambi
propambig.nu <- c()

#For every position,
dim(ncRTseqcaps)
for(i in 1:ncol(ncRTseqcaps)){

    #Make a table of all present nucleotides
    tmp <- table(subset(ncRTseqcaps[,i], ncRTseqcaps[,i] != "." & ncRTseqcaps[,i] != "N" ))
    
    #Calculate the proportion of all positions that are ambiguous
    propambig.nu[i] <- sum(tmp[grep('[ATCG]',names(tmp), invert = T)])/sum(tmp)
}

#Nuc level 
heterozyg.nuc <- function(x){
    bar <- x[x != "." & x != "N" ]
    #If not enough reads, return NA for this patient
    if(length(bar) <= 1){ return(NA) }
    #Otherwise, return the entropy
    tmp <- table(bar)
    #pi
    return(1-sum((tmp/sum(tmp))^2))
}

clonalDiv.nu <- matrix(data = NA, nrow = length(mults), ncol = ncol(ncRTseqcaps))
print(length(mults))
for(j in 1:length(mults)){
    #Record the clonal diversity
    if(j %% 10 == 0){ print(j) }
    clonalDiv.nu[j, ] <- apply(RTseqcaps[all$PtID == mults[j], ], 2, heterozyg.nuc)
}


clonalMeans.nu <- apply(clonalDiv.nu, 2, function(x){mean(x, na.rm = T)} ) 

#plot(propambig.nu, clonalMeans.nu, log = "xy")
cor.test(propambig.nu, clonalMeans.nu)

#AA level

propambig.aa <- c()
ambfreeentrop.aa <-c()

for(i in 1:ncol(ncaaalign)){
    tmp <- table(subset(ncaaalign[,i], ncaaalign[,i] != "." & ncaaalign[,i] != "X" ))
    propambig.aa[i] <- sum(tmp[grep('\\/',names(tmp))])/sum(tmp)
}

clonalDiv.aa <- matrix(data = NA, nrow = length(mults), ncol = ncol(ncaaalign))
#Here, we need a function that properly handles X amino acids
heterozyg.aa <- function(x){
    bar <- x[x != "X"]
    if(length(bar) <= 1){ return(NA) }
    tmp <- table(bar)
    return(1-sum((tmp/sum(tmp))^2))
}

for(j in 1:length(mults)){
    if(j %% 10 == 0){ print(j) }
    clonalDiv.aa[j, ] <- apply(aaalign[all$PtID == mults[j], ], 2, heterozyg.aa)    
}
clonalMeans.aa <- apply(clonalDiv.aa, 2, function(x){mean(x, na.rm = T)} )

#plot(propambig.aa, clonalMeans.aa, log = "xy")
cor.test(propambig.aa, clonalMeans.aa)

