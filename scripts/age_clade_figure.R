# See the original 'age_summary_regression.r' script for all exploratory commands!

args <- commandArgs(trailingOnly = T)

inFile = args[1]
outFile = args[2]
suppFigBoot = args[3]
suppFigClade = args[4]

#setwd("C:/Users/Emma/bash/github/2018_evd68_paneurope_analysis")
#inFile = "../enterovirus_d68/vp1/results/metadata_subgenotype_2018y.tsv"
#outFile = "figures/age_clade_plot-fig3.pdf"
#suppFigBoot = "figures/supp-BootstrapAgeDist.pdf"
#suppFigClade = "figures/supp-age_by_clade.pdf"

agedat <- read.csv(inFile, sep="\t", as.is=T)


someAge <- length(which(agedat$age_range3 != ""))
# have some kind of age data

ageRange1 <- length(which(agedat$age_range1 != ""))
#age range 1 data

ageRange2 <- length(which(agedat$age_range2 != ""))
#age range 2 data

exactAge <- length(which(!is.na(agedat$age)))
# have exact age data

someAgeSomeClade <- length(which(agedat$subgenogroup != "" & agedat$age_range3 != ""))
# have some kind of age data and some kind of clade

exactAgeSomeClade <- length(which(agedat$subgenogroup != "" & !is.na(agedat$age)))
# exact age data and some kind of clade

cladeA <- c("A", "A2", "A1")
cladeC <- c("C")
cladeB <- c("B", "B1", "B2", "B3")

cA_age <- agedat[which(agedat$subgenogroup %in% cladeA),]
cC_age <- agedat[which(agedat$subgenogroup %in% cladeC),]
cB_age <- agedat[which(agedat$subgenogroup %in% cladeB),]

numA <- nrow(cA_age)
numAAge <- length(which(cA_age$age_range3 != ""))
#  clade A have age data
numC <- nrow(cC_age)
numCAge <- length(which(cC_age$age_range3 != ""))
#  clade C
numB <- nrow(cB_age)
numBAge <- length(which(cB_age$age_range3 != ""))
#  clade B


age18plus <- length(which(agedat$age_range3 == ">=18y"))
# 18 or over
under18 <- length(which(agedat$age_range3 == "<18y"))
# under 18


cat(c (
paste("There are",nrow(agedat),"samples total"),
paste("  ",someAge,"have some kind of age data (age range 3)"),
paste("  ",ageRange2,"have age range 2 data"),
paste("  ",ageRange1,"have age range 1 data"),
paste("  ",exactAge,"have some exact age data"),
paste("  ",someAgeSomeClade,"have some kind of age data and have a clade assignment"),
paste("  ",exactAgeSomeClade,"have exact age data and have a clade assignment"),
paste("  ",numAAge,"/",numA,"in clade A have some kind of age data"),
paste("  ",numBAge,"/",numB,"in clade B have some kind of age data"),
paste("  ",numCAge,"/",numC,"in clade C have some kind of age data"),
paste("  ",age18plus,"(",round(age18plus/someAge*100,2),"%) of those with age data are 18 or older"),
paste("  ",under18,"(",round(under18/someAge*100,2),"%) of those with age data are under 18")
), sep="\n")


numB <- length(which(agedat$subgenogroup == "B"))
numBWithAge <- length(which(agedat$subgenogroup == "B" & agedat$age_range3 != ""))
numBover18 <- length(which(agedat$subgenogroup == "B" & agedat$age_range3 == ">=18y"))
numBunder18 <- length(which(agedat$subgenogroup == "B" & agedat$age_range3 == "<18y"))

#take out A1, A2, B1, B2, B3, C
wantClade <- c("A1", "A2", "B1", "B2", "B3", "C")
cMeta <- agedat[which(agedat$subgenogroup %in% wantClade),]
#only with ages
cMeta <- cMeta[which(!is.na(cMeta$age)),]

cat("\n\n\n")
cat(c(
paste("Subgenogroup 'B' (those not classed as B1, B2, or B3) contains only",numB,"sequences."),
paste(numBover18,"are 18 or over, and",numBunder18,"are under 18."),
paste("These",numB,"'B' sequences will be excluded from the regression analysis"),
paste("A sequences (not in A1 or A2) are also excluded"),
paste(""),
paste(nrow(cMeta),"(exact age, and clade (excluding 'B')) will be used for regression analysis\n")
), sep="\n")


#################################################
####   ACTUAL REGRESSION ANALYSIS AND AGE STATS
#################################################

nowAge <- cMeta$age
nowClade <- cMeta$subgenogroup

nowClade <- factor(nowClade, levels = c("B3", "B2", "B1", "A2", "A1", "C"))

numB3 <- length(which(cMeta$subgenogroup=="B3"))
numB2 <- length(which(cMeta$subgenogroup=="B2"))
numB1 <- length(which(cMeta$subgenogroup=="B1"))
numA2 <- length(which(cMeta$subgenogroup=="A2"))
numA1 <- length(which(cMeta$subgenogroup=="A1"))
numC <- length(which(cMeta$subgenogroup=="C"))

meanB3 <- mean(cMeta$age[which(cMeta$subgenogroup=="B3")])
meanB2 <- mean(cMeta$age[which(cMeta$subgenogroup=="B2")])
meanB1 <- mean(cMeta$age[which(cMeta$subgenogroup=="B1")])
meanA2 <- mean(cMeta$age[which(cMeta$subgenogroup=="A2")])
meanA1 <- mean(cMeta$age[which(cMeta$subgenogroup=="A1")])
meanC <- mean(cMeta$age[which(cMeta$subgenogroup=="C")])

CIB3 <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="B3")])/sqrt(length(which(cMeta$subgenogroup=="B3")))
CIB2 <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="B2")])/sqrt(length(which(cMeta$subgenogroup=="B2")))
CIB1 <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="B1")])/sqrt(length(which(cMeta$subgenogroup=="B1")))
CIA2 <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="A2")])/sqrt(length(which(cMeta$subgenogroup=="A2")))
CIA1 <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="A1")])/sqrt(length(which(cMeta$subgenogroup=="A1")))
CIC <- qnorm(0.975)*sd(cMeta$age[which(cMeta$subgenogroup=="C")])/sqrt(length(which(cMeta$subgenogroup=="C")))

iqrB3 <- IQR(cMeta$age[which(cMeta$subgenogroup=="B3")])
iqrB2 <- IQR(cMeta$age[which(cMeta$subgenogroup=="B2")])
iqrB1 <- IQR(cMeta$age[which(cMeta$subgenogroup=="B1")])
iqrA2 <- IQR(cMeta$age[which(cMeta$subgenogroup=="A2")])
iqrA1 <- IQR(cMeta$age[which(cMeta$subgenogroup=="A1")])
iqrC <- IQR(cMeta$age[which(cMeta$subgenogroup=="C")])

cat(c( "\n",
paste("B3:  n=",numB3,",\tmean age of",round(meanB3,2),"+-",round(CIB3,2),",\tIQR of",iqrB3,"            "),
paste("B2:  n=",numB2,",\tmean age of",round(meanB2,2),"+-",round(CIB2,2),",\tIQR of",iqrB2,"            "),
paste("B1:  n=",numB1,",\tmean age of",round(meanB1,2),"+-",round(CIB1,2),",\tIQR of",iqrB1,"            "),
paste("A2:  n=",numA2,",\tmean age of",round(meanA2,2),"+-",round(CIA2,2),",\tIQR of",iqrA2,"            "),
paste("A1:  n=",numA1,",\tmean age of",round(meanA1,2),"+-",round(CIA1,2),",\tIQR of",iqrA1,"            "),
paste("C:   n=",numC,", \tmean age of",round(meanC,2),"+-",round(CIC,2),",\tIQR of",iqrC,"            ")
),sep="\n")

cat("\n\n###########################\n")
cat("#### LINEAR REGRESSION ####\n")
cat("###########################\n")
summary(lm(nowAge~nowClade))
cat("###########################\n")

cols <- c("#67a0cb","#3E58CF","#DCAB3C","#E67932","#DC2F24","#009900")
wantyears <- c(2014, 2016, 2018)

year <- gsub("-XX", "", agedat$date)
year <- as.numeric(gsub("-[0-9]+", "", year))
agedat <- cbind(agedat, year)

########################################################################
########################################################################
#           PLOT AGE AND CLADE TOGETHER IN BAR PLOT
########################################################################
########################################################################

#interesting breakdown of % in age category by clade
#prop.table(table(cMeta$age_range2, cMeta$subgenogroup), margin=2)*100

#####
# get ages with unknown

evenYear <- agedat[which(agedat$year %in% wantyears),]
evenYear$age_range1 <- as.factor(evenYear$age_range1)
evenYear$age_range1 <- factor(evenYear$age_range1, levels(evenYear$age_range1)[c(5,2,4,3,6,1)])
counts <- table(evenYear$age_range1, evenYear$year)
counts <- prop.table(counts, margin=2)*100
ageCol <- c("#d7191c", "#fdae61", "#ffd11a", "#abd9e9", "#2c7bb6", "white")

#####
# make table of both clade and age
cladeC <- table(evenYear$subgenogroup, evenYear$year)
cladeC <- prop.table(cladeC, margin=2)*100
colClade <- c("#DC2F24", "#E67932", "#DCAB3C", "#3E58CF", "#67a0cb")
#barplot(cladeC, xlab="Year", ylab="Percent of samples", legend=T, col=colClade)

cladeC2 <- rbind(cladeC, matrix(0, nrow=6, ncol=3))
rownames(cladeC2) <- c(rownames(cladeC), rownames(counts))
colnames(cladeC2) <- paste(colnames(cladeC2), "Clade")

counts2 <- rbind(matrix(0, nrow=5, ncol=3), counts)
rownames(counts2) <- c(rownames(cladeC), rownames(counts))
colnames(counts2) <- paste(colnames(counts2), "Age")

cladeCounts <- cbind(counts2, cladeC2)
cladeCounts <- cladeCounts[,order(colnames(cladeCounts))]

############
# ACTUAL PLOT

#new age colours to try...
newAgeCol <- c("#7b3294", "#c2a5cf", "#b4c0b7", "#a6dba0", "#008837", "white")

#barplot(cladeCounts, xlab="Year", ylab="Percent of samples", legend=T,
#    col=c(ageCol,colClade))

pdf(outFile, 
    width=12, height=5)
layout(matrix(1:2,nrow=1),widths=c(3,4))
par(xpd=T)

denT = 4
plot(density(cMeta$age[which(cMeta$subgenogroup=="C")],    from=0, bw=denT), type="l", col=cols[6], lty=1, lwd=2, 
    zero.line=F, ylab="Density", xlab="Age", main="")
points(density(cMeta$age[which(cMeta$subgenogroup=="B3")], from=0, bw=denT), type="l", col=cols[5], lty=1, lwd=2)
points(density(cMeta$age[which(cMeta$subgenogroup=="B2")], from=0, bw=denT), type="l", col=cols[4], lty=2, lwd=2)
points(density(cMeta$age[which(cMeta$subgenogroup=="B1")], from=0, bw=denT), type="l", col=cols[3], lty=1, lwd=2)
points(density(cMeta$age[which(cMeta$subgenogroup=="A2")], from=0, bw=denT), type="l", col=cols[2], lty=1, lwd=2)
points(density(cMeta$age[which(cMeta$subgenogroup=="A1")], from=0, bw=denT), type="l", col=cols[1], lty=1, lwd=2)
legend("topright", levels(as.factor(cMeta$subgenogroup)), col=cols,
    lty=c(1,1,1,2,1,1), lwd=2)
text(-11,0.097, "A", cex=3)

bb <- barplot(cladeCounts, xlab="Year", ylab="Percent of samples", #legend=T,
    col=c(cols[1:5],newAgeCol), xlim=c(0,9),
    space=c(0.1, rep(c(0.2,0.4),2), 0.2), names=rep("",6))
legend(7.8, 110, c(rownames(counts)[1:5],"Unknown"), fill=newAgeCol, title="Age")
legend(7.8, 50, rownames(cladeC), fill=cols, title="Clade")
text(bb, -5, rep(c("Age", "Clade"),3))
text(c(mean(bb[1:2]), mean(bb[3:4]), mean(bb[5:6])), -13, c("2014", "2016", "2018"))
text(-.7,117, "B", cex=3)

par(xpd=F)
dev.off()

cat ("",
paste("Figure main age-clade figures has been written to",outFile)
, sep="\n")


###############################
#### BOOTSTRAPPING
###############################

cat("\nBootstrap running.... this may take a few moments.\n")

cntryChoices <- unique(cMeta$country)

pA2 <- c() #store p values from linear regression for A2 clade
pA1 <- c() # store p values from linear regression for A1 clade
mm <- c() #holds the bootstrapped samples

for (i in 1:100) {

    btrsp <- cMeta[which(cMeta$country=="not"),] #make an empty table
    #pick the countries we are going to include for this bootstrap, with replacement
    smp <- sample(cntryChoices, length(cntryChoices), replace=T)

    #add these countries' samples to our dataset
    for(cy in smp) {
        btrsp <- rbind(btrsp, cMeta[which(cMeta$country==cy),])
    }
    #reorder factors so that the linear regression table always in the right order!
    btrsp$subgenogroup <- factor(as.factor(btrsp$subgenogroup), levels=c("B3", "B2", "B1", "A2", "A1","C"))
    sy <- summary(lm(btrsp$age~btrsp$subgenogroup))
    #store the p values
    pA2 <- c(pA2, sy[[4]]["btrsp$subgenogroupA2",4]) #[4,4]) #A2
    pA1 <- c(pA1, sy[[4]]["btrsp$subgenogroupA1",4]) #sy[[4]][5,4]) #A1

    mm[[i]] <- btrsp[which(btrsp$subgenogroup %in% c("A2","B3","A1")),]
}

A2_sig <- length(which(pA2<0.05))
A1_sig <- length(which(pA1<0.05))

A2_bonf_sig <- length(which(pA2<(0.05/length(pA2))))
A1_bonf_sig <- length(which(pA1<(0.05/length(pA1))))

cat(c( "",
"100 Replicates performed",
"A2:",
paste("\tNumber < 0.05:",A2_sig),
paste("\tNumber < 0.05 after Bonf. correct:",A2_bonf_sig),
"A1:",
paste("\tNumber < 0.05:",A1_sig),
paste("\tNumber < 0.05 after Bonf. correct:",A1_bonf_sig)
), sep="\n")

pdf(suppFigBoot, 
    width=5, height=5)
denT = 4
plot(density(mm[[1]][which(mm[[1]]$subgenogroup=="B3"),"age"], from=0, bw=denT), type="l", col=cols[5],
    lwd=1, xlab="Age", ylab="Density", ylim=c(0,0.08), zero.line=F, 
    main="Age Distribution by Clade of 100 Bootstraps")
for (i in 1:100){
    points(density(mm[[i]][which(mm[[i]]$subgenogroup=="B3"),"age"], from=0, bw=denT), 
        col=cols[5], type="l", lwd=1)
    points(density(mm[[i]][which(mm[[i]]$subgenogroup=="A2"),"age"], from=0, bw=denT), 
        col=cols[2], type="l", lty=2)
    points(density(mm[[i]][which(mm[[i]]$subgenogroup=="A1"),"age"], from=0, bw=denT), 
        col=cols[1], type="l", lty=1)
}
legend("topright", c("B3","A2","A1"), col=c(cols[5], cols[2:1]),
    lty=c(1,2,1), lwd=2)
dev.off()

cat ("",
paste("Figure of the replicate density distribution has been written out to",suppFigBoot)
, sep="\n")

#############################################
# Compare if age is diff in different years
#############################################

year <- gsub("-XX", "", agedat$date)
year <- as.numeric(gsub("-[0-9]+", "", year))
agedat <- cbind(agedat, year)

evenYear <- agedat[which(agedat$year %in% wantyears),]

cat("\n\n###########################\n")
cat("#### AGE YEAR DIFFERENCE ####\n")
cat("###########################\n")
summary(lm((evenYear$age~as.factor(evenYear$year))))
cat("###########################\n")

###########################################
#### Supp graph - plot age_range1 % per clade
###########################################

cMeta$age_range1 <- as.factor(cMeta$age_range1)
cMeta$age_range1 <- factor(cMeta$age_range1, levels(cMeta$age_range1)[c(4,1,3,2,5)])

ageClade <- prop.table(table(cMeta$age_range1, cMeta$subgenogroup), margin=2)*100
ageCol2 <- c("#d7191c", "#fdae61", "#ffd11a", "#abd9e9", "#2c7bb6")

cladeNums <- c(
    length(which(cMeta$subgenogroup=="A1")),
    length(which(cMeta$subgenogroup=="A2")),
    length(which(cMeta$subgenogroup=="B1")),
    length(which(cMeta$subgenogroup=="B2")),
    length(which(cMeta$subgenogroup=="B3")),
    length(which(cMeta$subgenogroup=="C"))  )

pdf(suppFigClade, 
    width=7, height=5)
par(xpd=T)
bb2 <- barplot(ageClade, xlab="Clade", ylab="Percent of samples", legend=T,
    col=ageCol2, xlim=c(0,8.7), space=c(0.1, rep(0.2,5)), names=rep("",6) )

text(bb2, -5, c("A1", "A2", "B1", "B2", "B3", "C"))
text(bb2, -12, paste("N=",cladeNums, sep=""))
par(xpd=F)
dev.off()

cat ("",
paste("Age_range1 of each clade has been written out to",suppFigClade)
, sep="\n")

