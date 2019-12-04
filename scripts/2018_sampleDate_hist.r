
args <- commandArgs(trailingOnly = T)

inputFile = args[1]
outputFile = args[2]

#inputFile <- "C:/Users/Emma/bash/github/enterovirus_d68/data/20190611_Karolinska-region.csv"
#outputFile <- "C:/Users/Emma/bash/github/2018_evd68_paneurope_analysis/figures/supp-sampleHist.pdf"

swed_meta <- read.csv(inputFile, sep=",", as.is=T)

date_yr <- swed_meta$date
date_yr <- gsub("-XX", "", date_yr)
date_yr <- gsub("-[0-9]+", "", date_yr)

met2018 <- swed_meta[grep("2018", swed_meta$date), ]

# exclude ESP_022 as it had poor coverage in 1 fragment
met2018 <- met2018[grep("ESP_022", met2018$strain, invert=T), ]

#Now have 53, same as final count in paper.

#get hist for all to get breaks etc
bh <- hist(as.Date(met2018$date), breaks="weeks")

#x axis labels are (convert from days to dates):
brksAsDates <- as.Date(bh$breaks, origin="1970-01-01")

#get for each country
met2018$country <- as.factor(met2018$country)
sampCnts <- matrix(nrow=length(levels(met2018$country)), ncol=length(bh$counts), 0)
rownames(sampCnts) <- c(levels(met2018$country))
colnames(sampCnts) <- format(brksAsDates, format="%d %b")[1:length(brksAsDates)-1] #breaks are always 1 longer than counts!

for (cn in levels(met2018$country)) {
    cndat <- met2018[which(met2018$country==cn),"date"]
    bt <- hist(as.Date(cndat), breaks=as.Date(bh$breaks, origin="1970-01-01"))
    sampCnts[cn,] <- bt$counts
}

#library("viridis")
#cols <- plasma(5) #also viridis magma plasma inferno https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/

#these colours come from R library 'viridis' the 'plasma' one - but installing this is nuts.
#So I'm just stealing the Hex colors from using it on my own computer:
cols <-  c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF")

#barplot(sampCnts, beside=T, xlab="Week Of", ylab="# of Samples",
#    col=cols, legend=rownames(sampCnts),
#    main="Sample Date of 53 Included 2018 Samples")

pdf(outputFile, width=9, height=4.5)
x <- barplot(sampCnts, beside=T, xlab="Week Beginning", ylab="# of Samples",
    col=cols, legend=rownames(sampCnts),
    main="Sample Date of 53 Included 2018 Samples", xaxt="n")
labs <- colnames(sampCnts)
text(cex=1, x=colMeans(x)+1.25, y=-0.55, labs, xpd=T, srt=-45)
dev.off()
