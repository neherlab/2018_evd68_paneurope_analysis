
args <- commandArgs(trailingOnly = T)

inFile = args[1]
outFile = args[2]

#inFile = "C:/Users/Emma/bash/github/enterovirus_d68/vp1/results/metadata-ages.tsv"
#outFile = "C:/Users/Emma/bash/github/2018_evd68_paneurope_analysis/results/reduced_country_meta.tsv"

meta <- read.csv(inFile, as.is=T, sep="\t")

countries <- c("France", "Spain", "Germany", "Italy", "Sweden",
               "france", "spain", "germany", "italy", "sweden")
regions <-c("china", "europe", "north_america")

meta$country <- tolower(meta$country)
meta$region <- tolower(meta$region)


#country
new_country <- rep("rest_world", nrow(meta))

new_country[which(meta$region=="europe")] <- "rest_europe"

new_country[which(meta$country %in% countries)] <- meta$country[which(meta$country %in% countries)]

#region
new_region <- rep("rest_world", nrow(meta))

new_region[which(meta$region %in% regions)] <- meta$region[which(meta$region %in% regions)]

#add them
meta$country <- new_country
meta$region <- new_region

write.table(meta, outFile, sep="\t", quote=F, row.names=F)


