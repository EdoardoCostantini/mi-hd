
# Install BiocMangaer pacakge
# source: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
# source: https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
# if(!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

rm(list = ls())

library(Biobase)
library(GEOquery)

# open an existing GDS file (even if its compressed):
gds3289 <- getGEO(filename = "../data/GDS3289_full.soft.gz")
gds3289 <- getGEO(filename = "../data/GDS3289.soft.gz") # from disk
gds3289 <- getGEO(GEO = "GDS3289", AnnotGPL = TRUE) # from web

# source: https://www.ncbi.nlm.nih.gov/sites/GDSbrowser [GDS3289 data link]
#         https://pubmed.ncbi.nlm.nih.gov/17173048/     [paper]

# Get to know data
gds3289_dt <- Table(gds3289)
class(gds3289_dt)
dim(gds3289_dt)
names(gds3289_dt)
gds3289_dt$IDENTIFIER
sum(grepl("GSM", colnames(gds3289_dt)))
sum(grepl("GSM", gds3289_dt$IDENTIFIER))

biom_predictors <- c(
  which(grepl("FAM178A", gds3289_dt$IDENTIFIER )),
  which(grepl("IMAGE:813259", gds3289_dt$IDENTIFIER )),
  which(grepl("UGP2", gds3289_dt$IDENTIFIER ))
)

# For missing study
GEO_dt <- as.data.frame(t(gds3289_dt[, -c(1:2)]))

GEO_dt[1:5, 1:5]
colnames(GEO_dt) <- gds3289_dt$IDENTIFIER
GEO_dt[, 1:5]

colMeans(is.na(GEO_dt[, c("FAM178A", "IMAGE:813259", "UGP2")]))

colMeans(is.na(GEO_dt[, biom_predictors]))