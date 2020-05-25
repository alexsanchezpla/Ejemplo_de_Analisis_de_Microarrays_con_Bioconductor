
if (!require(GEOquery))
  { source("http://bioconductor.org/biocLite.R")
    biocLite("GEOquery")
  }
require(GEOquery)
listOfStudies <- c("45547", "62564", "3446", "73517", "16476")
for (studyID in listOfStudies){
  gse <- getGEO(paste0("GSE",studyID), GSEMatrix =TRUE, AnnotGPL=TRUE)
  # class(gse)
  # class(gse[[1]])
  eset <- gse[[1]]
  phenoDat <- pData(eset)
  dim(phenoDat)
  # colnames(phenoDat)
  require(xlsx)
  write.xlsx(phenoDat, file="phenoData.xlsx", sheetName=paste0("GSE",studyID), append=TRUE)
}

