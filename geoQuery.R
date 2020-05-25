
inVHIR <- TRUE
if (inVHIR) {
  if (!(require(httr))) install.packages("httr")
  set_config(
    use_proxy(url = "conf_www.ir.vhebron.net", port=8081)
  )
}
if (!require(GEOquery))
  { source("http://bioconductor.org/biocLite.R")
    biocLite("GEOquery")
  }
require(GEOquery)
gse <- getGEO("GSE16254")
names(gse)
a <- getGEOSuppFiles('GSE16254')