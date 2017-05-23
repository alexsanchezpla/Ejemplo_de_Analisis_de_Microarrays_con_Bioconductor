## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(fig.path='images/grafic', tidy=FALSE, cache=FALSE)

## ----librerias, echo=FALSE, message=FALSE--------------------------------
installifnot <- function (pckgName){
 if (!(require(pckgName, character.only=TRUE))) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(eval(pckgName), suppressUpdates=TRUE)
  }else{
    print(paste("Package", pckgName, "already installed", sep=" "))
  } 
}
installifnot("Biobase")
installifnot("affy")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("limma")
installifnot("hgu133a.db")
installifnot("annotate")
installifnot("annaffy")
installifnot("hwriter")
installifnot("gplots")
installifnot("GOstats")
installifnot("GSA")

## ----preparaDirectorios, eval=TRUE---------------------------------------
workingDir <-getwd()
dataDir <-file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"results")
celfilesDir <- file.path(workingDir,"celfiles")
setwd(workingDir)

## ----tomaMuestras--------------------------------------------------------
# muestras <- read.csv2(file.path(dataDir, "Asignacion_de_muestras_a_grupos.csv"),  head=T)
# misMuestras <- as.character (muestras$Sample)
# paraAnalisis <- c(sample(misMuestras[1:6], 5), 
#                  sample(misMuestras[7:20], 5),
#                  sample(misMuestras[23:48], 5))
# alAnalisis <-muestras[muestras$Sample %in% paraAnalisis,]
# write.table(alAnalisis, file=file.path(dataDir, "targets.txt"), 
#            sep="\t", row.names=FALSE, quote=FALSE)

## ----phenoData1, echo=FALSE, results='asis', eval=TRUE-------------------
require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"), 
  header = TRUE, row.names = 1) 
stopifnot(require(xtable))
x.big<-xtable(pData(my.targets)[,1:4],
    caption='Archivo targets.txt con la asignación a cada muestra de su condición experimental')
print(x.big,tabular.environment='longtable',floating=FALSE)

## ----affybatch.create, eval=TRUE-----------------------------------------
require(affy)
sampleInfo <- read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"), 
    header = TRUE, row.names = 1, sep="\t")
fileNames <- rownames(pData(sampleInfo))
rawData <- read.affybatch(filenames=file.path(celfilesDir,fileNames),
                          phenoData=sampleInfo)
show(rawData)

## ----preajustes, echo=F, eval=TRUE---------------------------------------
colores <- c(rep("yellow", 5), rep("blue", 5), rep("red", 5))
grupos <- pData(rawData)$Group
numSamples <- nrow(pData(rawData))
sampleNames <-paste( pData(rawData)$SampleIDs, grupos, sep=".")
colnames(exprs(rawData))<-sampleNames

## ----plotHist,echo=F,eval=TRUE-------------------------------------------
hist(rawData, main="Signal distribution", col=colores, lty=1:numSamples)
legend (x="topright", legend=sampleNames , col=colores, lty=1:numSamples)

## ----boxplot,  eval=TRUE-------------------------------------------------
boxplot(rawData, cex.axis=0.6, col=colores, las=2, names=sampleNames, 
        main="Signal distribution for selected chips")

## ----plotPCA, eval=TRUE--------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=0.8)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

## ----plotPCA2D, eval=TRUE------------------------------------------------
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="selected samples")

## ----distAnalisis, eval=TRUE---------------------------------------------
  manDist <-  dist(t(exprs(rawData))) 
  heatmap (as.matrix(manDist),  col=heat.colors(16))  

## ----plotDendro, echo=F,eval=TRUE----------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)

## ----arrayQuality, cache=TRUE, warning=FALSE-----------------------------
stopifnot(require(arrayQualityMetrics))
arrayQualityMetrics(rawData, 
                    intgroup = "Group",
                    outdir = file.path(resultsDir, "arrayQuality"), 
                    force=TRUE)

## ----normalization.rma, eval=TRUE----------------------------------------
stopifnot(require(affy))
eset_rma <- rma(rawData)    

## ----normBoxPlot, eval=TRUE----------------------------------------------
boxplot(eset_rma,main="RMA", names=sampleNames, cex.axis=0.7, col=colores,las=2)

## ----variancePlot--------------------------------------------------------
sds <- apply (exprs(eset_rma), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Global variation over all genes")

## ----filtraje------------------------------------------------------------
require(genefilter)
filtered <- nsFilter(eset_rma, require.entrez=TRUE,
         remove.dupEntrez=TRUE, var.func=IQR,
         var.cutoff=0.5, var.filter=TRUE,
         filterByQuantile=TRUE, feature.exclude="^AFFX")

## ----filtrado------------------------------------------------------------
names(filtered)
class(filtered$eset)
print(filtered$filter.log)
eset_filtered <-filtered$eset

## ----writeNormalized-----------------------------------------------------
shortNames<- paste(pData(eset_rma)$Group, pData(eset_rma)$SampleIDs,sep=".")
sum(rownames(pData(eset_filtered))!=colnames(exprs(eset_filtered)))
colnames(exprs(eset_filtered))<- colnames(exprs(eset_rma))<- shortNames
write.csv2(exprs(eset_rma), file.path(resultsDir, "Datos.Normalizados.csv2"))
write.csv2(exprs(eset_filtered), file.path(resultsDir, "Datos.Normalizados.Filtrados.csv2"))
save(eset_rma, eset_filtered, file=file.path(resultsDir, "datos.normalizados.Rda"))

## ------------------------------------------------------------------------
eset2 <- eset_filtered[,1:10]
pData(eset2)
dim(exprs(eset2))

## ----teststat------------------------------------------------------------
stopifnot(require(genefilter))
teststat <-rowttests(eset2, "Group")
head(teststat)

## ----sortStatistics------------------------------------------------------
topDown<-order(teststat$p.value)
ranked<-teststat[topDown,]
print(top10<-ranked[1:10,])

## ----volcano1------------------------------------------------------------
x<-ranked$dm; y<--log(ranked$p.value)
plot(x, y, xlab="Fold Change", ylab ="-logPval", 
     main="Volcano plot\nA vs B")
abline(v=-1);abline(v=1);
text (x[1:10], y[1:10],rownames(ranked)[1:10], cex=0.7)

## ----selectNaif----------------------------------------------------------
selectedNaif <-ranked[ranked$p.value<0.01,]
nrow(selectedNaif)

## ----adjustPvalues-------------------------------------------------------
stopifnot(require(multtest))
#procs <- c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH", "BY")
procs <- c("Bonferroni","BH")
adjPvalues <- mt.rawp2adjp(ranked$p.value, procs)
names(adjPvalues)
ranked.adjusted<-cbind(ranked[,c(1,3)], adjPvalues$adjp[,-1])
ranked.adjusted[1:10,]

## ----selectedAdjusted----------------------------------------------------
selectedAdjusted<-ranked.adjusted[ranked.adjusted$BH<0.01,]
nrow(selectedAdjusted)
selectedAdjusted2<-ranked.adjusted[ranked.adjusted$BH<0.05,]
nrow(selectedAdjusted2)
selectedAdjusted2<-ranked.adjusted[ranked.adjusted$BH<0.25,]
nrow(selectedAdjusted2)

## ----matDesign, eval=TRUE------------------------------------------------
design<-matrix(
            c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,1,1,1,1,1),
            nrow=15,byrow=F)
colnames(design)<-c("A", "B", "L")
rownames(design) <-  sampleNames 
print(design)

## ----setContrasts, eval=TRUE---------------------------------------------
require(limma)
cont.matrix <- makeContrasts (
      AvsB = B-A,
      AvsL = L-A,
      BvsL = L-B,
      levels=design)
print(cont.matrix)

## ----linearmodelfit,echo=F, eval=TRUE------------------------------------
require(limma)
fit<-lmFit(eset_filtered, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)

## ----topTables, echo=TRUE, eval=TRUE-------------------------------------
topTab_AvsB <- topTable (fit.main, number=nrow(fit.main), coef="AvsB", adjust="fdr")
topTab_AvsL <- topTable (fit.main, number=nrow(fit.main), coef="AvsL", adjust="fdr")
topTab_BvsL  <- topTable (fit.main, number=nrow(fit.main) , coef="BvsL", adjust="fdr")

## ----volcanos, results='asis', echo=FALSE, eval=FALSE--------------------
## for(i in 1:ncol(cont.matrix)){
##   compName <-colnames(cont.matrix)[i]
##   file=paste("volcanoPlot", compName, ".pdf", sep="")
##   pdf(file=file.path(workingDir, "images", file), paper="special", width=6, height=6)
##   volcanoplot(fit.main, coef=i, highlight=10, names=rownames(fit.main),
##             main=paste("Differentially expressed genes",compName, sep="\n"))
##   abline(v=c(-1,1))
##   dev.off()
##   cat("\\includegraphics{", file, "}\n\n", sep="")
## }

## ----CuantosGenes, echo=F, eval=FALSE------------------------------------
## cat("Numero de genes con un p--valor inferior a 0.05 en cada comparacion:\n")
## cat ("En la comparacion 'A vs B': ", sum(topTab_AvsB$adj.P.Val<=0.05),"\n")
## cat ("En la comparacion 'A vs L': ", sum(topTab_AvsL$adj.P.Val<=0.05),"\n")
## cat ("En la comparacion 'B vs L': ", sum(topTab_BvsL$adj.P.Val<=0.05),"\n")
## 
## cat("\nNumero de genes con un p--valor inferior a 0.01 en cada comparacion:\n")
## cat ("En la comparacion 'A vs B': ", sum(topTab_AvsB$adj.P.Val<=0.01),"\n")
## cat ("En la comparacion 'A vs L': ", sum(topTab_AvsL$adj.P.Val<=0.01),"\n")
## cat ("En la comparacion 'B vs L': ", sum(topTab_BvsL$adj.P.Val<=0.01),"\n")

## ----topGenesAvsB, echo=FALSE, results='asis', eval=FALSE----------------
## require(Biobase)
## stopifnot(require(xtable))
## AvsB10<-xtable(topTab_AvsB[1:10,1:6],
##     label='topTab_AvsB',
##     caption='10 genes más expresados diferencialmente en la comparacion A vs B')
## print(AvsB10, tabular.environment='longtable',floating=FALSE)

## ----topGenesAvsL,echo=FALSE,results='asis', eval=FALSE------------------
## require(Biobase)
## stopifnot(require(xtable))
## AvsL10<-xtable(topTab_AvsL[1:10,1:6],
##     label='topTab_AvsL',
##     caption='10 genes más expresados diferencialmente en la comparacion A vs L')
## print(AvsL10, tabular.environment='longtable',floating=FALSE)

## ----topGenesBvsL,echo=FALSE,results='asis', eval=FALSE------------------
## BvsL10<-xtable(topTab_BvsL[1:10,1:6],
##     label='topTab_BvsL',
##     caption='10 genes más expresados diferencialmente en la comparacion B vs L')
## print(BvsL10, tabular.environment='longtable',floating=FALSE)

## ----decideTests.1, echo=F, eval=TRUE------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.01, lfc=1)

## ----resumeDecideTests, eval=TRUE----------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ----decideTests.2, echo=F, eval=TRUE------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.05, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ----venn1, eval=TRUE----------------------------------------------------
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)

## ----anota1, eval=TRUE---------------------------------------------------
require(hgu133a.db)
hgu133a()

## ----annaffy-------------------------------------------------------------
require(annaffy)
genesSelected <- rownames(res.selected)
at <- aafTableAnn(genesSelected, "hgu133a.db")
saveHTML (at, file.path(resultsDir, "anotations.html"), 
          "Annotations for selected genes")

## ----htmlPages-----------------------------------------------------------
require(annotate)
require(hgu133a.db)
listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
for (i in 1:length(listOfTables)){
  # Seleccionamos la "topTable"
  topTab <- listOfTables[[i]]
  # Escogemos los grupos de sondas a incluir en la tabla
  whichGenes<-topTab["P.Value"]<0.05
  selectedIDs <- rownames(topTab)[whichGenes]
  # Los convertimos a identificadores Entrez ("EG") y a Gene Symbols
  genes<- getEG(selectedIDs, "hgu133a.db")
  simbols <-getSYMBOL(selectedIDs, "hgu133a.db")
  # Haremos la columna de Entrez sea hiperenlazable
  paraEnlace <- list (misgenes=genes)
  # Preparamos el data.frame con el que se creará el archivo de resultados
  otherNames = data.frame(selectedIDs, simbols, topTab[whichGenes,-1])
  names(otherNames) = c("Affy ID", "Gene Symbol", colnames(topTab)[-1])
  # Invocamos la función "htmlpage"
  comparison <- names(listOfTables)[i]
  htmlpage(paraEnlace, 
           filename =file.path(resultsDir, 
           paste("Selected Genes in comparison ",comparison,".html", sep="")) , 
           title = paste("Diff. expressed genes in comparison ", comparison, sep=""), 
           othernames = otherNames, 
           table.head = c("Entrez IDs", names(otherNames)),
           table.center = TRUE, 
           repository=list("en"))
}

## ----expressAndTop-------------------------------------------------------
require(annotate)

EntrezsA <- getEG (rownames(topTab_AvsB), annotation(eset_rma))
SymbolsA <- getSYMBOL (rownames(topTab_AvsB), annotation(eset_rma))
ExpressAndTop_AvsB <- cbind(SymbolsA, EntrezsA, topTab_AvsB, exprs(eset_filtered)[rownames(topTab_AvsB), 1:10])
write.csv2(ExpressAndTop_AvsB, file.path(resultsDir, "ExpressAndTop_AvsB.csv2"))

EntrezsA <- getEG (rownames(topTab_AvsL), annotation(eset_rma))
SymbolsA <- getSYMBOL (rownames(topTab_AvsL), annotation(eset_rma))
ExpressAndTop_AvsL <- cbind(SymbolsA, EntrezsA, topTab_AvsL, exprs(eset_filtered)[rownames(topTab_AvsL), c(1:5,11:15)])
write.csv2(ExpressAndTop_AvsL, file.path(resultsDir, "ExpressAndTop_AvsL.csv2"))

EntrezsA <- getEG (rownames(topTab_BvsL), annotation(eset_rma))
SymbolsA <- getSYMBOL (rownames(topTab_BvsL), annotation(eset_rma))
ExpressAndTop_BvsL <- cbind(SymbolsA, EntrezsA, topTab_BvsL, exprs(eset_filtered)[rownames(topTab_BvsL), 6:15])
write.csv2(ExpressAndTop_BvsL, file.path(resultsDir, "ExpressAndTop_BvsL.csv2"))


## ----prepareData, eval=TRUE----------------------------------------------
probeNames<-rownames(res)
probeNames.selected<-probeNames[sum.res.rows!=0]
exprs2cluster <-exprs(eset_rma)[probeNames.selected,]
colnames(exprs2cluster)<-sampleNames
color.map <- function(grupo) { 
  if (grupo=="A"){
    c<- "yellow" 
  }else{ 
    if (grupo=="B"){
      c<- "red"
    }else{
      c<- "blue"
   }
  }
return(c)}

## ----plotHeatMap1, eval=TRUE---------------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
heatmap(exprs2cluster, col=rainbow(100), ColSideColors=grupColors, cexCol=0.9)

## ----plotHeatMap2,  eval=TRUE--------------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
require("gplots")
heatmap.2(exprs2cluster, 
          col=bluered(75), scale="row",
          ColSideColors=grupColors, key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexCol=1)

## ----GOAnalysis----------------------------------------------------------
require(GOstats)
listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
for (i in 1:length(listOfTables)){
  # Seleccionamos la "topTable"
  topTab <- listOfTables[[i]]
  # Definimos el universo de genes: todos los que se han incluido en el análisis
  # EL programa trabaja con identificadores "entrez" y no admite duplicados
  
  entrezUniverse = unique(getEG(rownames(topTab), "hgu133a.db"))
  
  # Escogemos los grupos de sondas a incluir en el análisis
  # Este análisis trabaja bien con varios centenares de genes 
  # por lo que es habitual basarse en p-valores sin ajustar para incluirlos
  
  whichGenes<-topTab["adj.P.Val"]<0.05
  geneIds <-   unique(getEG(rownames(topTab)[whichGenes],"hgu133a.db"))
  
  # Creamos los "hiperparámetros" en que se basa el análisis
  GOparams = new("GOHyperGParams",
    geneIds=geneIds, universeGeneIds=entrezUniverse,
    annotation="org.Hs.eg.db", ontology="BP",
    pvalueCutoff=0.001, conditional=FALSE,
    testDirection="over")
  KEGGparams = new("KEGGHyperGParams",
    geneIds=geneIds, universeGeneIds=entrezUniverse,
    annotation="org.Hs.eg.db",  
    pvalueCutoff=0.01, testDirection="over")

  # Ejecutamos los análisis

  GOhyper = hyperGTest(GOparams)
  KEGGhyper = hyperGTest(KEGGparams)
  
# Creamos un informe html con los resultados
   comparison = names(listOfTables)[i]
   GOfilename =file.path(resultsDir, 
     paste("GOResults.",comparison,".html", sep=""))
   KEGGfilename =file.path(resultsDir, 
     paste("KEGGResults.",comparison,".html", sep=""))
  htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
  htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))
}

## ----preparaGenestes-----------------------------------------------------
require(annotate)
require(hgu133plus2.db)
geneSets <- read.csv(file.path(dataDir, "PBTs_all_affy.csv"))

## gene sets a comparar amb cadascuna de les toptable
gsEntrez <- list()
for (i in 1:dim(geneSets)[2]) {
  gs <- as.character(geneSets[,i])
  gs <- gs[gs != ""]
  gsE <- getEG(gs, data = "hgu133plus2.db")
  gsEntrez[[colnames(geneSets[i])]] <- unique(gsE[!is.na(gsE)])
}

gsSymbol <- list()
for (i in 1:dim(geneSets)[2]) {
  gs <- as.character(geneSets[,i])
  gs <- gs[gs != ""]
  gsSymb <- getSYMBOL(gs,data = "hgu133plus2.db")
  gsSymbol[[colnames(geneSets[i])]] <- unique(gsSymb[!is.na(gsSymb)])
}
genesetsNames <- names(gsEntrez)

## ----prepareGSEA---------------------------------------------------------
require(genefilter)
esetUnique <- featureFilter(eset_rma, require.entrez=TRUE, remove.dupEntrez=FALSE)
dim(exprs(esetUnique))
esetUnique <- featureFilter(eset_rma, require.entrez=TRUE, remove.dupEntrez=TRUE)
dim(exprs(esetUnique))
esetAB<- esetUnique[,pData(esetUnique)$Group %in% c("A","B")]
pData(esetAB)
xAB <- exprs(esetAB)
yAB <- ifelse(pData(esetAB)$Group=="A", 1, 2)
entrezs<-getEG(rownames(xAB), "hgu133a")
simbols <-getSYMBOL(rownames(xAB), "hgu133a")

## ----execGSA-------------------------------------------------------------
require(GSA)
GSA.obj<-GSA(xAB, yAB, genenames=entrezs, genesets=gsEntrez,  resp.type="Two class unpaired", nperms=100)
GSA.listsets(GSA.obj, geneset.names=genesetsNames, FDRcut=.5)

## ----GSAbyKEGG, eco=FALSE, message=FALSE---------------------------------
geneset.obj <- GSA.read.gmt(file.path(dataDir,"c2.cp.kegg.v5.1.entrez.gmt"))
GSA.obj2 <- GSA(xAB, yAB, genenames=entrezs, genesets=geneset.obj$genesets,  resp.type="Two class unpaired", nperms=100)
GSA.listsets (GSA.obj2, geneset.names=geneset.obj$geneset.names, FDRcut=.5)

## ----listaArchivos, echo=FALSE, results='asis', eval=TRUE----------------
require(gdata)
listaArchivos <- dir(resultsDir) 
data.frame(listaArchivos)

## ----listaArchivos2html,echo=FALSE, eval=TRUE----------------------------
require(hwriter)
hwrite(listaArchivos,file.path(resultsDir, "listaArchivos.html"))

