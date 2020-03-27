library(metagenomeSeq)

#reading a biom file
library(biomformat)
biom_file <-system.file("extdata", "min_sparse_otu_table.biom",package = "biomformat")
b <-read_biom(biom_file)
biom2MRexperiment(b)

#loading count data
dataDirectory <-system.file("extdata", package = "metagenomeSeq")
lung =loadMeta(file.path(dataDirectory, "CHK_NAME.otus.count.csv"))
dim(lung$counts)

#loading taxanomy
taxa =read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"),stringsAsFactors = FALSE)

#loading metadata
clin =loadPhenoData(file.path(dataDirectory, "CHK_clinical.csv"),tran = TRUE)
ord =match(colnames(lung$counts),rownames(clin))
clin = clin[ord, ]
head(clin[1:2, ])


#creating MRexperiment object
phenotypeData =AnnotatedDataFrame(clin)
phenotypeData

OTUdata =AnnotatedDataFrame(taxa)
OTUdata

obj =newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
obj

#lung data
data(lungData)
lungData

phenoData(obj)
head(pData(obj), 3)

#to acess feature information
featureData(obj)
head(fData(obj)[, -c(2, 10)], 3)

#to acess normalized counts matrix with MR counts function
head(MRcounts(obj[, 1:2]))

#subsetting MRexperiment class object
featuresToKeep =which(rowSums(obj) >= 100)
samplesToKeep =which(pData(obj)$SmokingStatus == "Smoker")
obj_smokers = obj[featuresToKeep, samplesToKeep]
obj_smokers

head(pData(obj_smokers), 3)

#Alternative method using normFactors
head(normFactors(obj))
normFactors(obj) <-rnorm(ncol(obj))
head(normFactors(obj))

#libsizemethods
head(libSize(obj))
libSize(obj) <-rnorm(ncol(obj))
head(libSize(obj))

#calculate normalization factors
data(lungData)
p =cumNormStatFast(lungData)
lungData =cumNorm(lungData, p = p)

#exportingdata
mat =MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, file =file.path(dataDirectory, "tmp.tsv"))

exportStats(lungData[, 1:5], file =file.path(dataDirectory,"tmp.tsv"))
head(read.csv(file =file.path(dataDirectory, "tmp.tsv"), sep = "\t"))

#using featurefit model for differential abundance testing
data(lungData)
lungData = lungData[, -which(is.na(pData(lungData)$SmokingStatus))]
lungData =filterData(lungData, present = 30, depth = 1)
lungData <-cumNorm(lungData, p = 0.5)
pd <-pData(lungData)
mod <-model.matrix( ~1 + SmokingStatus, data = pd)
lungres1 =fitFeatureModel(lungData, mod)
head(MRcoefs(lungres1))

##using fitZig

data(lungData)
controls =grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[, -controls]
rareFeatures =which(rowSums(MRcounts(lungTrim) > 0) < 10)
lungTrim = lungTrim[-rareFeatures, ]
lungp =cumNormStat(lungTrim, pFlag = TRUE, main = "Trimmed lung data")

lungTrim =cumNorm(lungTrim, p = lungp)

smokingStatus =pData(lungTrim)$SmokingStatus
bodySite =pData(lungTrim)$SampleType
normFactor =normFactors(lungTrim)
normFactor =log2(normFactor/median(normFactor) + 1)
`mod =model.matrix( ~smokingStatus + bodySite + normFactor)
settings =zigControl(maxit = 10, verbose = TRUE)
fit =fitZig(obj = lungTrim, mod = mod, useCSSoffset = FALSE,control = settings)

# maxit=1 is for demonstration purposes

settings =zigControl(maxit = 1, verbose = FALSE)
mod =model.matrix(~bodySite)
colnames(mod) =levels(bodySite)
# fitting the ZIG model
res =fitZig(obj = lungTrim, mod = mod, control = settings)

zigFit =slot(res, "fit")
finalMod =slot(res, "fit")$design
contrast.matrix =makeContrasts(BAL.A - BAL.B, OW - PSB, levels = finalMod)
fit2 =contrasts.fit(zigFit, contrast.matrix)
fit2 =eBayes(fit2)
topTable(fit2)
##exporting fits
taxa =sapply(strsplit(as.character(fData(lungTrim)$taxa), split = ";"),function(i){i[length(i)]})
head(MRcoefs(fit, taxa = taxa, coef = 2))

#log normal permutation test
coeffOfInterest = 2
res =fitLogNormal(obj = lungTrim, mod = mod, useCSSoffset = FALSE,B = 10, coef = coeffOfInterest)
adjustedPvalues =p.adjust(res$p, method = "fdr")
# extract the absolute 
foldChange =abs(res$fit$coef[, coeffOfInterest])
sigList =which(adjustedPvalues <= 0.05)
sigList = sigList[order(foldChange[sigList])]

head(taxa[sigList])

##feature specific

head(MRtable(fit, coef = 2, taxa = 1:length(fData(lungTrim)$taxa)))
patients =sapply(strsplit(rownames(pData(lungTrim)), split = "_"),function(i){i[3]})
pData(lungTrim)$patients = patientsclassIndex =list(smoker =which(pData(lungTrim)$SmokingStatus =="Smoker"))
classIndex =list(smoker =which(pData(lungTrim)$SmokingStatus =="Smoker"))
classIndex$nonsmoker =which(pData(lungTrim)$SmokingStatus =="NonSmoker")
otu = 779
plotOTU(lungTrim, otu = otu, classIndex, main = "Neisseria meningitidis")
x =fData(lungTrim)$taxa[otu]
otulist =grep(x,fData(lungTrim)$taxa)

plotGenus(lungTrim, otulist, classIndex, labs = FALSE, main = "Neisseria meningitidis")
lablist <-c("S", "NS")
axis(1, at =seq(1, 6, by = 1), labels =rep(lablist, times = 3))




