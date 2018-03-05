source("http://bioconductor.org/biocLite.R")
biocLite("affy")

library(affy)
affy.data = ReadAffy()
eset.mas5 = mas5(affy.data)

exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)
# Rename the column names if we want
colnames(exprSet.nologs) = c("brain.1", "brain.2", 
                             "fetal.brain.1", "fetal.brain.2",
                             "fetal.liver.1", "fetal.liver.2", 
                             "liver.1", "liver.2")
exprSet = log(exprSet.nologs, 2)
write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")

# Run the Affy A/P call algorithm on the CEL files we processed above
data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
data.mas5calls.calls = exprs(data.mas5calls)
# Print the calls as a matrix
write.table(data.mas5calls.calls, file="Su_mas5calls.txt", quote=F, sep="\t")

targetsFile <- "estrogen/estrogen.txt"
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
pData(pd)

raw <-ReadAffy(celfile.path = "estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw
boxplot(raw,col="red",las=2)
par(mfrow=c(2,1))
hist(log2(pm(raw[,1])),breaks=100,col="steelblue",main="PM",xlim=c(4,14))
hist(log2(mm(raw[,1])),breaks=100,col="steelblue",main="MM",xlim=c(4,14))
mva.pairs(pm(raw)[,1:4],plot.method="smoothScatter")
mva.pairs(pm(raw)[,5:8],plot.method="smoothScatter")
biocLite("affyPLM")
library(affyPLM)
plmset <- fitPLM(raw)
NUSE(plmset,las=2)
RLE(plmset,las=2)

bad <- ReadAffy(celfile.path = "estrogen/",filenames="bad.cel")
image(bad)

par(mfrow=c(2,4))
image(raw[,1])
image(raw[,2])
image(raw[,3])
image(raw[,4])
image(raw[,5])
image(raw[,6])
image(raw[,7])
image(raw[,8])

eset <- rma(raw)
eset
head(exprs(eset))
