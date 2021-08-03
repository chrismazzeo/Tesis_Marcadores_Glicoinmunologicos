#http://bioconductor.org/packages/release/bioc/vignettes/DaMiRseq/inst/doc/DaMiRseq.pdf
#feature selection
#imput raw counts + info extra por ejemplo datos cliniclos o inmuno infiltracion
# resultado: un set genes(los mas informativos) que van a servir como entrdad para modelos de prediccion


library(DaMiRseq)
library(readr)
source("./src/shared/functions/HELPER_Function.r")
source("./src/shared/functions/DEA_Function.r")

project.resourcePath = "./resources/mouse/LGALS1_KO-6KO/"  
samplesListPath = "./resources/mouse/LGALS1_KO-6KO/LGALS1_KO-6KO_sampleList.csv"
samplesList = HELPER_LoadCSV(samplesListPath)
htseq.rawCountMatrix = DESeq2_MergeHTSeqExpressionFiles(paste0(project.resourcePath,"htseq"),samplesList$HTSEQ_FILE)

samplestList.DaMir = samplesList
rownames(samplestList.DaMir) = samplesList$SAMPLE
samplestList.DaMir = samplestList.DaMir[,c(4,6)]

colnames(samplestList.DaMir)[1] = "class"
colnames(htseq.rawCountMatrix) = samplesList$SAMPLE
se = DaMiR.makeSE(x = htseq.rawCountMatrix, y = samplestList.DaMir)
data_norm <- DaMiR.normalization(se , minCounts=10, fSample=0.8,
                                 hyper = "yes", th.cv=3)


data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.9)
sv <- DaMiR.SV(data_filt, method = "be")
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=dim(sv)[2])
#assay(data_adjust[c(1:5), c(1:5, 21:25)])
DaMiR.Allplot(data_filt, colData(data_filt))
DaMiR.Allplot(data_adjust, colData(data_adjust))

set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.4)
data_reduced <- DaMiR.FReduct(data_reduced$data)
DaMiR.MDSplot(data_reduced, df)
df.importance <- DaMiR.FSort(data_reduced, df)
head(df.importance)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,
                                 n.pred = 10)
DaMiR.Clustplot(selected_features$data, df)
Classification_res <- DaMiR.EnsembleLearning(selected_features$data,
                                             classes=df$class, fSample.tr = 0.5,
                                             fSample.tr.w = 0.5, iter = 30)
