##############################
##############################
## ABRF 2017 Satellite workshop - Hands-on 4 : Statistical methods for high-throughput biology
## Date: March 25, 2017
## Created by Meena Choi and Ting Huang
##############################
##############################

load(file='Section4.RData')

##############################
##############################
## 1. MS Proteomics, iPRG intensity data
##############################
##############################

##############################
## Load MSstats
##############################

library(MSstats)
?MSstats

##############################
## Read data
##############################

# read skyline output
raw <- read.csv(file="iPRG_10ppm_2rt_15cut_nosingle.csv")

# read annotation table
annot <- read.csv("iPRG_skyline_annotation.csv", header=TRUE)
annot

# reformating and pre-processing for Skyline output.
quant <- SkylinetoMSstatsFormat(raw, annotation=annot)
head(quant)

##############################
## Processing data
##     Transformation = log2
##     Normalization = equalize median
##     Model-based run-level summarization (TMP) after imputation
##############################

quant.processed <- dataProcess(raw = quant, 
                               logTrans=2, 
                               normalization = 'equalizeMedians',
                               summaryMethod = 'TMP', 
                               MBimpute=TRUE,
                               censoredInt='0',
                               cutoffCensored='minFeature',
                               maxQuantileforCensored = 0.999)

# show the name of outputs
names(quant.processed)

# show reformated and normalized data.
# 'ABUNDANCE' column has normalized log2 transformed intensities.
head(quant.processed$ProcessedData)

# This table includes run-level summarized log2 intensities. (column : LogIntensities)
# Now one summarized log2 intensities per Protein and Run.
# NumMeasuredFeature : show how many features are used for run-level summarization.
#         If there is no missing value, it should be the number of features in certain protein.
# MissingPercentage : the number of missing features / the number of features in certain protein.
head(quant.processed$RunlevelData)

# show which summarization method is used.
head(quant.processed$SummaryMethod)

##############################
## Data visualization
##############################

dataProcessPlots(data = quant.processed, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 1,
                 address='iPRG_skyline_equalizeNorm_')

# It will generate profile plots per protein. It will take a while
# Please run at home. It takes a while.
dataProcessPlots(data = quant.processed, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 summaryPlot = TRUE,
                 originalPlot = FALSE,
                 address="iPRG_skyline_equalizeNorm_")

# Instead, make profile plot for only some.
dataProcessPlots(data = quant.processed, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P44015|VAC2_YEAST',
                 address="iPRG_skyline_equalizeNorm_P44015")


# sp|P55752|ISCB_YEAST
# sp|P44374|SFG2_YEAST
# sp|P44983|UTR6_YEAST
# sp|P44683|PGA4_YEAST
# sp|P55249|ZRT4_YEAST

# and first few proteins


# Please run at home. It takes a while.
dataProcessPlots(data = quant.processed, 
                 type="conditionplot", 
                 width=7, height=7,
                 address="iPRG_skyline_equalizeNorm_")

# Instead, make profile plot for only some.

# sp|P44015|VAC2_YEAST
# sp|P55752|ISCB_YEAST
# sp|P44374|SFG2_YEAST
# sp|P44983|UTR6_YEAST
# sp|P44683|PGA4_YEAST
# sp|P55249|ZRT4_YEAST

# and first few proteins


##############################
## Model-based comparison and adjustment for multiple testing
##############################

unique(quant.processed$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

test <- groupComparison(contrast.matrix=comparison, data=quant.processed)

names(test)

# Show test result
# Label : which comparison is reported.
# log2FC : estimated log2 fold change between conditions.
# adj.pvalue : adjusted p value by BH
# issue : detect whether this protein has any issue for comparison
#    such as, there is measurement in certain group, or no measurement at all.
# MissingPercentage : the number of missing intensities/total number of intensities 
#     in conditions your are interested in for comparison
# ImputationPercentage : the number of imputed intensities/total number of intensities 
#     in conditions your are interested in for comparison
head(test$ComparisonResult)

# After fitting linear model, residuals and fitted values can be shown.
head(test$ModelQC)

# Fitted model per protein
head(test$fittedmodel)

# save testing result as .csv file
Skyline.intensity.comparison.result <- test$ComparisonResult

write.csv(Skyline.intensity.comparison.result, 
          file='testResult_iprg_skyline.csv')

head(Skyline.intensity.comparison.result)
SignificantProteins <- Skyline.intensity.comparison.result[Skyline.intensity.comparison.result$adj.pvalue < 0.05 ,]
nrow(SignificantProteins)


##############################
## Visualization for testing result
##############################

groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="VolcanoPlot", 
                     address="testResult_iprg_skyline_")

groupComparisonPlots(data = Skyline.intensity.comparison.result, 
                     type = 'VolcanoPlot',
                     sig = 0.05, 
                     FCcutoff = 2^2, 
                     address = 'testResult_iprg_skyline_FCcutoff4_')

groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="Heatmap", 
                     address="testResult_iprg_skyline_")

# Please run at home. It takes a while.
groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="ComparisonPlot", 
                     address="testResult_iprg_skyline_")


##############################
## Verify the model assumption
##############################

# normal quantile-quantile plots
modelBasedQCPlots(data=test, type="QQPlots", 
                  width=5, height=5, 
                  address="iprg_skyline_")

# residual plots
modelBasedQCPlots(data=test, type="ResidualPlots", 
                  width=5, height=5, 
                  address="iprg_skyline_")



##############################
## Power calculation
##############################
test.power <- designSampleSize(data = test$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = TRUE,
                               numSample = 3)
test.power

designSampleSizePlots(data = test.power)



##############################
## Sample size calculation
##############################
samplesize <- designSampleSize(data = test$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = 0.9,
                               numSample = TRUE)
samplesize

designSampleSizePlots(data = samplesize)



##############################
## sample quantification
##############################
sampleQuant <- quantification(quant.processed)
head(sampleQuant)




##############################
## extra : quantile normalization
##############################
quant.processed.quantile <- dataProcess(raw = quant, 
                               logTrans=2, 
                               normalization = 'quantile',
                               summaryMethod = 'TMP', 
                               MBimpute=TRUE,
                               censoredInt='0',
                               cutoffCensored='minFeature',
                               maxQuantileforCensored = 0.999)

dataProcessPlots(data = quant.processed.quantile, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 1,
                 address='iPRG_skyline_quantile_')

dataProcessPlots(data = quant.processed.quantile, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 1,
                 address="iPRG_skyline_quantile_1_")

##############################
## extra : Limma
##############################

library(limma)

## use run-level summarized value from MSstats
input <- quant.processed$RunlevelData

## reformat
input2 <- dcast(Protein ~ originalRUN, data=input, value.var = 'LogIntensities')
head(input2)

## annotate protein id in rowname
rownames(input2) <- input2$Protein
input2 <- input2[, -1]

design <- model.matrix(~0+factor(c(1,1,1,2,2,2,3,3,3,4,4,4)))
colnames(design) <- c('Condition1', 'Condition2', 'Condition3', 'Condition4')
contrast.matrix <- makeContrasts(Condition2-Condition1, levels=design)

fit <- lmFit(input2, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

test.limma <- data.frame(Label='C2-C1',
                         log2FC=fit2$coefficients,
                         pvalue=fit2$p.value)
head(test.limma)
colnames(test.limma)[2] <- 'log2FC'
colnames(test.limma)[3] <- 'pvalue'
test.limma$adj.pvalue <- p.adjust(test.limma$pvalue, method="BH")
test.limma$Protein <- rownames(test.limma)

save(test.limma, file='test.limma.RData')


##############################
##############################
## 2. MS Proteomics, iPRG spectral count data
##############################
##############################


## Reformat
library(reshape2)

Y <- dcast(Protein ~ Run, data=iprg.count)
head(Y)

countData <- as.matrix(Y[,-c(1)])
protName <- as.character(Y[,1])      

##############################
## DESeq2: NB GLM + Wald test
##############################

library(DESeq2)

## Format the data for DESeq
## use M1 as reference in the design matrix
colData <- data.frame(M1=factor(c(rep(1, 3),rep(0, 9))), 
                      M2=factor(c(rep(0, 3), rep(1,3), rep(0, 6))),
                      M3=factor(c(rep(0, 6), rep(1,3),rep(0, 3))), 
                      M4=factor(c(rep(0, 9), rep(1,3))))

## Comparisons with 1st one  

## Design model ~ M2+M3+M4
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ M2+M3+M4)

## Test for differential expression with NB distribution and exact test
dds <- DESeq(dds)
dds # object of DESeqDataSet class 
resultsNames(dds)

## block automatic independent filtering after testing
c1vs2 <- results(dds, contrast=c('M2', '0', '1'), independentFiltering = FALSE)
c1vs3 <- results(dds, contrast=c('M3', '0', '1'), independentFiltering = FALSE)
c1vs4 <- results(dds, contrast=c('M4', '0', '1'), independentFiltering = FALSE)

rm(dds)


## Combine everything
results_count <- data.frame(rbind(c1vs2, c1vs3, c1vs4))
results_count$Label <- factor(c(rep(c("C2-C1"),dim(Y)[1]),rep(c("C3-C1"),
                                                              dim(Y)[1]),rep(c("C4-C1"),dim(Y)[1])))
results_count$Protein <- Y$Protein

results_count <- results_count[, c(8,7,1,2,3,4,5,6)]
head(results_count)
colnames(results_count)[4] <- 'log2FC'
colnames(results_count)[8] <- 'adj.pvalue'

## Save the result 
test.count.deseq2 <- results_count
save(test.count.deseq2, file='test.count.deseq2.RData')


##############################
## Venn Diagram
##############################

pdf('VennDiagram_intensity_iprg_msstat_p_or_adjp.pdf')
venn(list(adjusted.pvalue = test.msstats[test.msstats$adj.pvalue < 0.05, 'Protein'],
          pvalue = test.msstats[test.msstats$pvalue < 0.05, 'Protein']))
dev.off()

pdf('VennDiagram_count_iprg_adj.pvalue.pdf')
venn(list(ttest = test.count.ttest[test.count.ttest$adj.pvalue < 0.05, 'Protein'],
          Deseq2 = test.count.deseq2[test.count.deseq2$adj.pvalue < 0.05, 'Protein']))
dev.off()
