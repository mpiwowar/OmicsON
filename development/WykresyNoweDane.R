
# NOWE DANE
# Not Extended
# CCA
####################################################################################
####################################################################################
# TODO : Share knowladge with Sonia and Pani Monika
# Nowe dane - no extention
xNamesVector <- as.character(transcriptomicsInputData$HGNC)
yNamesVector <- as.character(lipidomicsInputData$Fatty_acids)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
# 0.3, 0.7

scalingFactor <- 1000
XDataFrame[seq(2,length(XDataFrame))] <- XDataFrame[seq(2,length(XDataFrame))] * scalingFactor
YDataFrame[seq(2,length(YDataFrame))] <- YDataFrame[seq(2,length(YDataFrame))] * scalingFactor
# 0.4, 0.7

noiceToData <- data.frame(matrix(0, nrow = dim(XDataFrame)[1], ncol = (dim(XDataFrame)[2]-1)))
for (i in seq(1, ncol(noiceToData))) {
    noiceToData[,i] <- round(rnorm(n = dim(XDataFrame)[1], sd = 3), digits = 0)
}
XDataFrame[seq(2,length(XDataFrame))] <- XDataFrame[seq(2,length(XDataFrame))] + noiceToData
noiceToData <- data.frame(matrix(0, nrow = dim(YDataFrame)[1], ncol = (dim(YDataFrame)[2]-1)))
for (i in seq(1, ncol(noiceToData))) {
    noiceToData[,i] <- round(rnorm(n = dim(YDataFrame)[1], sd = 3), digits = 0)
}
YDataFrame[seq(2,length(YDataFrame))] <- YDataFrame[seq(2,length(YDataFrame))] + noiceToData
# 0.3, 0.7

CcaResultsNoExtentionNewData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsNoExtentionNewData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "")

# NOTES : Time consumming!!! Our is selective.
# PCA
PlsResultsNoExtentionNewData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)

OmicsON::plotRegression(PLSResult = PlsResultsNoExtentionNewData)
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsNoExtentionNewData, nCols = 3, nRows = 1)


# Reactome Decoration
# Ensamble
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$genesSymbolsFromEnsemble)
yNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

CcaResultsReactomeEnsembleExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

CcaResultsReactomeEnsembleExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

# PCA
PlsResultsReactomeEnsembleExtentionNewData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)

OmicsON::plotRegression(PLSResult = PlsResultsReactomeEnsembleExtentionNewData)
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeEnsembleExtentionNewData, nCols = 3, nRows = 1)


# UniProt
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromUniProtFunctionalInteractions$genesSymbolsFromUniProt)
yNamesVector <- as.character(ontology2GenesSymboleFromUniProtFunctionalInteractions$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

CcaResultsReactomeUniProtExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeUniProtExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

CcaResultsReactomeUniProtExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeUniProtExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

# PCA
PlsResultsReactomeUniProtExtentionNewData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)
asd <- PlsResultsReactomeUniProtExtentionNewData$training$coefficients
pls::coefplot(PlsResultsReactomeUniProtExtentionNewData$training)
pls::compnames(PlsResultsReactomeUniProtExtentionNewData$training)
pls::respnames(PlsResultsReactomeUniProtExtentionNewData$training)
asd <- fitted(PlsResultsReactomeUniProtExtentionNewData$training)
plot(aaaaaaaa)
plot(pls::RMSEP(PlsResultsReactomeUniProtExtentionNewData$training["CHEBI:42504"]), legendpos = "topright")
OmicsON::plotRegression(PLSResult = PlsResultsReactomeUniProtExtentionNewData)
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeUniProtExtentionNewData, nCols = 3, nRows = 1)
aaaaa <- PlsResultsReactomeUniProtExtentionNewData$training
aaaaaaaa <- pls::RMSEP(PlsResultsReactomeUniProtExtentionNewData$training)
bbbbbb <- aaaaaaaa$val[]

aaaaaaaa[[1]]

aaaaaaaa$val <- aaaaaaaa$val[, 1:3, ]
plot(aaaaaaaa)
# STRING Decoration

# Narrow
# Ensamble
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsEnsemble$stringGenesSymbolsNarrow)
yNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsEnsemble$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

# PCA

# UniProt
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsUniProt$stringGenesSymbolsNarrow)
yNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsUniProt$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringUniProtExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringUniProtExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringUniProtExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")


# PCA


# Extended
# Ensamble
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsEnsemble$stringGenesSymbolsExpand)
yNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsEnsemble$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringEnsembleExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringEnsembleExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")


# PCA

# UniProt
# CCA
xNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsUniProt$stringGenesSymbolsExpand)
yNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsUniProt$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringUniProtExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.5, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringUniProtExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.3, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

# PCA



####################################################################################
####################################################################################
# Stare dane - no extention
xNamesVector <- as.character(transcriptomicsInputData$symbol)
yNamesVector <- as.character(lipidomicsInputData$ChEBI)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

DataUnderAnalysis <- list()
DataUnderAnalysis$OldData_CleanData$xData <- xNamesVector
DataUnderAnalysis$OldData_CleanData$yData <- yNamesVector

CcaResultsNoExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.8)
DataUnderAnalysis$OldData_Cca_CleanData_06_08$functionalGrouping <- CcaResultsNoExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_CleanData_06_08$correlationCutOff <- CcaResultsNoExtentionOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsNoExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "")

CcaResultsNoExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
DataUnderAnalysis$OldData_Cca_CleanData_06_07$functionalGrouping <- CcaResultsNoExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_CleanData_06_07$correlationCutOff <- CcaResultsNoExtentionOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsNoExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "")

PlsResultsNoExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0)
DataUnderAnalysis$OldData_Pls_CleanData_10_10$functionalGrouping <- PlsResultsNoExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_CleanData_10_10$correlationCutOff <- PlsResultsNoExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsNoExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsNoExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.8
)
DataUnderAnalysis$OldData_Pls_CleanData_06_08$functionalGrouping <- PlsResultsNoExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_CleanData_06_08$correlationCutOff <- PlsResultsNoExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsNoExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsNoExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7
)
DataUnderAnalysis$OldData_Pls_CleanData_06_07$functionalGrouping <- PlsResultsNoExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_CleanData_06_07$correlationCutOff <- PlsResultsNoExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsNoExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

# Stare dane - reactome extention - Ensemble
xNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$genesSymbolsFromEnsemble)
yNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

CcaResultsReactomeEnsembleExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.9, yCutoff = 0.5)
DataUnderAnalysis$OldData_Cca_Reactome_Ensemble_09_05$functionalGrouping <- CcaResultsReactomeEnsembleExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_Reactome_Ensemble_09_05$correlationCutOff <- CcaResultsReactomeEnsembleExtentionOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

CcaResultsReactomeEnsembleExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
DataUnderAnalysis$OldData_Cca_Reactome_Ensemble_06_07$functionalGrouping <- CcaResultsReactomeEnsembleExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_Reactome_Ensemble_06_07$correlationCutOff <- CcaResultsReactomeEnsembleExtentionOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")

PlsResultsReactomeEnsembleExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_10_10$functionalGrouping <- PlsResultsReactomeEnsembleExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_10_10$correlationCutOff <- PlsResultsReactomeEnsembleExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeEnsembleExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))


PlsResultsReactomeEnsembleExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.9, yCutoff = 0.5
)
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_09_05$functionalGrouping <- PlsResultsReactomeEnsembleExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_09_05$correlationCutOff <- PlsResultsReactomeEnsembleExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeEnsembleExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsReactomeEnsembleExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7,
    ncompValue = 7
)
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_06_07$functionalGrouping <- PlsResultsReactomeEnsembleExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_Ensemble_06_07$correlationCutOff <- PlsResultsReactomeEnsembleExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeEnsembleExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

# Stare dane - reactome extention - UniProt
xNamesVector <- as.character(ontology2GenesSymboleFromUniProtFunctionalInteractions$genesSymbolsFromUniProt)
yNamesVector <- as.character(ontology2GenesSymboleFromUniProtFunctionalInteractions$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

CcaResultsReactomeUniProtExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.8, yCutoff = 0.9
)
DataUnderAnalysis$OldData_Cca_Reactome_UniProt_08_09$functionalGrouping <- CcaResultsReactomeUniProtExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_Reactome_UniProt_08_09$correlationCutOff <- CcaResultsReactomeUniProtExtentionOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeUniProtExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")
CcaResultsReactomeUniProtExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7
)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeUniProtExtentionOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome")


PlsResultsReactomeUniProtExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_10_10$functionalGrouping <- PlsResultsReactomeUniProtExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_10_10$correlationCutOff <- PlsResultsReactomeUniProtExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeUniProtExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))


PlsResultsReactomeUniProtExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.8, yCutoff = 0.9
)
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_08_09$functionalGrouping <- PlsResultsReactomeUniProtExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_08_09$correlationCutOff <- PlsResultsReactomeUniProtExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeUniProtExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsReactomeUniProtExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7,
    ncompValue = 6
)
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_06_07$functionalGrouping <- PlsResultsReactomeUniProtExtentionOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_Reactome_UniProt_06_07$correlationCutOff <- PlsResultsReactomeUniProtExtentionOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsReactomeUniProtExtentionOldData,
                         nCols = 3, nRows = 2, lty = c(2))


# Stare dane - STRING extention - Ensemble - Expand
xNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsEnsemble$stringGenesSymbolsExpand)
yNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsEnsemble$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringEnsembleExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.7, yCutoff = 0.8)
DataUnderAnalysis$OldData_Cca_String_Expand_Ensemble_07_08$functionalGrouping <- CcaResultsStringEnsembleExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Expand_Ensemble_07_08$correlationCutOff <- CcaResultsStringEnsembleExtentionExpandOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

CcaResultsStringEnsembleExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
DataUnderAnalysis$OldData_Cca_String_Expand_Ensemble_06_07$functionalGrouping <- CcaResultsStringEnsembleExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Expand_Ensemble_06_07$correlationCutOff <- CcaResultsStringEnsembleExtentionExpandOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

PlsResultsStringEnsembleExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame
)
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_10_10$functionalGrouping <- PlsResultsStringEnsembleExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_10_10$correlationCutOff <- PlsResultsStringEnsembleExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringEnsembleExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsStringEnsembleExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.7, yCutoff = 0.8
)
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_07_08$functionalGrouping <- PlsResultsStringEnsembleExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_07_08$correlationCutOff <- PlsResultsStringEnsembleExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringEnsembleExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsStringEnsembleExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7
)
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_06_07$functionalGrouping <- PlsResultsStringEnsembleExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_Ensemble_06_07$correlationCutOff <- PlsResultsStringEnsembleExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringEnsembleExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))


# Stare dane - STRING extention - UniProt - Expand
xNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsUniProt$stringGenesSymbolsExpand)
yNamesVector <- as.character(ontology2GenesSymboleFromStringExpandFunctionalInteractionsUniProt$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringUniProtExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.7, yCutoff = 0.8)
DataUnderAnalysis$OldData_Cca_String_Expand_UniProt_07_08$functionalGrouping <- CcaResultsStringUniProtExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Expand_UniProt_07_08$correlationCutOff <- CcaResultsStringUniProtExtentionExpandOldData$correlationCutOff

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")
CcaResultsStringUniProtExtentionExpandOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
DataUnderAnalysis$OldData_Cca_String_Expand_UniProt_06_07$functionalGrouping <- CcaResultsStringUniProtExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Expand_UniProt_06_07$correlationCutOff <- CcaResultsStringUniProtExtentionExpandOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionExpandOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

PlsResultsStringUniProtExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0
)
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_10_10$functionalGrouping <- PlsResultsStringUniProtExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_10_10$correlationCutOff <- PlsResultsStringUniProtExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringUniProtExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsStringUniProtExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.7, yCutoff = 0.8
)
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_07_08$functionalGrouping <- PlsResultsStringUniProtExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_07_08$correlationCutOff <- PlsResultsStringUniProtExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringUniProtExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsStringUniProtExtentionExpandOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7
)
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_06_07$functionalGrouping <- PlsResultsStringUniProtExtentionExpandOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Expand_UniProt_06_07$correlationCutOff <- PlsResultsStringUniProtExtentionExpandOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringUniProtExtentionExpandOldData,
                         nCols = 3, nRows = 2, lty = c(2))


# Stare dane - STRING extention - Ensemble - Narrow
xNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsEnsemble$stringGenesSymbolsNarrow)
yNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsEnsemble$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0)
DataUnderAnalysis$OldData_Cca_String_Narrow_Ensemble_10_10$functionalGrouping <- CcaResultsStringEnsembleExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Narrow_Ensemble_10_10$correlationCutOff <- CcaResultsStringEnsembleExtentionNarrowOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

CcaResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")


PlsResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0
)
DataUnderAnalysis$OldData_Pls_String_Narrow_Ensemble_10_10$functionalGrouping <- PlsResultsStringEnsembleExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Narrow_Ensemble_10_10$correlationCutOff <- PlsResultsStringEnsembleExtentionNarrowOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringEnsembleExtentionNarrowOldData,
                         nCols = 3, nRows = 2, lty = c(2))

PlsResultsStringEnsembleExtentionNarrowOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7,
    ncompValue = 5
)
DataUnderAnalysis$OldData_Pls_String_Narrow_Ensemble_06_07$functionalGrouping <- PlsResultsStringEnsembleExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Narrow_Ensemble_06_07$correlationCutOff <- PlsResultsStringEnsembleExtentionNarrowOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringEnsembleExtentionNarrowOldData,
                         nCols = 3, nRows = 2, lty = c(2))


# Stare dane - STRING extention - UniProt - Narrow
xNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsUniProt$stringGenesSymbolsNarrow)
yNamesVector <- as.character(ontology2GenesSymboleFromStringNarrowFunctionalInteractionsUniProt$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
CcaResultsStringUniProtExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0)
DataUnderAnalysis$OldData_Cca_String_Narrow_UniProt_10_10$functionalGrouping <- CcaResultsStringUniProtExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Cca_String_Narrow_UniProt_10_10$correlationCutOff <- CcaResultsStringUniProtExtentionNarrowOldData$correlationCutOff
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

CcaResultsStringUniProtExtentionNarrowOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringUniProtExtentionNarrowOldData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "Reactome and STRING")

PlsResultsStringUniProtExtentionNarrowOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1.0, yCutoff = 1.0
)
DataUnderAnalysis$OldData_Pls_String_Narrow_UniProt_10_10$functionalGrouping <- PlsResultsStringUniProtExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Narrow_UniProt_10_10$correlationCutOff <- PlsResultsStringUniProtExtentionNarrowOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringUniProtExtentionNarrowOldData,
                         nCols = 3, nRows = 2, lty = c(2))


PlsResultsStringUniProtExtentionNarrowOldData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7,
    ncompValue = 6
)
DataUnderAnalysis$OldData_Pls_String_Narrow_UniProt_06_07$functionalGrouping <- PlsResultsStringUniProtExtentionNarrowOldData$functionalGrouping
DataUnderAnalysis$OldData_Pls_String_Narrow_UniProt_06_07$correlationCutOff <- PlsResultsStringUniProtExtentionNarrowOldData$correlationCutOff
OmicsON::plotRmsepForPLS(PLSResult = PlsResultsStringUniProtExtentionNarrowOldData,
                         nCols = 3, nRows = 2, lty = c(2))





write(DataUnderAnalysis$OldData_Cca_Reactome_Ensemble_06_07$functionalGrouping$xData, ncolumns = 1, file = "D:\\projects\\science\\OldData_Cca_Reactome_Ensemble_06_07functionalGroupingxData.txt")
save(DataUnderAnalysis, file = "D:\\projects\\science\\dataUnderAnalysis.data")
load(file = "D:\\projects\\science\\dataUnderAnalysis.data")

####################################################################################
####################################################################################
