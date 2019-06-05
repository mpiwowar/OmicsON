# library(CCA)
# data(nutrimouse)
# X=as.matrix(nutrimouse$gene[,1])
# Y=as.matrix(nutrimouse$lipid)
# dim(X)
# dim(Y)
# res.cc=cc(X,Y)
#
# plot(res.cc$cor,type="b")
# plt.cc(res.cc)


################## TEST NOWYCH DANYCH ##################
# PARAGRAPH: Set up mapping files
OmicsON::setUpReactomeMapping(
    ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt",
    Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt",
    UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")

# PARAGRAPH: Data input
pathToFileWithLipidomicsData <- paste("C:/Users/Cezary/Downloads", "fatty_acid_KDOdoC.txt",
                                      sep = "/")
# UWAGI do vignetty: No spaces in colnames and rownames.
# UWAGI do vignetty: Use . as separator in decimal, not come ,.
lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)

pathToFileWithTranscriptomicsData <- paste("C:/Users/Cezary/Downloads", "HGNC_KDOdoC.txt",
                                           sep = "/")
transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)

# PARAGRAPH: Decorate data by Reactome data
decoratedByReactome <- OmicsON::decorateByReactomeData(chebiMoleculesDf = lipidomicsInputData,
                                                       chebiIdsColumnName = "Fatty_acids",
                                                       organismTaxonomyId = '9606')


# PARAGRAPH:
# UWAGI do vignetty: Some accepted errors which are not stopping processing; Error in as.igraph.vs(graph, v) : Invalid vertex names
# UWAGI do vignetty: It takes long time! For 30 rows it took 30 min.
decoratedByStringBaseOnEnsembleIds <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactome,
    listOfEnsembleIdColumnName = 'ensembleIds')
# PROCESSING: It takes also long time.
decoratedByStringBaseOnUniProtIds <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactome,
    listOfEnsembleIdColumnName = 'uniProtIds')

# PARAGRAPH: Functional Ineractions DF
ontology2EnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactome,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'ensembleIds')
ontology2UniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    decoratedByReactome,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'uniProtIds')
ontology2GenesSymboleFromEnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    decoratedByReactome,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'genesSymbolsFromEnsemble')
ontology2GenesSymboleFromUniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    decoratedByReactome,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'genesSymbolsFromUniProt')
# PROCESSING: Very big set, 280 000 rows!!!
ontology2GenesSymboleFromStringExpandFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    decoratedByStringBaseOnEnsembleIds,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'stringGenesSymbolsExpand')
ontology2GenesSymboleFromStringNarrowFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    decoratedByStringBaseOnEnsembleIds,
    singleIdColumnName = 'ontologyId',
    idsListColumnName = 'stringGenesSymbolsNarrow')


# PARAGRAPH: CCA - Canonical Correlation Analysis
# UWAGI do kodu: Add this to CCA's function code and to PLS's function code.
commonColNames <- intersect(colnames(transcriptomicsInputData), colnames(lipidomicsInputData))
transcriptomicsInputData <- transcriptomicsInputData[
    ,c(colnames(transcriptomicsInputData)[1],commonColNames)]
lipidomicsInputData <- lipidomicsInputData[
    ,c(colnames(lipidomicsInputData)[1],commonColNames)]
# rownames(transcriptomicsInputData)
# rownames(lipidomicsInputData)

# UWAGI do vignetty: Manipulacja yCutoff i xCutoff bardzo ważna!!!
# UWAGI do vignetty: Using cutoff can handle errors like:
#       i)   NaNs produced
#       ii)  singular matrix 'a' in solve
#       iii) imaginary parts discarded in coercion
ccaResultsNarrow1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.55, yCutoff = 0.75)
ccaResultsNarrow2 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.5, yCutoff = 0.6)
# IMPORTANT: Compare and analyze two results with different cutoff!
par(mfrow=c(1,2))
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow1)
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow2)

# Same with Expand set (it will take time for operations)
# UWAGI do vignetty: "CCA is not possible - XDataFrame do not have specified names in xNamesVector"
#       Sprawdz dane wejściowe do funkcji, print na argumentach
ccaResultsExpand1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.5, yCutoff = 0.7)
ccaResultsExpand2 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.4, yCutoff = 0.5)

par(mfrow=c(1,2))
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsExpand1)
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsExpand2)

# Only with Recatome decorations.
ccaResultsEnsemble <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromEnsembleFunctionalInteractions$genesSymbolsFromEnsemble,
    yNamesVector = ontology2GenesSymboleFromEnsembleFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.5, yCutoff = 0.75)
ccaResultsUniprot <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromUniProtFunctionalInteractions$genesSymbolsFromUniProt,
    yNamesVector = ontology2GenesSymboleFromUniProtFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.5, yCutoff = 0.75)

par(mfrow=c(1,2))
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsEnsemble)
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsUniprot)

# SHOULD BE: Comparison of raw data without any extentions, what are the results.
# TODO

# PARAGRAPH: PLS - Partial Least Squares Regression
PLSResults <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)
# UWAGI do vignetty: Error in plot.new() : figure margins too large.
#       Powiększyć okno Plots. :)
# UWAGI do kodu: Plotowanie na jednym obszarze?
OmicsON::plotRmsepForPLS(PLSResults)
OmicsON::plotRegression(PLSResults, ncompValue = 10)






# PARAGRAPH:


# PARAGRAPH:


# PARAGRAPH:






