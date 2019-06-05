## ----global_options, include=TRUE, echo=FALSE----------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, error = TRUE)

## ---- results='asis', echo=TRUE------------------------------------------
    OmicsON::setUpReactomeMapping(ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt", 
                                  Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt", 
                                  UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")

## ---- results = 'asis', echo=TRUE----------------------------------------
pathToFileWithLipidomicsData <- system.file(package="OmicsON", 
                                            "extdata", "nm-lipidomics-min.txt")
lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(lipidomicsInputData[1:4, 1:7], caption = "Lipidomisc data")

## ---- results = 'asis', echo=TRUE----------------------------------------
pathToFileWithTranscriptomicsData <- system.file(package="OmicsON", 
                                                 "extdata", "nm-transcriptomics-min.txt")
transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)

## ---- results = 'asis'---------------------------------------------------
knitr::kable(transcriptomicsInputData[1:10, 1:7], caption = "Transcriptomics data")

## ---- echo=TRUE, results='hide'------------------------------------------
decoratedByReactome <- OmicsON::decorateByReactomeData(
    chebiMoleculesDf = lipidomicsInputData, 
    chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')

## ---- echo=TRUE, results='hide'------------------------------------------
decoratedByStringBaseOnEnsembleIds <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactome, 
    listOfEnsembleIdColumnName = 'ensembleIds')
decoratedByStringBaseOnUniProtIds <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactome, 
    listOfEnsembleIdColumnName = 'uniProtIds')

## ---- echo=TRUE, results='asis'------------------------------------------
as.character(decoratedByStringBaseOnEnsembleIds[4, "root"])
as.character(decoratedByStringBaseOnUniProtIds[4, "root"])

## ---- echo=TRUE, results='asis'------------------------------------------
decoratedByStringBaseOnEnsembleIds[4, "ensembleIds"][[1]]
decoratedByStringBaseOnUniProtIds[4, "uniProtIds"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
decoratedByStringBaseOnEnsembleIds[2, "stringGenesSymbolsNarrow"][[1]]
decoratedByStringBaseOnUniProtIds[2, "stringGenesSymbolsNarrow"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
ontology2EnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactome, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'ensembleIds')

ontology2UniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactome, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'uniProtIds')

ontology2GenesSymboleFromEnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactome, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'genesSymbolsFromEnsemble')

ontology2GenesSymboleFromUniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactome, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'genesSymbolsFromUniProt')

ontology2GenesSymboleFromStringExpandFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByStringBaseOnUniProtIds, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'stringGenesSymbolsExpand')

ontology2GenesSymboleFromStringNarrowFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByStringBaseOnUniProtIds, 
    singleIdColumnName = 'ontologyId', idsListColumnName = 'stringGenesSymbolsNarrow')

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ontology2EnsembleFunctionalInteractions, 6))
knitr::kable(head(ontology2UniProtFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromEnsembleFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromUniProtFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromStringExpandFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromStringNarrowFunctionalInteractions, 6))

## ---- echo=TRUE, results='asis'------------------------------------------
ccaResultsNarrow1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 1, yCutoff = 1)
ccaResultsNarrow2 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, xCutoff = 0.5, yCutoff = 1)

## ---- echo=TRUE, results='asis', error=TRUE------------------------------
ccaResultsExpand1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector =  ontology2GenesSymboleFromStringExpandFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)

## ---- echo=TRUE, results='asis'------------------------------------------
ccaResultsExpand1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector =  ontology2GenesSymboleFromStringExpandFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, 
    xCutoff = 0.7)

## ---- fig.show='hold', fig.width=7, fig.height=5, echo=TRUE, fig.align='center'----
par(mfrow=c(1,2))
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow1)
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow2)
par(mfrow=c(1,1))
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsExpand1)

## ---- echo=TRUE, results='hide'------------------------------------------
PLSResults <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
    yNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$root,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)

## ---- fig.show='hold', fig.width=7, fig.height=5, echo=TRUE--------------
OmicsON::plotRmsepForPLS(PLSResults)

## ---- fig.show='hold', fig.width=7, fig.height=5, echo=TRUE--------------
OmicsON::plotRegression(PLSResults, ncompValue = 10)

