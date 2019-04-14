## ----global_options, include=TRUE----------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, error = TRUE)

## ---- results = 'asis'---------------------------------------------------
    OmicsON::setUpReactomeMapping(ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt", 
                                  Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt", 
                                  UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")

## ---- results = 'asis'---------------------------------------------------
    pathToFileWithLipidomicsData <- system.file(package="OmicsON", "extdata", "nm-lipidomics.txt")
    lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
    lipidomicsInputDf <- head(lipidomicsInputData, 6)
    knitr::kable(lipidomicsInputDf[1:7], caption = "Lipidomisc data")
    
    pathToFileWithTranscriptomicsData <- system.file(package="OmicsON", "extdata", "nm-transcriptomics.txt")
    transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)
    transcriptomicsInputDf <- head(transcriptomicsInputData, 6)
    knitr::kable(transcriptomicsInputDf[1:7], caption = "Transcriptomics data")

## ---- echo=TRUE, results='asis'------------------------------------------
    dataToVignetteEvaluation <- lipidomicsInputData[lipidomicsInputData[,"ChEBI"] %in% c("CHEBI:28875","CHEBI:73705", "CHEBI:35465"),]
    decoratedByReactome <- OmicsON::decorateByReactomeData(chebiMoleculesDf = dataToVignetteEvaluation,
                                chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')

## ---- echo=TRUE, results='hide'------------------------------------------
    decoratedByReactome

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnEnsembleIds <- OmicsON::decorateByStringDbData(
        chebiIdsToReactomePathways = decoratedByReactome, listOfEnsembleIdColumnName = 'ensembleIds')
    decoratedByStringBaseOnUniProtIds <- OmicsON::decorateByStringDbData(
        chebiIdsToReactomePathways = decoratedByReactome, listOfEnsembleIdColumnName = 'uniProtIds')

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnUniProtIds

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnEnsembleIds[1, "root"]

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnEnsembleIds[1, "ensembleIds"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnEnsembleIds[1, "stringGenesSymbolsNarrow"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnUniProtIds[1, "root"]

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnUniProtIds[1, "uniProtIds"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
    decoratedByStringBaseOnUniProtIds[1, "stringGenesSymbolsNarrow"][[1]]

## ---- echo=TRUE, results='asis'------------------------------------------
    ontology2EnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByReactome, 
                                                                             singleIdColumnName = 'ontologyId', 
                                                                             idsListColumnName = 'ensembleIds')
    ontology2UniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByReactome, 
                                                                             singleIdColumnName = 'ontologyId', 
                                                                             idsListColumnName = 'uniProtIds')
    ontology2GenesSymboleFromEnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByReactome, 
                                                                             singleIdColumnName = 'ontologyId', 
                                                                             idsListColumnName = 'genesSymbolsFromEnsemble')
    ontology2GenesSymboleFromUniProtFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByReactome, 
                                                                             singleIdColumnName = 'ontologyId', 
                                                                             idsListColumnName = 'genesSymbolsFromUniProt')
    ontology2GenesSymboleFromStringExpandFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByStringBaseOnUniProtIds, 
                                                   singleIdColumnName = 'ontologyId', 
                                                   idsListColumnName = 'stringGenesSymbolsExpand')
    ontology2GenesSymboleFromStringNarrowFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(decoratedByStringBaseOnUniProtIds, 
                                                   singleIdColumnName = 'ontologyId', 
                                                   idsListColumnName = 'stringGenesSymbolsNarrow')

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ontology2EnsembleFunctionalInteractions, 6))
knitr::kable(head(ontology2UniProtFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromEnsembleFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromUniProtFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromStringExpandFunctionalInteractions, 6))
knitr::kable(head(ontology2GenesSymboleFromStringNarrowFunctionalInteractions, 6))

## ---- echo=TRUE, results='asis'------------------------------------------
ccaResultsNarrow.73705 <- OmicsON::makeCanonicalCorrelationAnalysis(
        xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
        yNamesVector = c("CHEBI:73705"),
        XDataFrame = transcriptomicsInputData,
        YDataFrame = lipidomicsInputData)

ccaResultsNarrow.73705.28875 <- OmicsON::makeCanonicalCorrelationAnalysis(
        xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
        yNamesVector = c("CHEBI:73705", "CHEBI:28875"),
        XDataFrame = transcriptomicsInputData,
        YDataFrame = lipidomicsInputData)

## ---- echo=TRUE, results='asis'------------------------------------------
ccaResultsExpand.73705 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector =  c("CHEBI:73705"),
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)

## ---- echo=TRUE, results='asis'------------------------------------------
ccaResultsExpand.73705 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = ontology2GenesSymboleFromStringExpandFunctionalInteractions$stringGenesSymbolsExpand,
    yNamesVector =  c("CHEBI:73705"),
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData, 
    xCutoff = 0.7)

## ---- echo=TRUE, fig.show='hold', fig.width=6, fig.height=6--------------
    OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsExpand.73705)
    OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow.73705)
    OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow.73705.28875)

## ---- echo=TRUE, results='hide'------------------------------------------
    PLSResults <- OmicsON::makePartialLeastSquaresRegression(
        xNamesVector = ontology2GenesSymboleFromStringNarrowFunctionalInteractions$stringGenesSymbolsNarrow,
        yNamesVector = c("CHEBI:73705","CHEBI:28875"),
        XDataFrame = transcriptomicsInputData,
        YDataFrame = lipidomicsInputData)


## ---- echo=TRUE, fig.show='hold', fig.width=6, fig.height=6--------------
    OmicsON::plotRmsepForPLS(PLSResults)

## ---- fig.show='hold', fig.width=6, fig.height=6-------------------------
    OmicsON::plotRegression(PLSResults, ncompValue = 10)

