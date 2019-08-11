CcaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 1, yCutoff = 1)

for (i in seq(from = 0.4,to=0.6, by=0.1)) {
    for (j in seq(from = 0.4,to=0.6, by=0.1)) {
        print("********************************************")
        print(paste("i = ", i, ",", "j = ", j))
        print("********************************************")
        CcaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
            xNamesVector = xNamesVector,
            yNamesVector = yNamesVector,
            XDataFrame = XDataFrame,
            YDataFrame = YDataFrame,
            xCutoff = i, yCutoff = j)
    }
}

for (i in seq(from = 0.2,to=0.4, by=0.1)) {
    for (j in seq(from = 0.2,to=0.4, by=0.1)) {
        print("********************************************")
        print(paste("i = ", i, ",", "j = ", j))
        print("********************************************")
        CcaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
            xNamesVector = xNamesVector,
            yNamesVector = yNamesVector,
            XDataFrame = XDataFrame,
            YDataFrame = YDataFrame,
            xCutoff = i, yCutoff = j)
    }
}# PCA

for (i in seq(from = 0.7,to=1.0, by=0.1)) {
    for (j in seq(from = 0.7,to=1.0, by=0.1)) {
        print("********************************************")
        print(paste("i = ", i, ",", "j = ", j))
        print("********************************************")
        CcaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
            xNamesVector = xNamesVector,
            yNamesVector = yNamesVector,
            XDataFrame = XDataFrame,
            YDataFrame = YDataFrame,
            xCutoff = i, yCutoff = j)
    }
}

i <- 0.6
for (j in seq(from = 0.5, to=1.0, by=0.1)) {
    print("********************************************")
    print(paste("i = ", i, ",", "j = ", j))
    print("********************************************")
    CcaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
        xNamesVector = xNamesVector,
        yNamesVector = yNamesVector,
        XDataFrame = XDataFrame,
        YDataFrame = YDataFrame,
        xCutoff = i, yCutoff = j)
}
