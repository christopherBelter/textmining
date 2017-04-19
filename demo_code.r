setwd("C:/Documents/Text mining")
source("pubmedXML.R")
mystopwords <- scan("stopwords.txt", what="varchar", skip=1)
library(reutils)
mquery <- "huntington disease[mh] AND 2013:2015[dp]"
pmids <- esearch(mquery, db = "pubmed", usehistory = TRUE)
pmids
rxml <- efetch(pmids, retmode = "xml", outfile = "huntingtonPubs_raw.xml")
cxml <- clean_api_xml(rxml, outfile = "huntingtonPubs_clean.xml")
theData <- extract_xml(cxml)
rm(pmids)
rm(cxml)
library(plyr)
mesh <- count(unlist(strsplit(theData$meshHeadings, "|", fixed = TRUE)))
mesh <- mesh[order(mesh$freq),]
tail(mesh, 25)
docs <- data.frame(theData$abstract)
pmid <- as.vector(theData$pmid)
library(tm)
library(slam)
corpDocs <- Corpus(DataframeSource(docs))
corpDocs <- tm_map(corpDocs, removePunctuation)
corpDocs <- tm_map(corpDocs, content_transformer(tolower))
corpDocs <- tm_map(corpDocs, removeWords, stopwords("English"))
corpDocs <- tm_map(corpDocs, removeWords, mystopwords)
corpDocs <- tm_map(corpDocs, stemDocument)
corpDocs <- tm_map(corpDocs, stripWhitespace)
dtm <- DocumentTermMatrix(corpDocs)
rownames(dtm) <- pmid
dtm <- dtm[row_sums(dtm) > 0,]
dtm
dtm <- removeSparseTerms(dtm, 0.99)
dtm
findFreqTerms(dtm, 400)
tDocs <- findMostFreqTerms(dtm)
tDocs[[5]]
tDocs[[52]]
tDocs$`26767207`
theTerms <- col_sums(dtm)
theTerms <- sort(theTerms)
tail(theTerms, 50)
assoc1 <- findAssocs(dtm, "cag", 0.2)
assoc1
assoc2 <- findAssocs(dtm, "protein", 0.2)
assoc2
assoc3 <- findAssocs(dtm, "gene", 0.2)
assoc3
assoc4 <- findAssocs(dtm, c("htt", "huntingtin"), 0.2)
assoc4
library(topicmodels)
seed <- list(79, 524, 1291, 706, 1044)
ldaOut <- LDA(dtm, 10, method = "Gibbs", control = list(nstart = 5, seed = seed, best = TRUE, burnin = 4000, iter = 2000, thin = 500))
docTopics <- as.matrix(topics(ldaOut))
write.csv(docTopics, file = "docstotopics.csv")
topicTerms <- as.matrix(terms(ldaOut, 15))
write.csv(topicTerms, file = "topicsToTerms.csv")
topicTerms
topicProbs <- as.data.frame(ldaOut@gamma)
write.csv(topicProbs, file = "topicProbabilities.csv")
topicList <- topics(ldaOut, threshold = 0.15)
topicList <- sapply(topicList, paste0, collapse = "|")
newData <- merge(theData, docTopics, by.x = "pmid", by.y = "row.names", all.x = TRUE)
newData <- merge(newData, topicList, by.x = "pmid", by.y = "row.names", all.x = TRUE)
write.csv(newData, file = "clustered_data.csv")