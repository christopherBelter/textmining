# Text mining in R

## About this repo

This repo provides files and sample code to do text mining in R. It includes code to download articles from PubMed using the reutils package, extract the resulting XML into a data frame, and then perform simple text mining and document clustering tasks using the tm and topicmodels packages. In addition to these packages, you will also need the XML, plyr, slam, and snowballC packages, as well as the pubmedXML.R and stopwords files included in this repo. Step by step information about what the code is doing is below, and the full code without comments is provided in the demo_code.r file in this repo. 

## How the code works

### Download the records

First, set your working directory and load the necessary helper files and packages.

    setwd("C:/Documents/Text mining")
    source("pubmedXML.R")
    mystopwords <- scan("stopwords.txt", what="varchar", skip=1)
    library(reutils)

Then save the PubMed search query you want to run. In this case I'll search for publications on Huntington's Disease published over a three-year time period.

    mquery <- "huntington disease[mh] AND 2013:2015[dp]"

Next, run the search against the PubMed API and save the pmids of the search results to NCBI's history server so you can download them later.

    pmids <- esearch(mquery, db = "pubmed", usehistory = TRUE)

Now you can get some basic information about the search results by

    pmids

Then download the full records of the search results in XML format and save them to a text file.

    rxml <- efetch(pmids, retmode = "xml", outfile = "huntingtonPubs_raw.xml")

Since the search returned more than 500 results, the reutils package downloaded them in batches of 500 and saved all of the batches to the outfile specified in the efetch() function. As a result, the outfile actually consists of three separate XML files pasted together, so you cannot parse it as it is. So, use the helper function clean_api_xml() in the pubmedXML.R file to remove the extra headers and footers.

    cxml <- clean_api_xml(rxml, outfile = "huntingtonPubs_clean.xml")

Now you can use the other helper function in the pubmedXML.R file to extract the XML data into a data frame in R.

    theData <- extract_xml(cxml)

Finally, remove the pmids object from the PubMed history server and the cxml object to free up some additional memory.

    rm(pmids)
    rm(cxml)

### Docuemnt preprocessing

Now that you have the documents in a data frame in R, you can do some basic count analyses using various functions in R. For example, you can use the plyr package to count and display the most frequently occurring MeSH headings in the publication set.

    library(plyr)
    mesh <- count(unlist(strsplit(theData$meshHeadings, "|", fixed = TRUE)))
    mesh <- mesh[order(mesh$freq),]
    tail(mesh, 25)

To do more detailed analyses, we'll create a document-term matrix from the publication set in which rows are documents and columns are terms that appear in their abstracts. First, we'll isolate the columns we're interested in using, namely the publications' abstracts and PMIDs.

    docs <- data.frame(theData$abstract)
    pmid <- as.vector(theData$pmid)

Then we'll load the relevant packages for doing text mining and working with the resulting sparse document-term matrix. 

    library(tm)
    library(slam)

The first step in working with the tm package is to create a document corpus, in this case using the docs object we just created.

    corpDocs <- Corpus(DataframeSource(docs))

Now before we can do anything meaningful with the publication abstracts, we need to do some preprocessing on those abstracts. We'll use the tm_map() function to apply a series of functions to each document in the corpus. The next lines remove punctuation, transform all of the characters to lowercase, remove stopwords like 'and', 'the', and 'of' in two steps (the first uses the base set of stopwords and the second uses a list of custom stopwords provided in this repo), stem the terms in each document to remove word endings like '-s' and '-ing', and finally strip any extra whitespace created in the previous operations.

    corpDocs <- tm_map(corpDocs, removePunctuation)
    corpDocs <- tm_map(corpDocs, content_transformer(tolower))
    corpDocs <- tm_map(corpDocs, removeWords, stopwords("English"))
    corpDocs <- tm_map(corpDocs, removeWords, mystopwords)
    corpDocs <- tm_map(corpDocs, stemDocument)
    corpDocs <- tm_map(corpDocs, stripWhitespace)

Finally, we can create the document-term matrix, set the rownames of the matrix to be the document PMIDs (so we can map the results back onto the original document set later on) and remove documents without abstracts to avoid errors later on.

    dtm <- DocumentTermMatrix(corpDocs)
    rownames(dtm) <- pmid
    dtm <- dtm[row_sums(dtm) > 0,]

As before, we can get some basic information about the document term matrix we've created by simply doing

    dtm

You'll note the high term sparsity in the matrix. This essentially means we have lots of terms that are only used in one document in the set. We can safely remove them to speed up our functions without losing any meaningful information, so we'll go ahead and do so and then see what the dtm looks like afterwards.

    dtm <- removeSparseTerms(dtm, 0.99)
    dtm

### Finding frequently occurring terms and term correlations

One of the most common tasks in text mining is to find the most frequently occurring terms in the document set. The tm package provides a couple of ways of doing this, but we can also use the slam package to do others. First, we can find the terms that appear at least x times in the entire data set using the findFreqTerms() function. We'll find the terms that appear at least 400 times in the data set.

    findFreqTerms(dtm, 400)

Another thing you can do is find the most frequently occurring terms in each individual document in the data set. We'll create a new object, tDcos, to store this information. 

    tDocs <- findMostFreqTerms(dtm)

To access the terms appearing in any given document, we can use the standard R syntax for accessing lists

    tDocs[[5]]
    tDocs[[52]]

or specify a document by it's PMID.

    tDocs$`26767207`

Finally, we can also create a list of terms that appear most frequently in the entire document set and actually get the number of times each term appears in all documents. You could probably use the base colSums() function, but I happen to use the col_sums() function from the slam package instead. Here, we'll create a list of all terms that appear in the dtm, sort them by occurrence frequency, and then print out the top 50. 

    theTerms <- col_sums(dtm)
    theTerms <- sort(theTerms)
    tail(theTerms, 50)

Another common task in text mining is to find associations between terms. That is, we can generate a list of terms that frequently co-occur with a specified term or set of terms. The function in tm is findAssocs(). In this example, we'll find term associations with the terms 'cag' (actually the genetic marker for Huntington's disease), 'protein', and 'gene' with a correlation of 0.2 or higher. 

    assoc1 <- findAssocs(dtm, "cag", 0.2)
    assoc1
    assoc2 <- findAssocs(dtm, "protein", 0.2)
    assoc2
    assoc3 <- findAssocs(dtm, "gene", 0.2)
    assoc3

We can also create a list of associated terms for a list of terms using R's c() function. In this example, we'll look for terms associated with either 'htt' or 'huntingtin', which are variant ways of referring to the specific gene associated with Huntington's Disease.

    assoc4 <- findAssocs(dtm, c("htt", "huntingtin"), 0.2)
    assoc4

### Cluster documents by topic
x

    library(topicmodels)
    seed <- list(79, 524, 1291, 706, 1044)

x

    ldaOut <- LDA(dtm, 10, method = "Gibbs", control = list(nstart = 5, seed = seed, best = TRUE, burnin = 4000, iter = 2000, thin = 500))

x

    docTopics <- as.matrix(topics(ldaOut))
    write.csv(docTopics, file = "docstotopics.csv")

x

    topicTerms <- as.matrix(terms(ldaOut, 15))
    write.csv(topicTerms, file = "topicsToTerms.csv")
    topicTerms

x

    topicProbs <- as.data.frame(ldaOut@gamma)
    write.csv(topicProbs, file = "topicProbabilities.csv")


x

    topicList <- topics(ldaOut, threshold = 0.15)
    topicList <- sapply(topicList, paste0, collapse = "|")

x

    newData <- merge(theData, docTopics, by.x = "pmid", by.y = "row.names", all.x = TRUE)
    newData <- merge(newData, topicList, by.x = "pmid", by.y = "row.names", all.x = TRUE)

x

    write.csv(newData, file = "clustered_data.csv")

Obviously there are lots of other things that we could do with this data set, but this gives you an idea of what's possible and some of the basic code to get you there. 
