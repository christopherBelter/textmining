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

Note: if you get an error message after the stemDocument command, it probably means you don't have the snowballC package installed, which the tm package uses to stem words. Just install the snowballC package and try the command again and it'll probably work. 

Finally, we can create the document-term matrix, set the rownames of the matrix to be the document PMIDs (so we can map the results back onto the original document set later on) and remove documents without abstracts to avoid future errors.

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

Another common task in text mining is to group similar documents into topics. There are lots of ways to do this, but in this example I'll use Latent Dirichlet Analysis, also referred to as LDA or topic modeling. Essentially, LDA works by first creating a set of term vocabularies (based on term co-occurrences in documents) for each of a prespecified number of topics, and then using those vocabularies to assign a probability that each document belongs to each topic. These probabilities allow documents to be assigned to multiple topics, which we'll do a bit later on.

The first step is to load the package and to set the seed for the LDA algorithm. The seed is a list of document numbers that the algorithm will use to start it's sampling process. The seed should be different for each document set, so choose 5 random numbers that are all less than the total number of documents in your data set. Since the document set here includes around 1300 papers, I'll set them all to less than that. 

    library(topicmodels)
    seed <- list(79, 524, 1291, 706, 1044)

With that done, we can now run the actual algorithm. Be careful with this step, because depending on your computer and the number of documents you have, it can take anywhere from 5 minutes to 48 hours to run. 

    ldaOut <- LDA(dtm, 10, method = "Gibbs", control = list(nstart = 5, seed = seed, best = TRUE, burnin = 4000, iter = 2000, thin = 500))

There's a lot in that command, so let me point out the relevant parts. The LDA() function calls the algorithm and the arguments set the parameters you want to use. In this case, I'm running it on the document term matrix I created in the tm package, I'm telling it to look for 10 topics in my document set, to use the Gibbs sampling method, and a number of other arguments in the 'control' argument. Most of this should stay the same from analysis to analysis, but you might want to change the number of topics, especially for larger data sets to obtain more or less granular topics. 

Once the algorithm is finished, there's a lot you can do with this ldaOut object. First, we'll extract the primary topic number for each document and save that to a .csv file.

    docTopics <- as.matrix(topics(ldaOut))
    write.csv(docTopics, file = "docstotopics.csv")

Next, we'll extract the top 15 terms that appear in each topic. That is, we'll create lists of the 15 terms that most frequently appear in each of the 10 topics the algorithm identified. We'll also save that list to a .csv file and print it to the console to have a look. These terms are essential for understanding what each topic is actually about and subsequently naming it.

    topicTerms <- as.matrix(terms(ldaOut, 15))
    write.csv(topicTerms, file = "topicsToTerms.csv")
    topicTerms

Next, we'll extract the full set of probabilities that each document belongs to each topic and save it to a .csv file.

    topicProbs <- as.data.frame(ldaOut@gamma)
    write.csv(topicProbs, file = "topicProbabilities.csv")

One thing that we can do with these probabilities is to generate a list of topics to which each document has a high probability of belonging. In this case, we'll extract a list of the topics to which each document has at least a 15% probability of belonging, and then paste the list for each document together to create a vector of term lists.

    topicList <- topics(ldaOut, threshold = 0.15)
    topicList <- sapply(topicList, paste0, collapse = "|")

Finally, we'll merge the both the primary topic and the topic probability list back into our original publications data frame and write the clustered data to a .csv file for later use.

    newData <- merge(theData, docTopics, by.x = "pmid", by.y = "row.names", all.x = TRUE)
    newData <- merge(newData, topicList, by.x = "pmid", by.y = "row.names", all.x = TRUE)
    write.csv(newData, file = "clustered_data.csv")

Obviously there are lots of other things that we could do with this data set, but this gives you an idea of what's possible and some of the basic code to get you there. 
