# Text mining in R

## About this repo

This repo provides files and sample code to do text mining in R. It includes code to download articles from PubMed using the rentrez package, extract the resulting XML into a data frame, and then perform simple text mining and document clustering tasks using the tm and text2vec packages. In addition to these packages, you will also need the XML, plyr, and snowballC packages, as well as the pubmedXML.R and stopwords files included in this repo. Step by step information about what the code is doing is below, and the full code without comments is provided in the demo_code.r file in this repo. 

Note that this is an updated version of my previous text mining workflow. After working with the previous workflow for several years, I have switched to a different set of packages to do the publication retrieval and document clustering tasks. The new version results in a more reliable XML file to work with and a substantial improvement in the amount of time these tasks take to run. The previous workflow is still available in this repo as README_old.md. 

## How the code works

### Download the records

First, set your working directory and load the necessary helper files and packages.

    setwd("C:/Documents/Text mining")
    source("pubmedXML.R")
	mystopwords <- scan("stopwords.txt", what="varchar", skip=1)
    library(rentrez)

Then save the PubMed search query you want to run. In this case I'll search for all publications in PubMed on Zika virus.

    mquery <- "zika virus[mh]"

Next, run the search against the PubMed API and save the pmids of the search results to NCBI's history server so you can download them later.

    zika_ids <- entrez_search(db = "pubmed", term = mquery, use_history = TRUE)

The resulting object, zika_ids, gives you all of the information you need to download the results. For example, you can see how many records were returned by the search by doing

    zika_ids$count

You can then download the full records of the search results in XML format and save them to a text file.

    zika <- entrez_fetch(db = "pubmed", web_history = zika_ids$web_history, rettype = "xml")

The entrez_fetch() function will automatically download the first 10,000 search results for any query. If you want to retrieve more than 10,000 records, you can run the function multiple times using the optional "retstart" argument to request subsequent sets of 10,000 records.  

Once the download process is complete, you can use the extract_xml() helper function in the pubmedXML.R file to extract the XML data into a data frame in R.

    pubs <- extract_xml(zika)

Finally, remove the search results from the PubMed history server, save the XML file, and remove the XML object from your session to free up some additional memory.

    rm(zika_ids)
    readr::write_file(zika, "zikaPubs.xml")
    rm(zika)

### Docuemnt preprocessing

Now that you have the documents in a data frame in R, you can do some basic count analyses using various functions in R. For example, we can use the plyr package to count and display the most frequently occurring MeSH headings in the publication set.

    library(plyr)
    mesh <- count(unlist(strsplit(pubs$meshHeadings, "|", fixed = TRUE)))
    mesh <- mesh[order(mesh$freq),]
    tail(mesh, 15)

To do more detailed analyses on the text of these documents, we need to do some preprocessing to prepare the documents for analysis. First, we'll create a data frame with the document ID (we'll use the PMID for this) and a combined title and abstract field that contains all of the relevant text for us to work with.

    abstracts <- data.frame(pmid = pubs$pmid, text = paste(pubs$articletitle, pubs$abstract, sep = ". "), stringsAsFactors = FALSE)

Then we'll load the tm package. 

    library(tm)

Now before we can do anything meaningful with these documents, we need to do some preprocessing on the text using a series of functions from the tm package. The next lines remove punctuation, transform all of the characters to lowercase, remove stopwords like 'and', 'the', and 'of' in two steps (the first uses the base set of stopwords and the second uses a list of custom stopwords provided in this repo), stem the terms in each document to remove word endings like '-s' and '-ing', and finally strip any extra whitespace created in the previous operations.

    abstracts$text <- removePunctuation(abstracts$text)
    abstracts$text <- tolower(abstracts$text)
    abstracts$text <- removeWords(abstracts$text, stopwords("English"))
    abstracts$text <- removeWords(abstracts$text, mystopwords)
    abstracts$text <- stemDocument(abstracts$text)
    abstracts$text <- stripWhitespace(abstracts$text)

Note: if you get an error message after the stemDocument command, it probably means you don't have the snowballC package installed, which the tm package uses to stem words. Just install the snowballC package and try the command again and it'll probably work. 

The final proprocessing steps are to load the text2vec package and create an iterator function to work with our specific document set. 

    library(text2vec)
    it <- itoken(abstracts$text, pmid = abstracts$pmid) 

### Finding frequently occurring terms and term associations

One of the most common tasks in text mining is to find the most frequently occurring terms in the document set. The text2vec package does this with the create_vocabulary() function, which essentially creates a list of the unique terms found in all documents in the data set, along with the number of times each term appears in the data set and the number of documents in which each term appears.

    v <- create_vocabulary(it)
	tail(v, 10)

You can also prune the vocabulary to remove terms that don't appear in at least 0.1% of the documents in the data set. 

    v <- prune_vocabulary(v, doc_proportion_min = 0.001)
	head(v, 10)

Or, you could prune the vocabulay to remove terms that do not appear in at least 4 documents. 

    v <- prune_vocabulary(v, doc_count_min = 4)
	head(v, 10)

Another common task in text mining is to find associations between terms. That is, we want to find terms that frequently occur together or near each other in our documents. To do these kinds of analyses, we'll create a term co-occurrence matrix (or tcm) that calculates the number of times two terms appear within x terms of each other in our document set. We first create a vectorizer to iterate over our documents and then create the tcm itself. 

    vectorizer <- vocab_vectorizer(v)
    tcm <- create_tcm(it, vectorizer, skip_grams_window = 5)

The skip_grams_window argument specifies the maximum distance between terms, or the maximum number of terms that can separate any two terms. Setting it to 5 means that we want to count co-occurrences of terms that appear within 5 terms of each other in our document set. Setting it to 3 also works well. 

We can now get some information about the tcm using some basic R functions. 

    dim(tcm)
	str(tcm)

You'll notice that this doesn't look like a traditional matrix, because the text2vec package stores matrices using structures from the Matrix package to save memory. That does help speed up the processing time, but it means we can't use traditional R subsetting methods to work with it. So we have to be creative.

One thing we might want to do is find the most frequently co-occurring terms in our document set. With a traditional matrix, we'd just use [] to subset it, but with the different structure, we need to do something different. 

    highC <- which(tcm@x > 50)
    t1 <- tcm@i[highC]
    t2 <- tcm@j[highC]
    assoc <- data.frame(term1 = tcm@Dimnames[[1]][t1 + 1], term2 = tcm@Dimnames[[2]][t2 + 1], weight = tcm@x[highC])

Essentially, we first create a list of the matrix index values that are greater than 50, pull the row (here, i) and column (j) numbers for those values, and then create a data frame by pulling the correct terms and values from the matrix. As before, we can have a look at the results.

    head(assoc, 10) 

We can also pull associations with specific terms that we might be interested in. Since these documents are about Zika, we might be interested in terms frequently co-occurring with the term 'microcephaly.' 

    microcephali <- c(tcm[grep("microcephali\\b", tcm@Dimnames[[1]]),], tcm[,grep("microcephali\\b", tcm@Dimnames[[2]])])
    microcephali <- microcephali[microcephali > 0]
    microcephali <- sort(microcephali)
    tail(microcephali, 10)

First we found all rows and columns in the matrix where the term 'microcephaly' occurs, removed all of the 0 values, and then sorted the results to find the most frequent terms. We can do the same thing with the term 'mosquito'

    mosquito <- c(tcm[grep("mosquito\\b", tcm@Dimnames[[1]]),], tcm[,grep("mosquito\\b", tcm@Dimnames[[2]])])
    mosquito <- mosquito[mosquito > 0]
    mosquito <- sort(mosquito)
    tail(mosquito, 10)

or any other term we might be interested in. 

### Cluster documents by topic

Another common task in text mining is to group similar documents into topics. There are lots of ways to do this, but in this example I'll use Latent Dirichlet Analysis, also referred to as LDA or topic modeling. Essentially, LDA works by first creating a set of term vocabularies (based on term co-occurrences in documents) for each of a prespecified number of topics, and then using those vocabularies to assign a probability that each document belongs to each topic. These probabilities allow documents to be assigned to multiple topics, which we'll do a bit later on.

The first step is to create a document-term matrix where rows are documents in the data set and columns are terms that appear in all documents. The values are the number of times each term appears in each document.

    dtm <- create_dtm(it, vectorizer, type = "dgTMatrix")
	str(dtm)

Once again, the structure of the matrix isn't what we're used to in R, but it is what's required for the text2vec functions to work. The function syntax in the text2vec package is also unlike anything else I've seen in R, so fair warning. 

Now that we have the dtm created, we can set the parameters of the LDA algorithm run by creating an LDA model that we'll fit to our data set. 

    lda_model <- LDA$new(n_topics = 10, doc_topic_prior = 0.1, topic_word_prior = 0.01)
	
The most important thing about that line of code is the n_topics argument, which specifies the number of topics you want the algorithm to find. The other arguments control how the model runs, but in practice I never change them. 

With that done, we can now run the actual algorithm on our data set. 

doc_topic_dist <- lda_model$fit_transform(x = dtm, n_iter = 2000, convergence_tol = 0.000001, n_check_convergence = 100)

Again, the arguments are the important things to point out here. The x argument is what we're running the algorithm on (i.e. our data) and the n_iter argument specifies the number of sampling iterations to run. In practice I set this to 2000, but the algorithm rarely runs that many iterations. Instead, it will stop early if the additional iterations don't improve the model's performance. 

You can control that early stopping with the remaining arguments. The convergence_tol argument sets a threshold at which the algorithm will stop. I typically set this extremely low because in practice I find that more iterations leads to better results. The final argument, n_check_convergence, is how frequently you want the algorithm to check whether the additional iterations increase the model's performance. 

Once the algorithm is finished, the result is a document-topic matrix giving you the probability that each document belongs to each topic. We can get an idea of what these topics are by extracting and printing out the most frequently occurring terms in each topic. These terms are essential for understanding what each topic is actually about and subsequently naming it.

    modelTerms <- lda_model$get_top_words(n = 15, lambda = 1.0)
    modelTerms

The n argument specifies the number of terms we want to get for each topic. The lambda parameter is an interesting one that I honestly don't fully understand. My sense is that a value of 1.0 means that you want the most frequent terms in each topic overall, whereas setting it to values below 1.0 means that you want terms that are more specific to the documents in that topic. That is, to words that rarely appear in other topics. 

In practice, I usually run this twice: once with lanbda = 1.0 and again with lambda = 0.4 to try and identify terms that are unique to each topic. So:

    terms2 <- lda_model$get_top_words(n = 15, lambda = 0.4)
    terms2
	
The two term lists are usually similar, but there are definitely differences between them that can be helpful in understanding what the topics are actually about. 

Next, we want to assign documents to topics using the probability matrix (doc_topic_dist) that we just created. Usually I assign each document to the topic that it has the highest probability of belonging to. 

    docTopics <- data.frame(pmid = rownames(doc_topic_dist), topic = apply(doc_topic_dist, 1, which.max))
	
But we can also generate a list of topics to which each document has a high probability of belonging. Since, in the LDA approach, documents can belong to multiple topics, we can represent these topic overlaps by extracting all of the topics to which a document has a certain probability of belonging. In this case, we'll extract a list of the topics to which each document has at least a 15% probability of belonging.

    tlist <- data.frame(pmid = rownames(doc_topic_dist), tList = sapply((sapply(1:nrow(doc_topic_dist), function(x) which(doc_topic_dist[x,] > 0.15))), paste, collapse = ";"))

Finally, we'll merge the both the primary topic and the topic probability list back into our original publications data frame and write the clustered data to a .csv file for later use.

    pubs <- merge(pubs, docTopics, by = "pmid", all.x = TRUE)
    pubs <- merge(pubs, tlist, by = "pmid", all.x = TRUE)
    write.csv(pubs, file = "clusteredData.csv", row.names = FALSE)
    
We can then get the number of documents per topic with the plyr count() function, as before.

    pubsPerTopic <- count(pubs$topic)
    pubsPerTopic

Obviously there are lots of other things that we could do with this data set, but this gives you an idea of what's possible and some of the basic code to get you there. 
