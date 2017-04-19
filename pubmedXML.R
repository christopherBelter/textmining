## clean pubmed XML returned from either the reutils or rentrez packages and save the cleaned XML to a new file
clean_api_xml <- function(infile, outfile) {
	theData <- readChar(infile, file.info(infile)$size, useBytes = TRUE)
	theData <- gsub("<?xml version=\"1.0\" ?>", "", theData, fixed = TRUE)
	theData <- gsub("<!DOCTYPE PubmedArticleSet PUBLIC \"-//NLM//DTD PubMedArticle, 1st January 2017//EN\" \"https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_170101.dtd\">", "", theData, fixed = TRUE, useBytes = TRUE)
	theData <- gsub("<PubmedArticleSet>", "", theData, fixed = TRUE)
	theData <- gsub("</PubmedArticleSet>", "", theData, fixed = TRUE)
	theData <- gsub("<U\\+\\w{4}>", "", theData) ## note: with some files this doesn't catch everything; potial issue with <OtherAbstract> tags especially
	theData <- paste("<?xml version=\"1.0\" ?>", "<!DOCTYPE PubmedArticleSet PUBLIC \"-//NLM//DTD PubMedArticle, 1st January 2017//EN\" \"https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_170101.dtd\">", "<PubmedArticleSet>", theData, sep = "\n")
	theData <- paste(theData, "</PubmedArticleSet>")
	theData <- iconv(theData, to = "UTF-8", sub = "")
	writeLines(theData, outfile, sep = " ")
	return(theData)
}

## extract a data frame from the cleaned XML
## Note: does not handle <pubmedBookArticle> documents
extract_xml <- function(theFile) {
	library(XML)
	newData <- xmlParse(theFile)
	records <- getNodeSet(newData, "//PubmedArticle")
	pmid <- xpathSApply(newData,"//MedlineCitation/PMID", xmlValue)
	doi <- lapply(records, xpathSApply, ".//ELocationID[@EIdType = \"doi\"]", xmlValue)
	doi[sapply(doi, is.list)] <- NA
	doi <- unlist(doi)
	authLast <- lapply(records, xpathSApply, ".//Author/LastName", xmlValue)
	authLast[sapply(authLast, is.list)] <- NA
	authInit <- lapply(records, xpathSApply, ".//Author/Initials", xmlValue)
	authInit[sapply(authInit, is.list)] <- NA
	authors <- mapply(paste, authLast, authInit, collapse = "|")
	## affiliations <- lapply(records, xpathSApply, ".//Author/AffiliationInfo/Affiliation", xmlValue)
	## affiliations[sapply(affiliations, is.list)] <- NA
	## affiliations <- sapply(affiliations, paste, collapse = "|")
	year <- lapply(records, xpathSApply, ".//PubDate/Year", xmlValue) 
	year[sapply(year, is.list)] <- NA
	year <- unlist(year)
	articletitle <- lapply(records, xpathSApply, ".//ArticleTitle", xmlValue) 
	articletitle[sapply(articletitle, is.list)] <- NA
	articletitle <- unlist(articletitle)
	journal <- lapply(records, xpathSApply, ".//ISOAbbreviation", xmlValue) 
	journal[sapply(journal, is.list)] <- NA
	journal <- unlist(journal)
	volume <- lapply(records, xpathSApply, ".//JournalIssue/Volume", xmlValue)
	volume[sapply(volume, is.list)] <- NA
	volume <- unlist(volume)
	issue <- lapply(records, xpathSApply, ".//JournalIssue/Issue", xmlValue)
	issue[sapply(issue, is.list)] <- NA
	issue <- unlist(issue)
	pages <- lapply(records, xpathSApply, ".//MedlinePgn", xmlValue)
	pages[sapply(pages, is.list)] <- NA
	pages <- unlist(pages)
	abstract <- lapply(records, xpathSApply, ".//Abstract/AbstractText", xmlValue)
	abstract[sapply(abstract, is.list)] <- NA
	abstract <- sapply(abstract, paste, collapse = "|")
	meshHeadings <- lapply(records, xpathSApply, ".//DescriptorName", xmlValue)
	meshHeadings[sapply(meshHeadings, is.list)] <- NA
	meshHeadings <- sapply(meshHeadings, paste, collapse = "|")
	grantAgency <- lapply(records, xpathSApply, ".//Grant/Agency", xmlValue)
	grantAgency[sapply(grantAgency, is.list)] <- NA
	grantAgency <- sapply(grantAgency, paste, collapse = "|")
	grantAgency <- sapply(strsplit(grantAgency, "|", fixed = TRUE), unique)
	grantAgency <- sapply(grantAgency, paste, collapse = "|")
	grantNumber <- lapply(records, xpathSApply, ".//Grant/GrantID", xmlValue)
	grantNumber[sapply(grantNumber, is.list)] <- NA
	grantNumber <- sapply(grantNumber, paste, collapse = "|")
	grantCountry <- lapply(records, xpathSApply, ".//Grant/Country", xmlValue)
	grantCountry[sapply(grantCountry, is.list)] <- NA
	grantCountry <- sapply(grantCountry, paste, collapse = "|")
	grantCountry <- sapply(strsplit(grantCountry, "|", fixed = TRUE), unique)
	grantCountry <- sapply(grantCountry, paste, collapse = "|")
	ptype <- lapply(records, xpathSApply, ".//PublicationType", xmlValue)
	ptype[sapply(ptype, is.list)] <- NA
	ptype <- sapply(ptype, paste, collapse = "|")
	theDF <- data.frame(pmid, doi, authors, year, articletitle, journal, volume, issue, pages, abstract, meshHeadings, grantAgency, grantNumber, grantCountry, ptype, stringsAsFactors = FALSE)
	return(theDF)
}