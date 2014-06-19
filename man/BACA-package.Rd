\name{BACA-package}
\alias{BACA-package}
\alias{BACA}
\docType{package}
\title{
A graphical tool to visualize and compare functional annotations retrieved from DAVID knowledgebase.
}
\description{
DAVID Bioinformatics Resources (DAVID) at \url{http://david.abcc.ncifcrf.gov} is the most popular tool in functional annotation and enrichment analysis. It provides an integrated biological knowledgebase and tools to systematically extract relevant biological terms (e.g., GO terms, KEGG pathways) associated with a given gene list. After submitting a gene list, DAVID annotation tool finds the most enriched annotations and presents them in a table format. This table contains many different types of data, such as text, numbers, bars and hyperlinks, that can be hard tor read and compare when multiple enrichment results are available. The BACA package tries to address these issues by providing a novel R-based graphical tool to concisely visualize DAVID annotations and show how they change across different experimental conditions. This R package uses some functions available in the R package RDAVIDWebService available at \url{http://www.bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html}.
}
\details{
\tabular{ll}{
Package: \tab BACA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-19\cr
License: \tab GPL (>=2) \cr
}
The package BACA provides three different functions: DAVIDsearch, BBplot and Jplot. DAVIDsearch: to call the Functional Annotation Tool of DAVID. BBplot: to build a grid where each row represents an enriched annotation found by DAVID and each column the condition/treatment where that annotation was highlighted. While, each cell reports a bubble indicating the number of genes enriching the correpsonding annotation and the state of these genes in terms of down- and up-regulation (default setting: green = "down" - red = "up"). Jplot: to make a table/matrix with colored boxes. The color of each box indicated the Jaccard index based similarity score calculated between the genes sets associated with two different functional annotations.  
}
\author{
Vittorio Fortino <vittorio.fortino@ttl.fi>
}
\examples{
  data(result.kegg)
 
  bbplot.kegg <- BBplot(result.kegg, max.pval = 0.05, min.ngenes = 10, 
                    name.com = c("Cond.1_12h","Cond.1_24h","Cond.2_12h","Cond.2_24h"), 
                    labels = c("down", "up"), colors = c("#009E73", "red"), 
                    title = "BBplot - KEGG", print.term = "description")
  bbplot.kegg
  ##ggsave("KEGG_terms.tiff", width=6, height=4, scale=2, dpi=200)
  
  jplot.kegg <- Jplot(result.kegg[[4]], result.kegg[[2]], max.pval = 0.05, min.ngenes = 10) 

  jplot.kegg
  ##ggsave("j_plot.tiff", width=6, height=4, scale=2, dpi=200)
}
