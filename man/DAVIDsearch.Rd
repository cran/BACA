\name{DAVIDsearch}
\alias{DAVIDsearch}
\title{
Get functional annotation charts of gene lists using DAVID.
}
\description{
Uploads multiple gene lists to DAVID, then perform an automated enrichment analysis based on a given database/resource (i.e., GO terms, KEGG, etc.) for each gene list. 
Returns DAVIDFunctionalAnnotationChart objects, one for each gene list. Requires a DAVID user id first
\url{http://david.abcc.ncifcrf.gov/webservice/register.htm}. 
Among the optional input paramters, can use the easeScore (or P-value) to do thresholding on results. However, we suggest to get DAVID to return all possible annotations 
despite non-significant P-values (easeScore = 1) and perform your own thresholding when using the plot functions.

}
\usage{
DAVIDsearch(gene.lists, david.user, idType = "AFFYMETRIX_3PRIME_IVT_ID", 
			listType = "Gene", easeScore = 1, annotation = "KEGG_PATHWAY", 
			species = NA)
}
\arguments{
  \item{gene.lists}{
   (Required). Lists of character vectors representing the gene lists. 
   To compare the enrichment of differential genes specify for each list the up- and down- regulated genes, separately.
}
  \item{david.user}{
   (Required). Character vector to specify the email of a given registered DAVID; it is necessary to use the DAVID web service.
}
  \item{idType}{
   (Optional). Character with the type of submitted ids. 
   Default value is "AFFYMETRIX_3PRIME_IVT_ID".
}
  \item{listType}{
   (Optional). Character with the type of list (Gene, Background). Default value is "Gene".
}   
  \item{easeScore}{
   (Optional).  Numeric with the EASE score. Default value is 1 (see DAVID Help page).
}
  \item{annotation}{
   (Optional). Character vector with the category name to use in the functional annotation analysis. Default value is "KEGG_PATHWAY".
}
  \item{species}{
   (Optional). Character vector with the species to use. Default value is NA. This means that DAVID will use the specie found for the submitted gene list.
}
}
\details{
 The following information are printed out during the the querying process: 
 For the list "name.list" you have:
   - Number of genes loaded = #.
   - Number of genes mapped/annotated in DAVID <inDavid> = #.
   - Number of unmapped genes <unmappedIds> = #.
   - Species involved.
}
\value{
 List of DAVIDFunctionalAnnotationChart objects, one for each specified list of genes.
}
\author{
Vittorio Fortino
}
\seealso{
\code{\link{BBplot}}
\code{\link{Jplot}}
}
\examples{
  #data(gene.lists.ex)
  #result.kegg <- DAVIDsearch(gene.lists.ex, david.user = "vittorio.fortino@ttl.fi") 
}
