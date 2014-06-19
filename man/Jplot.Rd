\name{Jplot}
\alias{Jplot}
\title{
To compute and visualize Jaccard similarity score between two enrichment analysis results.  
}
\description{
Comparing two enrichment analysis results, it reveals the relationships existings between annotations found by using two different gene lists. The Jplot() function takes in input two different DAVIDFunctionalAnnotationChart objects (1 and 2) and provides in output a table/matrix with colored boxes. Each box reports the Jaccard index-based similarity score computed between the gene sets enriching two functional annotations. 
}
\usage{
Jplot(david.obj.1, david.obj.2,  max.pval = 0.01, min.ngenes = 5, title = "Jplot", 
       print.term ="full")
}
\arguments{
  \item{david.obj.1}{
   (Required). First DAVIDFunctionalAnnotationChart objects.
}
  \item{david.obj.2}{
   (Required). Second DAVIDFunctionalAnnotationChart objects.
}
  \item{max.pval}{
   (Optional). Numeric with the p-value that must be present in order to select the most significant enriched annotations. Default values is 0.01.
}
  \item{min.ngenes}{
   (Optional). Integer with the minimum number of genes (greater or equal) that must be present in an enriched annotations. Default values is 5.
}
  \item{title}{
   (Optional). Character with the name of the plot. Default valuse is "Bubble plot". 
}
  \item{print.term}{
   (Optional). Character vector to specify the term used to indicate the annotation: 'full', 'name' or 'description'. The default value is "full". 
}
}
\value{
Return a ggplot2 graph.
}
\author{
Vittorio Fortino
}
\note{
The two DAVIDFunctionalAnnotationChart objects should be associated to two different enrichment analysis results. 
}

\seealso{
\code{\link{DAVIDsearch}} 
\code{\link{BBplot}} 
} 
\examples{
  data(result.kegg)
  jplot.kegg <- Jplot(result.kegg[[4]], result.kegg[[2]], max.pval = 0.05, min.ngenes = 10) 
  jplot.kegg
}
