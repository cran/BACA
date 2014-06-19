## ----message=FALSE, results='hide'---------------------------------------
library(BACA)

## ------------------------------------------------------------------------
data(gene.lists.ex)

## ------------------------------------------------------------------------
str(gene.lists.ex)

## ------------------------------------------------------------------------
result.kegg <- DAVIDsearch(gene.lists.ex, david.user = "vittorio.fortino@ttl.fi", idType="ENTREZ_GENE_ID", annotation="KEGG_PATHWAY")

## ------------------------------------------------------------------------
david.obj <- DAVIDWebService$new(email="vittorio.fortino@ttl.fi")
getAllAnnotationCategoryNames(david.obj)

## ------------------------------------------------------------------------
getIdTypes(david.obj)

## ------------------------------------------------------------------------
bbplot.kegg <- BBplot(result.kegg, max.pval = 0.05, min.ngenes = 10, 
                    name.com = c("Cond.1_12h","Cond.1_24h","Cond.2_12h","Cond.2_24h"), 
                    labels = c("down", "up"), colors = c("#009E73", "red"), 
                    title = "BBplot - KEGG")

## ----, fig.width=8, fig.height=8-----------------------------------------
bbplot.kegg

## ------------------------------------------------------------------------
ggsave("KEGG_terms.tiff", width=6, height=4, scale=2, dpi=200)

## ------------------------------------------------------------------------
jplot.kegg <- Jplot(result.kegg[[4]], result.kegg[[2]], max.pval = 0.05, min.ngenes = 10)

## ----, fig.width=12, fig.height=6----------------------------------------
jplot.kegg

