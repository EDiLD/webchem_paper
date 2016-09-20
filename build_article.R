require(knitr)
require(tools)
knit('revision1/article.Rnw')
texi2pdf('revision1/article.tex')
