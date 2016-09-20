require(knitr)
require(tools)
knit('article.Rnw')
texi2pdf('article.tex')
