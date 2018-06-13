all : index.html

index.html : index.Rmd data_cleaning.knit.md QAanalysis.knit.md DE_SVA.knit.md FE_SVA.knit.md DEanalysis_limma-voom_noSVA.knit.md FunctionaEnrichment_noSVA.knit.md
	Rscript -e "rmarkdown::render('$<')"

data_cleaning.knit.md : data_cleaning.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

QAanalysis.knit.md : QAanalysis.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

DE_SVA.knit.md : DE_SVA.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

FE_SVA.knit.md : FE_SVA.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

DEanalysis_limma-voom_noSVA.knit.md : DEanalysis_limma-voom_noSVA.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

FunctionaEnrichment_noSVA.knit.md : FunctionaEnrichment_noSVA.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"
