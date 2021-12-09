<!-- badges: start -->
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![R-CMD-check](https://github.com/mmollina/viewpoly/workflows/R-CMD-check/badge.svg)](https://github.com/mmollina/viewpoly/actions)
[![codecov](https://codecov.io/gh/mmollina/viewpoly/branch/main/graph/badge.svg)](https://codecov.io/gh/mmollina/viewpoly)
<!-- badges: end -->
  
# VIEWpoly

`VIEWpoly` is a shiny app and R package for visualizing and exploring results from [polyploid computational tools](https://www.polyploids.org/) using an interactive graphical user interface. The package allows users to directly upload output files from [polymapR](https://cran.r-project.org/web/packages/polymapR/index.html), [MAPpoly](https://cran.r-project.org/web/packages/mappoly/index.html) , [polyqtlR](https://cran.r-project.org/web/packages/polyqtlR/index.html) , [QTLpoly](https://cran.r-project.org/web/packages/qtlpoly/index.html), 
[diaQTL](https://github.com/jendelman/diaQTL) and genomic assembly, variants, annotation and alignment files. VIEWpoly uses [shiny](https://cran.r-project.org/web/packages/shiny/index.html), [golem](https://cran.r-project.org/web/packages/golem/index.html), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [plotly](https://cran.r-project.org/web/packages/plotly/index.html), and [JBrowseR](https://cran.r-project.org/web/packages/JBrowseR/index.html) libraries to graphically display the QTL profiles, positions, alleles estimated effects, progeny individuals containing specific haplotypes and their breeding values. It is also possible to access marker dosage and parental phase from the linkage map. If genomic information is available, the corresponding QTL positions are interactively explored using JBrowseR interface, allowing the search for candidate genes. It also provides features to download specific information into comprehensive tables and images for further analysis and presentation.

### Quick Start

You can run `VIEWpoly` locally installing the package and accessing the graphical interface through a web browser:

```{r}
# install.packages("devtools")
devtools::install_github("mmollina/viewpoly")
viewpoly::run_app()
```

The `Input data` tab has options for diverse types of inputs. You can upload directly outputs from:

* [MAPpoly](https://cran.r-project.org/web/packages/mappoly/index.html)
* [polymapR](https://cran.r-project.org/web/packages/polymapR/index.html)
* [polyqtlR](https://cran.r-project.org/web/packages/polyqtlR/index.html)
* [QTLpoly](https://cran.r-project.org/web/packages/qtlpoly/index.html)
* [diaQTL](https://github.com/jendelman/diaQTL)
* CSV, TSV or TSV.GZ standard formats

To relate the genetic maps and QTL analysis with genomic information, it is also required:

* FASTA reference genome

It is optional to upload also: 

* GFF3 annotation file
* BAM or CRAM alignment file
* VCF file
* bigWig file

### Documentation

Access the [tutorial](). 

We present the app main features in the video bellow:

<iframe width="600" height="315"
src="https://www.youtube.com/embed/yqWX86uT5jM">
</iframe> 

