<!-- badges: start -->
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/viewpoly)](https://cran.r-project.org/package=viewpoly)
[![R-universe PolyVerse Status Badge](https://polyploids.r-universe.dev/badges/viewpoly)](https://polyploids.r-universe.dev/badges/viewpoly)
[![codecov](https://codecov.io/github/mmollina/viewpoly/branch/main/graphs/badge.svg)](https://codecov.io/github/mmollina/viewpoly)
[![CRAN_monthly_downloads](https://cranlogs.r-pkg.org/badges/viewpoly)](https://cranlogs.r-pkg.org/badges/viewpoly)
<!-- badges: end -->
  
# VIEWpoly <img src="https://user-images.githubusercontent.com/7572527/145726577-7b01d48b-ca1d-446b-b9c8-aff8c3c9877b.png" align="right" width="230"/>

`VIEWpoly` is a shiny app and R package for visualizing and exploring results from [polyploid computational tools](https://www.polyploids.org/) using an interactive graphical user interface. The package allows users to directly upload output files from [polymapR](https://CRAN.R-project.org/package=polymapR), [MAPpoly](https://CRAN.R-project.org/package=mappoly) , [polyqtlR](https://CRAN.R-project.org/package=polyqtlR), [QTLpoly](https://CRAN.R-project.org/package=qtlpoly), 
[diaQTL](https://github.com/jendelman/diaQTL) and genomic assembly, variants, annotation and alignment files. VIEWpoly uses [shiny](https://CRAN.R-project.org/package=shiny), [golem](https://CRAN.R-project.org/package=golem), [ggplot2](https://CRAN.R-project.org/package=ggplot2), [plotly](https://CRAN.R-project.org/package=plotly), and [JBrowseR]( https://CRAN.R-project.org/package=JBrowseR) libraries to graphically display the QTL profiles, positions, alleles estimated effects, progeny individuals containing specific haplotypes and their breeding values. It is also possible to access marker dosage and parental phase from the linkage map. If genomic information is available, the corresponding QTL positions are interactively explored using JBrowseR interface, allowing the search for candidate genes. It also provides features to download specific information into comprehensive tables and images for further analysis and presentation.

### Quick Start

The quickest way of accessing `VIEWpoly` is [here](https://cris-taniguti.shinyapps.io/viewpoly/). However, our shinyapps.io does not upload files larger than 1GB. If you have larger datasets, you will need to run `VIEWpoly` locally.

You can run `VIEWpoly` locally installing the package and accessing the graphical interface through a web browser. To use the stable version, please install the package from CRAN:

```{r}
install.packages("viewpoly")
viewpoly::run_app()
```

If you want to use the latest development version, go ahead and install `VIEWpoly` from our Github repository:

```{r}
# install.packages("devtools")
devtools::install_github("mmollina/viewpoly")
viewpoly::run_app()
```

NOTE: Windows users may need to install the `Rtools` before compiling the package from source (development version).

The `Input data` tab has options for several types of inputs. You can upload directly outputs from:

* [MAPpoly](https://CRAN.R-project.org/package=mappoly)
* [polymapR](https://CRAN.R-project.org/package=polymapR)
* [polyqtlR](https://CRAN.R-project.org/package=polyqtlR)
* [QTLpoly](https://CRAN.R-project.org/package=qtlpoly)
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

Access the [tutorial](https://cristianetaniguti.github.io/viewpoly_vignettes/VIEWpoly_tutorial.html). 

We also presented the app main features in [this video](https://www.youtube.com/watch?v=OBt_jebhfeY)

If you would like to contribute to develop `VIEWpoly`, please check our [Contributing Guidelines](https://cristianetaniguti.github.io/viewpoly_vignettes/Contributing_guidelines.html).

### References

Taniguti CH, Gesteira GS, Lau J, Pereira GS, Zeng ZB, Byrne D, Riera-Lizarazu O, Mollinari M. VIEWpoly: a visualization tool to integrate and explore results of polyploid genetic analysis. Submitted.

Mollinari M, Garcia AAF. 2019. “Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models.” G3: Genes, Genomes, Genetics 9 (10): 3297-3314. doi:10.1534/g3.119.400378.

Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB. 2020. “Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population.” Genetics 215 (3): 579-595. doi:10.1534/genetics.120.303080.

Amadeu RR, Muñoz PR , Zheng C, Endelman JB. 2021."QTL mapping in outbred tetraploid (and diploid) diallel populations." Genetics 219 (3), iyab124, https://doi.org/10.1093/genetics/iyab124

Bourke PM , van Geest G, Voorrips RE, Jansen J, Kranenburg T, Shahin A, Visser RGF , Arens P, Smulders MJM , Maliepaard C. 2018."polymapR—linkage analysis and genetic map construction from F1 populations of outcrossing polyploids." Bioinformatics, 34 (20): 3496–3502, https://doi.org/10.1093/bioinformatics/bty371

Bourke PM, Voorrips RE, Hackett CA, van Geest G, Willemsen JH, Arens P, Smulders MJM, Visser RGF, Maliepaard C. 2021."Detecting quantitative trait loci and exploring chromosomal pairing in autopolyploids using polyqtlR." Bioinformatics, 37 (21): 3822–3829, https://doi.org/10.1093/bioinformatics/btab574

### Acknowledgment

VIEWpoly project is supported by the USDA, National Institute of Food and Agriculture (NIFA), Specialty Crop Research Initiative (SCRI) project [‘‘Tools for Genomics-Assisted Breeding in Polyploids: Development of a Community Resource’’](https://www.polyploids.org/)  and by the Bill & Melinda Gates Foundation under the Genetic Advances and [Innovative Seed Systems for Sweetpotato project (SweetGAINS)](https://cgspace.cgiar.org/handle/10568/106838).
