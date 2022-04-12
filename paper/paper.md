---
title: 'VIEWpoly: a visualization tool to integrate and explore results of polyploid
  genetic analysis'
tags:
- R
- Shiny
- Linkage map
- QTL
- Polyploid
date: "18 February 2021"
output: pdf_document
authors:
- name: Cristiane Hayumi Taniguti
  affiliation: 1
- name: Gabriel Gesteira
  affiliation: 2
- name: Jeekin Lau
  affiliation: 1
- name: Guilherme da Silva Pereira
  affiliation: 3
- name: Zhao-Bang Zeng
  affiliation: 2
- name: David Byrne
  affiliation: 1
- name: Oscar Riera-Lizarazu
  affiliation: 1
- name: Marcelo Mollinari
  affiliation: 2
bibliography: paper.bib
affiliations:
- name: Department of Horticultural Sciences, Texas A&M University, College Station,
    TX, USA
  index: 1
- name: Bioinformatics Research Center, Department of Horticultural Sciences, North
    Carolina State University, Raleigh, NC, USA
  index: 2
- name: Department of Agronomy, Federal University of Viçosa, Brazil
  index: 3
---

# Summary

Advancements in computational tools for linkage and quantitative trait loci (QTL) analysis in autopolyploids have allowed the in-depth exploration and understanding of their genetics. Despite providing advanced methods, these tools may be challenging to use when interpreting and connecting results to other sources of genomic information and when attempting to efficiently apply them in practical applications. `VIEWpoly` is an R package for visualizing, exploring, and integrating results from different polyploid computational tools using an interactive graphical user interface. The package allows users to explore linkage and QTL analysis results and integrates these with genomic information in a genome browser, facilitating the search for candidate genes. In addition, it provides features to download comprehensive tables and graphics for further analysis and presentation. `VIEWpoly` is freely available as an R package at https://CRAN.R-project.org/package=viewpoly. 

# Statement of need

In recent years, genetic analysis in autopolyploid crops has been gaining attention due to the flourishing of high-throughput genotyping technologies, which deliver massive amounts of DNA sequences. This technological advance prompted the development of many computational tools to perform an in-depth analysis of these complex genomes. Programs such as `TetraploidSNPmap` [@Hackett2017], `polymapR` [@Bourke2018], and `MAPpoly` [@Mollinari2019; @Mollinari2020] use dosage-based markers to build genetic maps and obtain the linkage phases in full-sib mapping populations. `TetraOrigin` [@Zheng2016], `PolyOrigin` [@Zheng2021], `MAPpoly` [@Mollinari2019], and `polyqtlR` [@Bourke2021] provide the probabilistic haplotype inheritance profiles of individuals in mapping populations in terms of their respective parental haplotypes. These probabilities are used to estimate the position and genetic effects of quantitative trait loci (QTL). Using this information, researchers can benefit from the inferred genome-wide genotype and phenotype relationship to make informed breeding decisions. Software programs for performing such analysis in autopolyploids include `GWASpoly` [@Rosyara2016], `QTLpoly` [@Pereira2020], `polyqtlR` [@Bourke2021], and `diaQTL` [@Amadeu2021]. Despite the increasing availability of computational and analytical tools focused on autopolyploid species, the need for programming skills to comprehensively organize, integrate, display results of these programs still imposes an obstacle to exploring and using this information in practical situations. We present `VIEWpoly`, a user-friendly R package that allows easy integration, visualization, and exploration of upstream genetic and genomic analysis focused on polyploid organisms. `VIEWpoly` is a self-contained interactive R package built as `Shiny` modules using `golem` [@Fay2022] framework. It includes documentation, data examples, and comprehensive video tutorials to guide users in easy-to-follow steps to extract information from polyploid computational tools in a seamless manner. It can be deployed locally or on servers.

![`VIEWpoly` conceptual flowchart highlighting the inputs and the three modules available: `VIEWqtl`, `VIEWgenome`, and `VIEWmap`.\label{fig:flowchart}](viewpoly-flowchart.png)

# Software implementation

`VIEWpoly` comprises three main modules arranged in an intuitive interface that guides the user through the analysis \autoref{fig:flowchart}. In its current version (V 0.1.1), the package supports `polymapR` [@Bourke2018] and `MAPpoly` [@Mollinari2019; @Mollinari2020] results to display genetic maps. The QTL analysis results can be imported from `QTLpoly` [@Pereira2020], `polyqtlR` [@Bourke2021], and `diaQTL` [@Amadeu2021] packages. Maps and QTL input files can be provided in the RData format or standard text formats (CSV or TSV). `VIEWpoly` is also compatible with genome-related files supported by `JBrowseR`, including FASTA, BED, GFF3, BAM, CRAM, VCF, Wiggle, and bigWig files [@Hershberg2021]. Once the data is uploaded, the user can export the information in a single compressed file to expedite the upload process in further analysis.

## `VIEWqtl`: the QTL Browser

The `VIEWqtl` module is designed to explore information from the QTL analysis. The first graphic contains the QTL profile, where users can display single or multiple linkage groups and phenotypes. The QTL peak positions are represented by triangles at the bottom of the graphic and can be interactively selected to explore the genetic effects of parental haplotypes. It is possible to compare pleiotropic and linked QTL effects simultaneously and select individuals carrying alleles with specific effects. In this module, users can also find a table with a summary of selected QTL characteristics, and progeny breeding values will be displayed if results from `QTLpoly` are available.
The information shown in graphics and tables may change according to the programs used in the upstream analysis. For example, when `QTLpoly` is used, the QTL profile will show the negative logarithm of the P values (LOP = −log10 P) associated with the score statistics used to test for QTL significance. When `polyqtlR` is used, the logarithm of the odds (LOD) ratio is shown instead. If `diaQTL` is used, the significance curve can be evaluated using the Deviance Information Criterion (DIC). Other results derived from QTL analysis, such as heritabilities, allele effects, breeding values, and confidence intervals, are available and shown according to the upstream method used. All information is displayed in tables that can be downloaded in the CSV, EXCEL, and PDF formats, and interactive graphics in static versions can be exported in the PNG, TIFF, JPEG, and PDF formats.

## `VIEWgenome`: the Genome Browser

This module lets users explore the relationship between QTL regions and a genome sequence with its features such as gene annotations, transcripts, and domains. After selecting a range in the genetic map, the `JBrowseR` interface [@Hershberg2021] displays all tracks identified in the uploaded files. The display will be set to meet the corresponding genomic region according to the selected map region. The region can be queried and explored while all related information is automatically updated. Users can select the genomic sequence content and all annotated features for any given genomic range, including QTL peaks and their supporting intervals. A scatterplot of the genetic versus the physical positions highlights the relationship between the genetic map, the reference genome, the QTL positions and intervals, and the recombination rate along the chromosomes.

## `VIEWmap`: the Map Browser

In the `VIEWmap` module, users can interactively select a genomic range based on a QTL profile and explore specific regions of the map. Given a selected region, it is possible to interactively explore the estimated haplotypes and inspect their marker positions, allele dosages, genetic distances, and sequence base contents for parental homologs. Users can optionally generate a table with the haplotypes and download them in the CSV, EXCEL, and PDF formats. The module highlights all linkage groups in a graphical representation with their respective marker positions and estimated recombination fractions, allowing for a direct comparison between homology groups. The module also shows a summarized table with information about the genetic maps, including their estimated sizes, densities, and genetic and genomic composition. Both figures and tables can be downloaded by users in standard static formats, such as CSV, EXCEL, PNG, TIFF, JPEG, and PDF.

# Conclusions and perspectives

`VIEWpoly` facilitates the integration, exploitation, and use of phenotypic, genetic, and genomic information provided by upstream analytical tools, leveraging the adoption of genomic-based tools for polyploid species. Its easy-to-use interface assists the decision-making process in practical applications, such as marker-assisted selection in breeding programs. Due to its modular nature, the package can be expanded and updated to incorporate further developments in diploid and polyploid genetic software. We expect to provide continuous integration and exploitation of novel features and tools, such as the analysis of multi-parental and multi-generation populations.

#	Acknowledgments

We acknowledge Tessa Hochhaus for voicing-over the tutorial video.

# Funding

This work was supported by the USDA, National Institute of Food and Agriculture (NIFA), Specialty Crop Research Initiative (SCRI) project ‘‘Tools for Genomics-Assisted Breeding in Polyploids: Development of a Community Resource’’ (Award No. 2020-51181-32156), and by the Bill & Melinda Gates Foundation under the Genetic Advances and Innovative Seed Systems for Sweetpotato project (SweetGAINS) (grant number OPP1213329).

# References
