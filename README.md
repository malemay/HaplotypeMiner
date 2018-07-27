# HaplotypeMiner

## What is HaplotypeMiner

*HaplotypeMiner* is an R package developed for exploring allelic diversity at genes of interest in a plant breeding context. The program minimally takes as input a dataset of SNP markers generated through various methods (e.g. genotyping-by-sequencing [GBS] or SNP arrays) and the genomic position of a gene of interest, and outputs a set of possible haplotypes defined by the genotypes of a reduced number of neighboring SNPs. The kinship and structure of the population assessed can also by used as input to *HaplotypeMiner* to yield more robust results.

The haplotyping model of *HaplotypeMiner* implies the following steps:
- SNP markers in a window of a given size surrounding the central position of the gene are extracted and optionally filtered according to user-specified settings.
- Linkage disequilibrium (LD) blocks are identified on each side (5' and 3') of the gene center. One SNP is selected from each block in order to avoid redundancy and reduce the dataset to a set of informative SNPs.
- Pairs of SNPs that are in significant LD across the gene center are selected and used for defining haplotypes.
- Every unique combination of alleles at the SNPs selected during the preceding step is considered as a haplotype and is output by the program.

The underlying assumption of *HaplotypeMiner* is that even though GBS or SNP array datasets may not provide a comprehensive view of the variation in gene-coding or other functionally relevant sequences, they may provide polymorphic markers that are in LD with variants having a functional impact. The model implemented by *HaplotypeMiner* assumes that if two markers are in LD *across* the central position of the gene, then they are even more likely to be in LD with variants located in the gene sequence. By identifying unique combinations of such marker pairs, our hope is therefore to identify haplotypes that have a one-to-one correspondence with functionally relevant alleles of genes of interest. *HaplotypeMiner* is therefore primarily intended as a tool to allow plant breeders to assess allelic diversity at specific genes in a germplasm collection and assist them in making decisions.

## Installation

*HaplotypeMiner* can be installed directly from within R by calling `devtools::install_github("malemay/HaplotypeMiner")`. This will directly fetch the package from the GitHub and install it on your computer. This requires the package [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) to be installed on your computer. This package is available from CRAN through the usual `ìnstall.packages()` interface.

One can also clone this repository from GitHub using `git clone https://github.com/malemay/HaplotypeMiner.git` and then install the source package from within R using `install.packages("HaplotypeMiner", repos = NULL, type = "source")`.

*HaplotypeMiner* requires a few [Bioconductor](https://www.bioconductor.org/) packages to be installed on your computer. You can install them by running the following commands in R:

```r
{
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GenomeInfoDb", "snpStats", "SummarizedExperiment", "VariantAnnotation"))
}
```

An unresolved issue requires package `snpStats` to be explicitly loaded with `library(snpStats)` in addition to `library(HaplotypeMiner)` whenever *HaplotypeMiner* is to be used. This does not, however, impact the normal functioning of the package.

## Typical usage

Using *HaplotypeMiner* will typically involve three steps for most users of the package:

* Creation of an object that stores analysis parameters using the helper function `haplo_params`. These parameters include the location of input files, the variant filtering parameters, as well as the various parameters used in the definition of haplotypes. Arguments to `haplo_params` can be obtained by running `args(haplo_params)` and the documentation of the different parameters can be consulted by running `?haplo_params`.
* Launching the analysis with the function `haplo_selection`, using the parameters object generated at the previous step as input.
* Automated generation of text files and figures describing the results of the analysis by applying the function `haplo_output` to the output of the call to `haplo_selection`.

Users interested in gaining more control over the workflow and the output of the analysis can also break the tasks performed by `haplo_selection` and `haplo_output` into different parts, as the functions called by these can also be called on their own. A vignette describing the use of these functions in more detail will be added to the package shortly. Similarly, detailed documentation has yet to be added for a few functions; this will be completed as soon as possible. 

## Notes

If you use *HaplotypeMiner* as part of your work, please cite it as numerous hours have been invested in its development. A manuscript describing the underlying model of *HaplotypeMiner* and results obtained with soybean will be submitted for publication as a research paper shortly. Until then, the recommended citation for the package can be obtained by running the command `citation("HaplotypeMiner")` from within R.

This software is provided without any guarantee. *HaplotypeMiner* has only been thoroughly tested with soybean. We do not know to what extent the model will hold for polyploid species, mainly outcrossing species, or species in which linkage disequilibrium decays more rapidly with physical distance than in soybean. If you test *HaplotypeMiner* with a different species, we would be happy to know about the results and provide some advice if needed.

Issues, bugs reports and questions can be shared on the GitHub page of the project or addressed to the package maintainer, Marc-André Lemay (e-mail address in the package DESCRIPTION).

