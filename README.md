<p align="center">
<img src="images/FluentVCF.png" width="1000" />
</p>


**FluentVCF** is Bioconductor compliant, R package to
determine the mutation type for a set of single nucleotide variants in a given
genome.

# Installation
Run the following lines in R to install FluentVCF:

```
if (!require("devtools")) {install.packages("devtools")}
devtools::install_github("mircomacchi/FluentVCF")
```

Then load the package:

```
library(FluentVCF)
```

# Introduction

Variant Call Format (VCF) is a specification for storing genotype data
in a tab-delimited file format. Below is a high-level diagram of a
typical bioinformatics pipeline that produces a VCF file:

<img src="images/VCFlux.png" width="679" />

# Features

FluentVCF is a package which takes:

1.  a set of mutations (single nucleotide variants) in VCF format,
2.  the corresponding reference genome (e.g., human genome hg38) and
3.  a parameter “context_length” (see below) and determines for each
    mutation the corresponding mutation type as follows:
    
    The mutation type is “UP\[REF\>ALT\]DOWN” where:
    - “REF\>ALT” is the single nucleotide variant from REF base to ALT base, e.g., “C\>T” 
    - “UP” is upstream base from the reference genome (depending on the user parameter “context_length”) 
    - “DOWN” is one or more downstream bases from the reference genome (same user parameter)

The package provides an additional function, “mut_summary”, which
summarizes the mutation types for the set of mutations into a count
table (number of mutations per mutation type).

Moreover, graphical outputs are integrated.

# Requirements

FluentVCF relies on R packages to efficiently analyze and visualize variant call format (VCF) data. These dependencies include VariantAnnotation for VCF data extraction and annotation, GenomicRanges for genomic coordinate-based operations, dplyr for streamlined data manipulation, ggplot2 for creating informative plots, and RColorBrewer for aesthetically pleasing color palettes. 

### Dependencies installation
```
# Define the list of packages to be installed
required_packages <- c("VariantAnnotation", "GenomicRanges", "dplyr", "ggplot2", "RColorBrewer")

# Check if each package is already installed, and install if not
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Load the installed packages
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# You can now use FluentVCF with the required packages loaded

```

