DNA Methylation and Expression Network Analysis with SMITE
================================

##### Author: Samantha Schaffner
##### Date: Sept 16, 2021

## Overview of SMITE

The Significance-based Modules Integrating the Transcriptome and Epigenome (SMITE) package ( [Wijetunga et al. 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1477-3) ) can be used to integrate multiple types of 'omics data (primarily focused on DNA methylation and expression) with gene scoring and network analysis. The objectives of SMITE are to determine which genes have changes to DNA methylation and/or expression between experimental groups, and then examine with network analysis how the encoded proteins interact with one another. The full SMITE workflow will return protein-protein interaction "modules" containing genes which have altered DNA methylation and/or expression, i.e. genes that are hypothesized to be co-regulated by multiple 'omics. 

SMITE builds upon a previous method called Functional Epigenetic Modules (FEM; [Jiao et al. 2014](https://academic.oup.com/bioinformatics/article/30/16/2360/2748243) ) by introducing a much wider array of customization in terms of which gene features you are interested in, how important each feature is, and whether or not you want to specify a directional relationship between DNA methylation and expression.


## Main steps in workflow

**1. Gather DNA methylation and expression information**

The first step is to create a "PvalueAnnotation" object which can store the p-values from different 'omics analyses and genome annotation information.

**2. Combine site-specific DNA methylation p-values to obtain gene-level p-values**

The second step performs data reduction on the CpG-level DNA methylation data, with the goal of obtaining a unit that is meaningful to combine with expression. SMITE will combine DNA methylation p-values across gene features (default promoter and gene body).

**3. Score genes**

At this stage, DNA methylation and expression p-values are combined to create a single gene-level score for each gene. The DNA methylation data can be weighted by gene feature, and the directional relationship between DNA methylation at each feature and expression can be optionally specificied.

**4. Find regulatory modules**

The next stage annotates gene scores onto pre-existing protein-protein interaction networks, which can be retrived from a database like REACTOME or STRING. Each gene represents a node in the network, and connections between genes are edges. A spin-glass algorithm is applied to find sub-networks of highly scoring genes.

**5. Module interpretation**

Modules can be plotted, and optionally you can perform pathway enrichment on the modules.


## Data for tutorial

We will be using DNA methylation and gene expression data from Lund human mesencephalic (LUHMES) cells, a fetal midbrain cell line that is commonly differentiated to dopaminergic neurons to model Parkinson's disease in vitro. We'll be looking at the effect of overexpressing alpha-synuclein, a gene implicated in Parkinson's disease. DNA methylation data was generated with the Illumina EPIC BeadChip microarray (Kobor lab), and gene expression data was generated with RNA-seq (Dr. Tiago Outeiro, University of Goettingen). For more information, see [Paiva et al., 2017](https://academic.oup.com/hmg/article/26/12/2231/3084502) and [Schaffner et al., 2021](https://www.biorxiv.org/content/10.1101/2021.06.12.448150v1) (preprint).

Although additional 'omics (DNA hydroxymethylation and ChIP-seq) exist for these cells, we'll focus on two for the sake of the tutorial. See the [code for Schaffner et al., 2021](https://github.com/samschaf/LUHMES) for a more detailed example of integrating several 'omics with SMITE.
