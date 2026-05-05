# AD_Entropy_MegaAnalysis
Transcriptional Network Entropy (TNE) as a Macroscopic Order Parameter in Alzheimer's Disease
Project Overview
This repository provides the formal computational implementation for the study: "Transcriptional network entropy as an order parameter for the pathological brain in Alzheimer’s disease".

The framework transitions from a reductionist gene-centric perspective toward a systems-level characterization of neurodegeneration. By applying principles of Statistical Mechanics and Information Theory, we utilize the degree-weighted Shannon entropy of the Protein-Protein Interaction (PPI) network to quantify the loss of regulatory coherence in the AD transcriptome. This approach identifies a large-magnitude elevation of network disorder, consistent with a first-order-like phase transition in the cortical regulatory landscape.

Theoretical Framework
The core of this methodology is the TNE-OP (Transcriptional Network Entropy - Order Parameter) framework. Unlike standard differential expression analysis, TNE-OP treats the transcriptome as a physical system where:

Network Topology: Information flow is constrained by a PPI backbone (STRINGdb v12.0).

Entropy Engine: Local and global disorder are calculated using a probability distribution of transcriptional activity weighted by node degree connectivity.

Regulatory Decoherence: The model specifically identifies the non-coding RNA (lncRNA) compartment as a primary driver of system-level "thawing" or loss of transcriptional order.

Repository Architecture
The pipeline is structured to ensure absolute reproducibility of the findings reported in the manuscript:

scripts/01_preprocessing.R: Multi-cohort integration, VST normalization, and empirical Bayes batch correction (limma/sva).

scripts/02_entropy_calculation.R: Mathematical core for calculating network entropy parameters and permutation testing.

scripts/03_wgcna_multi_biotype.R: Weighted gene co-expression network construction including coding and non-coding species.

scripts/04_differential_expression.R: Statistical quantification of transcriptional changes via DESeq2.

scripts/05_visualization.R: High-resolution generation of figures for publication (Vector format).

metadata/: Curated sample descriptors and clinical metadata for all discovery and validation cohorts.

Reproducibility and Data Access
The scripts are designed to be autonomous. They include routines to automatically retrieve raw count matrices from the Gene Expression Omnibus (GEO) for the following accessions:

Discovery Cohorts: GSE125050, GSE125583, GSE132177, GSE153071.

Validation Cohort: GSE132903.

System Requirements
Environment: R version 4.2.0 or higher.

Key Libraries: DESeq2, WGCNA, limma, STRINGdb, sva, ggplot2.

Hardware Note: Due to the computational intensity of the PPI network weighting and permutation tests (1,000+ iterations), a minimum of 16GB RAM is recommended.

Intellectual Property and Citation
This framework represents an original contribution to the field of systems neurobiology and physical biology.

Author: Juan M. Córdoba.
Affiliation: Department of Bioinformatics and Genomics, Universidad del Valle.
License: MIT License.

If you utilize this code, the TNE-OP framework, or the underlying theoretical concepts, you are required to cite the primary research article:

Córdoba, J. M. (2026). Transcriptional network entropy as an order parameter for the pathological brain in Alzheimer’s disease. Nature Communications (Submitted).

Contact
For technical inquiries regarding the entropy engine or requests for collaboration on the expansion of the TNE-OP framework, please contact:
juan.manuel.cordoba@correounivalle.edu.co
