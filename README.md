# Analysis of *Mapping of cis-regulatory variants by differential allelic expression analysis identifies candidate causal variants and target genes of 41 breast cancer risk loci*

This repository contains the code and data used in the analysis presented in "Mapping of cis-regulatory variants by differential allelic expression analysis identifies candidate causal variants and target genes of 41 breast cancer risk loci (https://www.medrxiv.org/content/10.1101/2022.03.08.22271889v2)", published on medRxiv. The study aims to unravel genetic risk factors for breast cancer by employing a unique approach—integrating genome-wide differential allelic expression (DAE) analysis with genome-wide association studies (GWAS). Through this, the research identifies variants with regulatory potential in normal breast tissue. It maps these variants (daeQTLs) and intersects them with GWAS data to pinpoint candidate risk regulatory variants (risk-daeQTLs) in active regulatory regions. The study also reveals new candidate target genes within these loci, shedding light on the intricate genetic regulatory landscapes associated with breast cancer. As a validation, the researchers functionally characterize specific causal variants, demonstrating the potential of DAE analysis for understanding breast cancer genetic risk.

## Citation

Please cite the following paper when using the code or data from this repository:


## Abstract

Genome-wide association studies (GWAS) have identified hundreds of risk loci for breast cancer, but identifying causal variants and candidate target genes remains challenging. Since most risk loci fall in active gene regulatory regions, we developed a novel approach to identify variants with greater regulatory potential in the disease’s tissue of origin. Using genome-wide differential allelic expression (DAE) analysis on microarray data from 64 normal breast tissue samples, we mapped over 54K variants associated with DAE (daeQTLs). We then intersected these with GWAS data to reveal candidate risk regulatory variants and analyzed their cis-acting regulatory potential. We found 122 daeQTLs in 41 loci in active regulatory regions that are in strong linkage disequilibrium with risk-associated variants (risk-daeQTLs). We also identified 65 new candidate target genes in 29 of these loci for which no previous candidates existed. As validation, we identified and functionally characterized five candidate causal variants at the 5q14.1 risk locus targeting the *ATG10* and *ATP6AP1L* genes, likely acting via modulation of alternative transcription and transcription factor binding. Our study demonstrates the power of DAE analysis and daeQTL mapping to understand breast cancer genetic risk, including in complex genetic regulatory landscapes. It additionally provides a genome-wide resource of variants associated with DAE for future functional studies.

## Repository Contents

- `daeQTL_analysis/`: Directory containing all the scripts and code files for reproducing the daeQTL analyses in the paper.
- `dae_analysis/`: Directory containing all the scripts and code files for reproducing the DAE analyses in the paper.

## Environment Setup

This project requires R and several scientific computing libraries indicated at the beginning of the analysis scripts.


## License

MIT

## Contact

For any questions or issues related to the code, please open an issue in this repository or contact the authors directly at atmaia@ualg.pt or jgxavier@ualg.pt.

