---
title: 'nf-gwas-pipeline: A Nextflow Genome-Wide Association Study Pipeline'
tags:
  - R
  - Nextflow
  - GWAS
  - gene-based analysis
  - longitudianl GWAS
authors:
  - name: Zeyuan Song
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Anastasia Gurinovich
    affiliation: 1
  - name: Anthony Federico
    affiliation: "2","4"
  - name: Stefano Monti
    affiliation: "2","4"
  - name: Paola Sebastiani
    affiliation: "3"
affiliations:
 - name: Department of Biostatistics, Boston University School of Public Health, 801 Massachusetts Avenue 3rd Floor, Boston, MA 02218, USA
   index: 1
 - name: Section of Computational Biomedicine, Boston University School of Medicine, 72 East Concord St., Boston, MA 02218, USA
   index: 2
 - name: Institute for Clinical Research and Health Policy Studies, Tufts Medical Center, 800 Washington Street, Boston, MA 02111, USA
   index: 3
 - name: Bioinformatics Program, Boston University, 24 Cummington Mall, Boston, MA 02215, USA
   index: 4
date: 21 November 2020
bibliography: paper.bib

---

# Summary

A tool for conducting Genome-Wide Association Study (GWAS) in a systematic, automated and reproducible manner is overdue. We developed an automated GWAS pipeline by combining multiple analysis tools – including bcftools, vcftools, the R packages SNPRelate/GENESIS/GMMAT and ANNOVAR – through Nextflow, which is a portable, flexible, and reproducible reactive workflow framework for developing pipelines. The GWAS pipeline integrates the steps of data quality control and assessment and genetic association analyses, including analysis of cross-sectional and longitudinal studies with either single variants or gene-based tests, into a unified analysis workflow. The pipeline is implemented in Nextflow, dependencies are distributed through Docker, and the code is publicly available on Github.

# Introduction

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
