---
title: "Matrix Factorization with applications in single-cell RNA-seq"
author: "Yongjin Park"
date: "05 February 2026"
classoption: "aspectratio=169"
fontsize: 12pt
bibliography: singleCellFactorization.bib
output:
    beamer_presentation:
        theme: Madrid
        keep_md: true
        keep_tex: true
        latex_engine: pdflatex
        slide_level: 2
header-includes:
  - \usepackage{cancel}
  - \usepackage{tikz}
  - \usetikzlibrary{fit}
  - \definecolor{UBCblue}{rgb}{0.04706, 0.13725, 0.26667}
  - \usecolortheme[named=UBCblue]{structure}
  - \setbeamertemplate{frametitle}{\color{UBCblue}\bfseries\insertframetitle\par\vskip-6pt\hrulefill}
  - \setbeamercolor{title in head/foot}{bg=white,fg=gray}
  - \setbeamercolor{section in head/foot}{bg=white,fg=gray}
  - \setbeamercolor{author in head/foot}{bg=white,fg=gray}
  - \setbeamercolor{date in head/foot}{bg=white,fg=gray}
  - \setbeamertemplate{page number in head/foot}{}
  - \setbeamertemplate{frame numbering}[none]
  - \setbeamercolor{alerted text}{bg=yellow}
  - \AtBeginSection[]{\begin{frame}\frametitle{Today's lecture}{\Large\tableofcontents[currentsection]}\end{frame}}
  - |
    \makeatletter
    \def\ps@titlepage{%
      \setbeamertemplate{footline}{}
    }
    \addtobeamertemplate{title page}{\thispagestyle{titlepage}}{}
    \makeatother
    \include{toc}
---

<!--
Rules: 

1. After editing this file, run `make` to rebuild the slides: make lect08-single-cell-technology-slides.pdf

2. Use data-beans and senna (legume) Rust-based CLI tools (already installed)

3. Small data sets: we can show both Seurat and Scanpy how we can process

4. Do not mix R and python in the same slide/frame (it's confusing)

5. Prefer hdf5 backend since zarr takes too many files.

6. Prefer onslide rather than only in latex beamer

7. use \large for discussion (text-driven) slides then close it with \normalsize

8. Figure width never exceeds 6 and height never 3

9. Rasterize plots with too many dots with dpi=300 using `ggrastr`

10. use pdflatex preferrably...

-->




## Roadmap of STAT/BIOF/GSAT540 2026

_data generation_ followed by _inference_ or validation

\vfill

\large

:::::: {.columns}
::: {.column width=.45}

3. Formulating biological questions

4. Hypothesis testing

5. RNA-seq data

6. Inference focusing on RNA-seq

7. Experimental design

:::
::: {.column width=.45}

8. Single-cell genomics

9. \textbf{\color{magenta} Matrix factorization} (Today)

10. Genome-wide association studies

11. Multiple hypothesis testing

:::
::::::

\normalsize


## Recap: Overview of single-cell data analysis

\centerline{\includegraphics[width=\textwidth]{img/scRNA_overview/Kharchenko_1.pdf}}

\vfill

\flushright
\small
[@Kharchenko2021-wr]

## Overview of single-cell data analysis cont'd

\centerline{\includegraphics[width=\textwidth]{img/scRNA_overview/Kharchenko_2.pdf}}

\vfill

\flushright
\small
[@Kharchenko2021-wr]


## Learning objectives

- Learn Principal Component Analysis (PCA) and Singular Value Decomposition (SVD)

- Understand matrix factorization and its applications in single-cell analysis

- Introduction to useful tools based on matrix factorization

	- `glmpca`, `rliger`, `fastTopic`, `MOFA`, ...


# What is Principal Component Analysis?




## Let's start with multivariate linear regression

**Linear model for gene $g$:**

$$Y_{ig} = \beta_{1g} \cdot X_{i1} + \beta_{2g} \cdot X_{i2} + \beta_{3g} \cdot X_{i3} + \epsilon_{ig},\, \epsilon_{ig} \sim N(0, \sigma^2_g)$$

:::::: {.columns}
::: {.column width=.4}

\begin{center}
\input{texfrag/graphical-model-linear-regression.tex}
\end{center}

:::
::: {.column width=.45}

* $X_{ij}$: covariate $j$ for sample $i$
* $\beta_{jg}$: effect of covariate $j$ on gene $g$
* $Y_{ig}$: expr. of gene $g$ in sample $i$

:::
::::::

:::
::::::

## A design matrix = a column space view

We can rewrite the linear model as:

$$\mathbf{y}_g = \beta_{0g} \mathbf{x}_0 + \beta_{1g} \mathbf{x}_1 + \beta_{2g} \mathbf{x}_2 + \beta_{3g} \mathbf{x}_3 + \boldsymbol{\epsilon}_g$$

where $\mathbf{x}_j$ is the $j$-th column of $\mathbf{X}$

\onslide<2>{

\textbf{\color{teal} What is a linear regression?}

Finding the "best" \textbf{\color{teal} linear combination} of {\color{magenta} column vectors} to approximate $\mathbf{y}_g$

\begin{itemize}
\item Each column $\mathbf{x}_j$ represents a covariate pattern across all samples
\item The coefficients $\beta_{jg}$ tell us how much each pattern contributes
\item The fitted values $\hat{\mathbf{y}}_g = \mathbf{X}\hat{\boldsymbol{\beta}}_g$ lie in the \textbf{\color{magenta} column space} of $\mathbf{X}$
\end{itemize}

}

## Linear regression = projection onto linear space

If the whole prediction space can be made up by a linear combination of two $n$ vectors $\mathbf{x}_{1}$ and $\mathbf{x}_{2}$...

\begin{center}
\input{texfrag/column-space-projection.tex}
\end{center}

$\hat{\mathbf{y}}_g$ is the **projection** of $\mathbf{y}_g$ onto the column space of $\mathbf{X}$

## Linear regression = projection onto some space

:::::: {.columns}
::: {.column width=.4}

\begin{center}
\input{texfrag/graphical-model-linear-regression-latent.tex}
\end{center}

:::
::: {.column width=.45}

**What if our covariates are unknown?**

$$Y_{ig} = U_{i1} V_{1g} + U_{i2} V_{2g} + X_{i3} V_{3g} + \epsilon_{ig}$$

\vspace{0.3cm}

* $U_{ij}$: **latent** factor $j$ for sample $i$
* $V_{jg}$: loading of gene $g$ on factor $j$
* $Y_{ig}$: observed expression of gene $g$ in sample $i$

:::
::::::

## Many linear regression models on many genes/columns

\begin{center}
\input{texfrag/graphical-model-multi-gene.tex}
\end{center}

## Many linear regression models on many genes/columns

\begin{center}
\input{texfrag/graphical-model-matrix-factorization.tex}
\end{center}


## The goal is to find hidden patterns

:::::: {.columns}
::: {.column width=.5}



 \scriptsize



 \normalsize

A data matrix was generated by the following:

$$\mathbf{x}_{i} = \mathbf{u}_{i} V^{\top} + \epsilon_{i}$$

for all $i$ with $\epsilon_{i} \sim \mathcal{N}\!\left(\mathbf{0},I\right)$


 \scriptsize



 \normalsize


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-8-1} \end{center}


}

 \normalsize

:::
::: {.column width=.5}


 \scriptsize



 \normalsize

Can we reverse engineer the $U$ and $V$ matrices from the data $X$?

$$X \to U V^{\top}$$


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-10-1} \end{center}


}

 \normalsize

:::
::::::


## Let's look at the PDAC 10x Genomics data

\vfill

:::::: {.columns}
::: {.column width=.5}

\small
```bash
PDAC_TISSUE_1/
|-- barcodes.tsv.gz   (9 KB)
|-- features.tsv.gz   (303 KB)
+-- matrix.mtx.gz     (17 MB)
```
\normalsize

\vspace{0.5cm}

* **`barcodes.tsv.gz`** - Cell IDs
* **`features.tsv.gz`** - Genes
* **`matrix.mtx.gz`** - Counts (sparse)

:::
::: {.column width=.45}

**Matrix Market format**

* Sparse matrix standard
* Memory efficient
* Compatible with Seurat, Scanpy

:::
::::::

## Loading 10x data with Seurat (R) \hfill \includegraphics[height=1em]{img/logos/seurat_logo.jpg}


 \scriptsize


```
## An object of class Seurat 
## 19273 features across 1410 samples within 1 assay 
## Active assay: RNA (19273 features, 0 variable features)
##  1 layer present: counts
```

 \normalsize

Extract the count matrix for later use:


 \scriptsize


```
## [1] 19260  1410
```

 \normalsize


## Let's just use top 200 most variable genes and cells


 \scriptsize



 \normalsize

Top 200 genes

:::::: {.columns}
::: {.column width=.25}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-13-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-14-1} \end{center}


}

 \normalsize

:::
::: {.column width=.65}

* We will call such a high-dimensional matrix $X$ ($m \times n$)

* X has $m$=19,260 rows (transcripts/genes/features)

* X has $n$=1,410 columns (cells/#data points)

* The rows were log-transformed and scaled by `scale()` for visualization

* Each cell is a 19,260-dimensional vector!

* Each gene is a 1,410-dimensional vector...

:::
::::::


## 1-dimensional representations/summary of data matrix

:::::: {.columns}
::: {.column width=.25}

Top 200 genes


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-15-1} \end{center}

 \normalsize

:::
::: {.column width=.65}


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-16-1} \end{center}


}

 \normalsize


 \scriptsize

\onslide<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-17-1} \end{center}


}

 \normalsize

:::
::::::


## Can you see some common patterns?


 \scriptsize



 \normalsize

:::::: {.columns}
::: {.column width=.5}


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-19-1} \end{center}

 \normalsize

:::
::: {.column width=.5}

Common patterns will become more apparent after applying PCA.

:::
::::::

## How do we recover "common" patterns from data?

Principal Component Analysis

:::::: {.columns}
::: {.column width=.45}

::: {.block}

### [@Pearson1901-on]

* Projection [of the original data] that minimizes the projection cost between the original and projected

* The cost = mean squared distance

:::

:::
::: {.column width=.45}

::: {.block}

### [@Hotelling1933-aa]

* Orthogonal projection of data into a lower-dimensional **[principal]** sub-space,

* such that the total **variation of the projected data** is maximized

:::

:::
::::::

## Let's see how much the first eigenvector can explain

:::::: {.columns}
::: {.column width=.35}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-20-1} \end{center}


}

 \normalsize


 \scriptsize



 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-22-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-23-1} \end{center}


}

 \normalsize

:::
::: {.column width=.65}


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-24-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-25-1} \end{center}


}

 \normalsize

:::
::::::


## Top two eigen-vectors

:::::: {.columns}
::: {.column width=.35}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-26-1} \end{center}


}

 \normalsize


 \scriptsize



 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-28-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-29-1} \end{center}


}

 \normalsize

:::
::: {.column width=.65}


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-30-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-31-1} \end{center}


}

 \normalsize

:::
::::::

## Top three eigen-vectors

:::::: {.columns}
::: {.column width=.35}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-32-1} \end{center}


}

 \normalsize


 \scriptsize



 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-34-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-35-1} \end{center}


}

 \normalsize

:::
::: {.column width=.65}


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-36-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-37-1} \end{center}


}

 \normalsize

:::
::::::

##

\huge

What is Principal Component Analysis?

\large

(a more mathematical definition)

## PCA: Consider this multivariate linear model

\only<1>{
$$\left(\begin{array}{l}
X_{1}\\
X_{2}\\
\vdots \\
X_{p}\\
\end{array}
\right) =
\left(
\begin{array}{l l l}
U_{11} & \cdots & U_{1k}\\
U_{21} & \cdots & U_{2k}\\
\vdots & \vdots & \vdots \\
U_{p1} & \cdots & U_{pk}\\
\end{array}
\right)
\left(
\begin{array}{l}
V_{1}\\
\vdots\\
V_{k}
\end{array}
\right)
+
\left(
\begin{array}{l}
\epsilon_{1}\\
\epsilon_{2}\\
\vdots\\
\epsilon_{p}
\end{array}
\right)$$

or

$\mathbf{x} = U \mathbf{v} + \epsilon$}

\only<2>{For many columns, $X = U V + E$, or
$$\left(\begin{array}{l l l}
X_{11} & \cdots & X_{1n}\\
X_{21} & \cdots & X_{2n}\\
\vdots & & \vdots \\
X_{p1} & \cdots & X_{pn}\\
\end{array}
\right) =
\left(
\begin{array}{l l l}
U_{11} & \cdots & U_{1k}\\
U_{21} & \cdots & U_{2k}\\
\vdots & \vdots & \vdots \\
U_{p1} & \cdots & U_{pk}\\
\end{array}
\right)
\left(
\begin{array}{l l l}
V_{11} & \cdots & V_{1n}\\
\vdots & & \vdots \\
V_{k1} & \cdots & V_{kn} \\
\end{array}
\right)
+ \cdots$$}

* If we **knew** $U$, we would be able to solve weights $V$.

* And vice versa...

## PCA: What is a projection matrix?

Suppose we knew $V$. For each gene (row) $g$, we solve the following regression.

\begin{eqnarray*}
\mathbf{x}_{g}^{\top} &\sim& V \mathbf{u}_{g}^{\top} + \epsilon I\\
\onslide<2->{\left(\begin{array}{l}
X_{g1} \\
\vdots \\
X_{gn}
\end{array}\right)
&\sim&
\left(\begin{array}{l l l}
V_{11} & \cdots & V_{1k} \\
\vdots & \vdots & \vdots \\
V_{n1} & \cdots & V_{nk} \\
\end{array}\right)
\left(\begin{array}{l}
U_{g1} \\
\vdots \\
U_{gk}
\end{array}\right) +
\left(\begin{array}{l}
\epsilon \\
\vdots \\
\epsilon
\end{array}\right)}
\end{eqnarray*}

Least square solution:

\only<1-2>{$$\underset{\mathbf{u}_{g}}{\min} \, \|\mathbf{x}_{g}^{\top} - V \mathbf{u}_{g}^{\top}\|$$}

\only<3>{$$\hat{\mathbf{u}}_{g}^{\top} = (V^{\top}V)^{-1} V^{\top}\mathbf{x}_{g}^{\top}$$}

## PCA: What is a projection matrix?

Plugging

$$\hat{\mathbf{u}}_{g}^{\top} = (V^{\top}V)^{-1} V^{\top}\mathbf{x}_{g}^{\top}$$

into the following prediction model:

\begin{eqnarray*}
\hat{\mathbf{x}}_{g}^{\top} &=& V \hat{\mathbf{u}}_{g}^{\top} \\
\only<2>{
&=& V (V^{\top}V)^{-1} V^{\top}\mathbf{x}_{g}^{\top} \\ }
\only<3>{
\hat{\mathbf{x}}_{g} &=& \mathbf{x}_{g} \underbrace{V (V^{\top}V)^{-1} V^{\top}}_{\textsf{$n{\times}n$ projection matrix}} }
\end{eqnarray*}

## PCA: Projection can be done on the other side

A multivariate regression model for sample (column) $j$:

$$\mathbf{x}_{j} \sim U \mathbf{v}_{j} + \epsilon$$

The least-square solution for the $V$:

$$\hat{\mathbf{v}}_{j} = (U^{\top}U)^{-1} U^{\top} \mathbf{x}_{j}$$

Then we have

$$\hat{\mathbf{x}}_{j} = \underbrace{U (U^{\top}U)^{-1} U^{\top}}_{\textsf{\color{red} $p{\times}p$ projection matrix}} \mathbf{x}_{j}$$

## PCA: maximizing total variance of the projected

::: {.block}

### Constrained optimization

Given the rank-1 factorization
$$\hat{X} = \mathbf{u} \mathbf{v}^{\top}$$
and
the least square solution $\hat{\mathbf{v}}$
$$\hat{\mathbf{v}}^{\top} = \frac{\mathbf{u}^{\top}}{\mathbf{u}^{\top}\mathbf{u}} X$$

We want to maximize the variance of this projected vector $\mathbf{v}$
$$\hat{\mathbf{v}}^{\top} \hat{\mathbf{v}}
= \frac{\mathbf{u}^{\top}}{\mathbf{u}^{\top}\mathbf{u}} X
X^{\top} \frac{\mathbf{u}}{\mathbf{u}^{\top}\mathbf{u}}
= \mathbf{u}^{\top}XX^{\top}\mathbf{u}$$
where we assume $\mathbf{u}$ is a unit vector.

:::

## PCA is an eigen value problem

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### PCA

Letting the covariance matrix $\hat{\Sigma} = XX^{\top}/(p-1)$, we want to find a unit vector $\mathbf{v}$ by

$$\max \mathbf{v}^{\top} \hat{\Sigma}\mathbf{v}$$

subject to $\mathbf{v}^{\top}\mathbf{v} = 1$.

:::

:::
::: {.column width=.45}

::: {.block}
### Eigen value problem

Given the covariance matrix $\hat{\Sigma}$, we can resolve an eigen-value $\lambda$ and the corresponding eigen-vector $\mathbf{v}$ such that

$$\hat{\Sigma} \mathbf{v} = \lambda \mathbf{v}$$

:::

:::
::::::


## SVD: another equivalent method for PCA

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Singular Value Decomposition

SVD identifies three matrices of $X$:

$$X = U D V^{\top}$$

where both $U$ and $V$ vectors are orthonormal,
namely,

- $U^{\top}U = I$, $\mathbf{u}_{k}^{\top}\mathbf{u}_{k}=1$ for all $k$,

- $V^{\top}V = I$, $\mathbf{v}_{k}^{\top}\mathbf{v}_{k} = 1$ for all $k$.

:::

:::
::: {.column width=.45}

::: {.block}
### Covariance by SVD

Covariance across the columns (samples)

$$X^{\top}X/p = V D^{2} V^{\top}/p$$

Covariance across the rows (genes)

$$XX^{\top}/n = U D^{2} U^{\top}/n$$

:::

*Remark*: The covariance formulas above assume data is centered (mean-subtracted).

:::
::::::

## SVD: another equivalent method for PCA

We can confirm the equivalent relations by multiplying singular vectors to the covariance matrix:

\small

\begin{eqnarray*}
\onslide<1->{
\underbrace{\frac{\left( X^{\top}X \right)}{p-1}}_{\textsf{\color{blue} sample covariance}} \mathbf{v}_{1}
&\propto&
\left(\mathbf{v}_{1}, \mathbf{v}_{2},\ldots, \mathbf{v}_{k}\right)
\left(
\begin{array}{l l l l}
D_{1}^{2} & 0 & \ldots & \ldots \\
0 & D_{2}^{2} & 0 & \ldots \\
0 & \ldots & \ddots & 0 \\
0 & \ldots & 0 & D_{k}^{2} \\
\end{array} \right)
\left(
\begin{array}{l}
\mathbf{v}_{1}^{\top}\\
\mathbf{v}_{2}^{\top}\\
\vdots\\
\mathbf{v}_{k}^{\top}
\end{array}
\right)
\mathbf{v}_{1} \\
}
\only<2>{
&\propto&
\left(\mathbf{v}_{1}, \mathbf{v}_{2},\ldots, \mathbf{v}_{k}\right)
\left(
\begin{array}{l l l l}
D_{1}^{2} & 0 & \ldots & \ldots \\
0 & D_{2}^{2} & 0 & \ldots \\
0 & \ldots & \ddots & 0 \\
0 & \ldots & 0 & D_{k}^{2} \\
\end{array} \right)
\left(
\begin{array}{l}
1 \\
0 \\
\vdots\\
0
\end{array}
\right) \\
}
\only<3>{
&\propto&
\left(\mathbf{v}_{1}, \mathbf{v}_{2},\ldots, \mathbf{v}_{k}\right)
\left(
\begin{array}{l}
D_{1}^{2} \\
0 \\
0 \\
0 \\
\end{array} \right)}
\only<4>{
&=&
\underbrace{\frac{D_{1}^{2}}{p-1}}_{\textsf{\color{red}eigenvalue}}
\underbrace{\mathbf{v}_{1}}_{\textsf{\color{red}eigenvector}}
}
\end{eqnarray*}

## Run SVD to find principal components


 \scriptsize


``` r
svd.out <- rsvd::rsvd(x.sorted, k = 7)
U <- svd.out$u; D <- diag(svd.out$d); V <- svd.out$v
```

 \normalsize

:::::: {.columns}
::: {.column width=.25}


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-39-1} \end{center}

 \normalsize

:::
::: {.column width=.65}


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/unnamed-chunk-40-1} \end{center}

 \normalsize

:::
::::::


## Summary: PCA

\large

* Matrix factorization = alternating linear regressions 

* PCA is a matrix factorization

* PCA can be found by Eigen value/vector problem

* PCA can be found by Singular Value Decomposition

\normalsize



## Generalized Linear Model PCA to model count data

:::::: {.columns}
::: {.column width=.45}

\large

- scRNA-seq data are **counts** (discrete, non-negative)

$$\log(\mu_{ij}) = (UV^{\top})_{ij}$$

$$Y_{ij} \sim \text{Poisson}(\mu_{ij})$$

or

$$Y_{ij} \sim \text{NB}(\mu_{ij}, \phi)$$

\normalsize

:::
::: {.column width=.45}


 \scriptsize


``` r
library(glmpca)

## Run GLM-PCA on count matrix
res <- glmpca(X, L = 5, fam = "nb")

factors <- res$factors
loadings <- res$loadings
```

 \normalsize

:::
::::::

## GLM-PCA vs standard PCA

\large

**When to use GLM-PCA?**

\normalsize

- Working with **raw counts** (not normalized)
- Want to avoid log-transformation artifacts
- More computationally intensive than standard PCA

\vspace{0.5cm}

**Reference:** [@Townes2019-rg]


# Non-negative Matrix Factorization



## Why Non-negative Matrix Factorization (NMF)?

:::::: {.columns}
::: {.column width=.45}

**PCA (SVD) can have negative values:**

E.g., For some gene $g$:
We can have
\begin{eqnarray*}
U_{g1} + U_{g2} &=& 0 \\
\text{if } U_{g1} &=& - U_{g2} < 0 \\
\text{or } U_{g1} &=& U_{g2} = 0 \\
\vdots
\end{eqnarray*}

:::
::: {.column width=.45}

**NMF constraint:**

\large

$$X \approx UV$$ 

where $U, V \geq 0$

\normalsize

:::
::::::

## Classic NMF in Genomics (metagenes)

:::::: {.columns}
::: {.column width=.45}

\large

**@Brunet2004-fk:**

\normalsize

- Discovers **"metagenes"** (coherent patterns)
- Non-negative gene/sample combinations
- Only additive combinations (vs. PCA's +/- components)

:::
::: {.column width=.45}

- Metagenes = pathways/processes
- Genes are "on" or "off" (not negative)
- Sparse representation
- Natural for expression data

\vspace{0.3cm}

**Applications:**

- Cancer subtype discovery
- Gene expression pattern analysis
- Identifying co-regulated gene modules

:::
::::::

\vfill

\flushright
\small
@Brunet2004-fk

## 

<!---------------------------->
<!-- will add some figures  -->
<!---------------------------->

## NMF for single-cell RNA-seq

\large

**Why NMF is particularly useful for scRNA-seq:**

\normalsize

- Gene expression is naturally non-negative (counts)
- Identifies **gene programs** (co-expressed gene sets)
- Discovers **cell states** and gradients
- Integrates datasets across batches/conditions
- More interpretable than PCA for biological processes

\vspace{0.5cm}

**Two popular tools:** `rliger` and `fasttopic`



## `rliger`: Integrative NMF for single-cell data

:::::: {.columns}
::: {.column width=.25}

\includegraphics[width=\textwidth]{img/logos/rliger_logo.png}

:::
::: {.column width=.6}

- Joint NMF across multiple datasets
- Shared factors ($W$) capture common biology
- Dataset-specific factors ($V_i$) capture batch effects
- Quantile normalization aligns factor loadings

$$X^{(d)} \approx H^{(d)} (W + V^{(d)})^{\top}$$

$H^{(d)}$: cell factor loadings for dataset $d$

$W$: shared gene factors, $V^{(d)}$: dataset-specific

:::
::::::


@Welch2019-jl; @Gao2021-vy; @Kriebel2022-fe 


## `rliger`: Running integrative NMF


 \scriptsize



 \normalsize





 \scriptsize


``` r
library(rliger)

## Load three PDAC tissue samples
tissues <- c("PDAC_TISSUE_1", "PDAC_TISSUE_2", "PDAC_TISSUE_3")

data.list <- lapply(tissues, function(tissue) {
  data <- Read10X(data.dir = paste0("../data/GSE155698_PDAC_Steele2020/", tissue, "/"))
  ## Remove MT-genes
  mt.genes <- grep("^MT-", rownames(data), value = TRUE)
  data[!rownames(data) %in% mt.genes, ]
})
names(data.list) <- tissues
```

 \normalsize

## `rliger`: Create object and preprocess


 \scriptsize


``` r
## Create liger object with three datasets
liger_obj <- createLiger(data.list)

## Normalize and select variable genes
liger_obj <- normalize(liger_obj)
liger_obj <- selectGenes(liger_obj)
liger_obj <- scaleNotCenter(liger_obj)
```

 \normalsize

## `rliger`: Run integrative NMF


 \Large


``` r
## Run integrative NMF
liger_obj <- optimizeALS(liger_obj, k = 20)
```

 \normalsize

## `rliger`: Align datasets and cluster


 \scriptsize


``` r
## Quantile normalization for alignment
liger_obj <- quantile_norm(liger_obj)

## Clustering on aligned factors
liger_obj <- louvainCluster(liger_obj, resolution = 0.4)

## UMAP visualization
liger_obj <- runUMAP(liger_obj, n_neighbors = 30)
```

 \normalsize

## `rliger`: UMAP by dataset


 \scriptsize


``` r
## UMAP colored by dataset
plotDatasetDimRed(liger_obj)
```



\begin{center}\includegraphics{../Fig/single-cell-factorization/liger-umap-dataset-1} \end{center}

 \normalsize

## `rliger`: UMAP by cluster


 \scriptsize


``` r
## UMAP colored by cluster
plotClusterDimRed(liger_obj)
```



\begin{center}\includegraphics{../Fig/single-cell-factorization/liger-umap-cluster-1} \end{center}

 \normalsize

## `fastTopics`: Topic modeling for scRNA-seq

\large

- Topic model: each cell $F$ = **mixture of topics**
- each topic $L$ = distribution over genes

$$Y_{gj} = \sum_{k} L_{gk} F_{kj}$$

\normalsize

@Carbonetto2021-ft; @Carbonetto2023-gd




 \scriptsize


``` r
library(fastTopics)

## Use the same PDAC tissue 1 count matrix (cells x genes)
## Select top 5000 most variable genes
gene_var <- apply(X, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]
counts <- t(X[top_genes, ])
dim(counts)
```

 \normalsize

## `fastTopics`: Fit topic model


 \Large


``` r
## Fit topic model with k=10 topics
fit <- fit_topic_model(counts, k = 10)
```

 \normalsize

## `fastTopics`: Topic proportions per cell


 \scriptsize


``` r
## L matrix: cells x topics (topic proportions)
head(round(fit$L, 3))
```

```
##                    k1    k2    k3    k4    k5   k6    k7    k8    k9   k10
## AAACGAAAGTGGAAAG-1  0 0.000 0.000 0.006 0.000 0.00 0.000 0.407 0.000 0.586
## AAACGAAGTAGGGTAC-1  0 0.000 0.000 0.002 0.000 0.55 0.000 0.004 0.292 0.152
## AAACGAAGTCATAGTC-1  0 0.028 0.000 0.523 0.000 0.00 0.161 0.000 0.225 0.062
## AAAGAACCATTAAAGG-1  0 0.114 0.761 0.000 0.011 0.00 0.000 0.050 0.023 0.040
## AAAGGATTCGGCTTGG-1  0 0.233 0.000 0.003 0.635 0.00 0.000 0.000 0.000 0.128
## AAAGGGCAGTAGCAAT-1  0 0.140 0.650 0.001 0.000 0.00 0.009 0.182 0.005 0.014
```

 \normalsize

## `fastTopics`: Top genes per topic


 \scriptsize


``` r
## F matrix: genes x topics (gene weights)
top_genes <- apply(fit$F, 2, function(x) {
  names(sort(x, decreasing = TRUE)[1:10])
})
top_genes[, 1:3]
```

```
##       k1       k2       k3      
##  [1,] "PRSS1"  "MALAT1" "FTL"   
##  [2,] "REG3A"  "RPS27"  "SPP1"  
##  [3,] "CTRB1"  "RPL41"  "FTH1"  
##  [4,] "PNLIP"  "RPL10"  "B2M"   
##  [5,] "REG1A"  "RPL13"  "CTSD"  
##  [6,] "CPB1"   "EEF1A1" "TMSB4X"
##  [7,] "CPA1"   "TPT1"   "APOE"  
##  [8,] "CLPS"   "B2M"    "HLA-B" 
##  [9,] "CTRB2"  "RPS18"  "CTSB"  
## [10,] "CELA3A" "RPS12"  "ACTB"
```

 \normalsize

## `fastTopics`: Structure plot


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-structure-1} \end{center}

 \normalsize

## `fastTopics`: Volcano plot — genes driving each topic

\only<1>{Topic 1: PRSS1, REG3A, CTRB1 — pancreatic acinar cells}
\only<2>{Topic 2: MALAT1, RPL41, RPL10 — housekeeping/ribosomal}
\only<3>{Topic 3: FTL, SPP1, CTSD — macrophage/myeloid markers}
\only<4>{Topic 4: RARRES1, PLA2G2A, RPL18A — ductal/regenerating}
\only<5>{Topic 5: TMSB4X, GNLY, GZMA — cytotoxic T/NK cells}
\only<6>{Topic 6: RARRES1, FTH1, HLA-B — ductal/epithelial}
\only<7>{Topic 7: S100A6, RPLP1, RPS19 — epithelial}
\only<8>{Topic 8: S100A9, S100A8, LYZ — neutrophils/monocytes}
\only<9>{Topic 9: S100A6, FTH1, TMSB10 — proliferating}
\only<10>{Topic 10: MALAT1, NEAT1, HLA-A — antigen presentation/stress}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-1-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-2-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-3-1} \end{center}


}

 \normalsize


 \scriptsize

\only<4>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-4-1} \end{center}


}

 \normalsize


 \scriptsize

\only<5>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-5-1} \end{center}


}

 \normalsize


 \scriptsize

\only<6>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-6-1} \end{center}


}

 \normalsize


 \scriptsize

\only<7>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-7-1} \end{center}


}

 \normalsize


 \scriptsize

\only<8>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-8-1} \end{center}


}

 \normalsize


 \scriptsize

\only<9>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-9-1} \end{center}


}

 \normalsize


 \scriptsize

\only<10>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/fasttopic-volcano-10-1} \end{center}


}

 \normalsize


# How can we use factorization results?

## Gene set enrichment analysis with PanglaoDB




 \scriptsize


``` r
library(fgsea)

## Read PanglaoDB cell-type marker gene sets
panglao.db <- read.panglao("Pancreas")

## Prepare gene weights from fastTopics F matrix
## Each column = topic, rows = genes
.beta <- apply(fit$F, 2, scale)
rownames(.beta) <- colnames(counts)

## Sort genes by topic assignment
.genes <- unique(unlist(panglao.db))
.dt <- sort.beta(.beta, genes.selected = .genes)

## Run fgsea per topic
.gsea <- run.beta.gsea(.dt, panglao.db)
```

 \normalsize

## GSEA results: which topics correspond to cell types?


 \scriptsize


\begin{center}\includegraphics{../Fig/single-cell-factorization/gsea-plot-1} \end{center}

 \normalsize

## GSEA: Gene weights across topics


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/gsea-beta-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/gsea-combined-1} \end{center}


}

 \normalsize

## `rliger`: GSEA on shared gene loadings (W matrix)




 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/liger-gsea-plot-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/single-cell-factorization/liger-gsea-combined-1} \end{center}


}

 \normalsize

## Wrap-up

\large

- **PCA/SVD**: Find orthogonal axes of maximum variance; fast, good for QC
- **NMF**: Non-negative, parts-based decomposition; interpretable gene programs
- **`rliger`**: Integrative NMF across multiple datasets (batch correction)
- **`fastTopics`**: Scalable topic model; cells as mixtures of gene programs
- **GSEA on factors**: Link learned topics to known cell types (e.g., PanglaoDB + `fgsea`)

\normalsize

##

## Reference {.allowframebreaks}
