---
title: "Regression for Regulatory Genomics"
author: "Yongjin Park"
date: "26 February 2026"
classoption: "aspectratio=169"
fontsize: 12pt
bibliography: regulatory-genomics.bib
output:
    beamer_presentation:
        theme: Madrid
        keep_md: true
        keep_tex: true
        latex_engine: pdflatex
        slide_level: 2
header-includes:
  - \usepackage{cancel}
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




## Learning objectives

By the end of this lecture, you will be able to:

\large

1. Apply regression methods to genomic data integration problems

2. Understand fine-mapping and variable selection in regulatory genomics

3. Integrate GWAS and regulatory data using TWAS and colocalization

## Roadmap: From biology to methodology

\large

:::::: {.columns}
::: {.column width=.48}

**Lecture 12 (Biology)**

* Regulatory genomics data types

* ChIP-seq, ATAC-seq, DNA methylation

* Biological interpretation

:::
::: {.column width=.48}

**Today (Methodology)**

* Polygenic score model

* eQTLs and chromatin accessibility

* Regression for genomic integration

* Integrative analysis methods

:::
::::::

## We will discuss multivariate linear regression




## Genetic data as an anchor variable

\large

$$\text{Genetic variants} \to \text{Molecular Phenotype} \to \text{Disease Phenotype}$$


\normalsize


## Using 1000 Genomes data via `bigsnpr`


 \large


``` r
library(bigsnpr)
download_1000G("../data/genotype/")
```

 \normalsize

\vfill


 \scriptsize



 \normalsize


 \scriptsize


``` r
data <- snp_attach(.bk.file)
str(data, max.level=1, strict.width = "cut")
```

```
## List of 3
##  $ genotypes:Reference class 'FBM.code256' [package "bigstatsr"] with 16 fields
##   ..and 26 methods, of which 12 are  possibly relevant
##  $ fam      :'data.frame':	2490 obs. of  6 variables:
##  $ map      :'data.frame':	1664852 obs. of  6 variables:
##  - attr(*, "class")= chr "bigSNP"
```

 \normalsize







 \scriptsize



 \normalsize


 \scriptsize



 \normalsize



# Polygenic risk prediction (G $\to$ phenotype)

## Recap: Notations

\large

* $X$: ${n} \times {p}$ genotype matrix of ${n}$ samples/individuals and ${p}$ variants

* $Y$: ${n} \times {1}$ phenotype vector of ${n}$ samples/individuals

## Let's do GWAS to find the disease-associated variants

A variant-by-variant association t-test for a variant $j$:

\begin{eqnarray*}
H_{0}:
&\,& 
\mathbb{E}\!\left[X_{ij}|Y_{i}=1\right] = \mathbb{E}\!\left[X_{ij}|Y_{i}=0\right] \\
& & \textsf{vs.} \\
\, H_{1}:&\,& \mathbb{E}\!\left[X_{ij}|Y_{i}=1\right] \neq \mathbb{E}\!\left[X_{ij}|Y_{i}=0\right]
\end{eqnarray*}

* The average genotype is the same between the case and control under the null.

<!-- ## For one variant's test -->
<!-- ```{r echo=T, size="large"} -->
<!-- j <- 1 -->
<!-- t.test(X[Y == 0, j], X[Y == 1, j]) -->
<!-- ``` -->

## GWAS statistics over all the variants


 \large


``` r
## library(matrixTests)
.gwas <- col_t_welch(X[Y == 0, ], X[Y == 1, ])
```

 \normalsize

* For each pair of columns in the two matrices, the function performs a t-test and summarizes all the results into a list of vectors.


 \scriptsize


``` r
names(.gwas)
```

```
##  [1] "obs.x"       "obs.y"       "obs.tot"     "mean.x"      "mean.y"     
##  [6] "mean.diff"   "var.x"       "var.y"       "stderr"      "df"         
## [11] "statistic"   "pvalue"      "conf.low"    "conf.high"   "mean.null"  
## [16] "alternative" "conf.level"
```

 \normalsize


## Human genetic variation recap

\centerline{\includegraphics[height=.7\textheight]{img/altshuler_daly_lander_GWAS.pdf}}

\vfill

[@Altshuler2008-wr]

## Human genetic variation recap

\centerline{\includegraphics[height=.8\textheight]{img/altshuler_daly_lander_LD.pdf}\vspace{-16pt}}

[@Altshuler2008-wr]

## Manhattan plot: a summary of all the GWAS p-values


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-9-1} \end{center}

 \normalsize


## The covariance of $X$ between the variants (columns)


 \scriptsize



 \normalsize

:::::: {.columns}
::: {.column width=.5}


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-11-1} \end{center}

 \normalsize

:::
::: {.column width=.5}


For a genotype matrix $X$ ${(n{\times}p)}$,



::: {.block}
### Linkage disequilibrium matrix

A ${p\times{p}}$ matrix:

$$\hat{R} = \frac{1}{n}X^{\top}X$$

column $\mathbf{x}_{j}$ standardized (mean 0, SD 1).

:::

\centerline{$R_{ij} = \frac{1}{n} \sum_{r=1}^{n} X_{ri} X_{rj}$}

:::
::::::

## The covariance of $X$ between the individuals (rows)


 \scriptsize



 \normalsize

:::::: {.columns}
::: {.column width=.5}


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-13-1} \end{center}

 \normalsize

:::
::: {.column width=.5}

For a genotype matrix $X$ ${(n{\times}p)}$


::: {.block}
### Genetic relatedness matrix

An ${n\times{n}}$ matrix:

$$\hat{K} = \frac{1}{p}XX^{\top}$$

each row $\mathbf{x}_{i}$ standardized (mean 0, SD 1).

:::

\centerline{$K_{ij} = \frac{1}{p} \sum_{g=1}^{p} X_{ig} X_{jg}$}

:::
::::::



## What's next after GWAS?

\Large

- GWAS until 2010s: heavy focuses on **mapping**

    - GWAS map: genetic variants $\to$ a phenotype

    - Stringent genome-wide p-value cutoff

    - Study design, meta analysis

- NHGRI-EBI GWAS Catalog: [`https://www.ebi.ac.uk/gwas/`](https://www.ebi.ac.uk/gwas/)

##

\huge

Let's take a look at the GWAS Catalog:

[`https://www.ebi.ac.uk/gwas/`](https://www.ebi.ac.uk/gwas/)


## What's next after GWAS?

\Large

- GWAS since 2010s: more emphasis on **prediction**

	- Can we turn GWAS results to a prediction model?

	- Can we understand the mechanisms?

    - Machine learning, data integration, causal inference



## Polygenic score models predict disease prevalence based on genome

- Poly + genic = many genes (units of hereditary components)

\vfill

::: {.block}
### A (linear) polygenic score (PGS)

$$Y_{i} = \sum_{j} X_{ij} \beta_{j}$$

for an individual $i$.

:::

\onslide<2>{Knowing $\beta$'s, we can predict:
$${\color{teal}Y^{\star}}_{i} \gets \sum_{j} {\color{magenta} X_{ij}^{\star}} \hat{\boldsymbol{\beta}}_{j}$$}



## Example: PGS for breast cancer occurrence

:::::: {.columns}
::: {.column width=.7}

\includegraphics[width=\linewidth]{img/PGS_Yang_Pharoah_NR_Cancer_2023_BRCA.pdf}

:::
::: {.column width=.3}

\normalsize

* Strong heritability (proportion of disease risk variation explained by genetics): 35% - 80%

* Prediction is more accurate than lifestyle and other environmental factors

* More powerful if combined with rare risk factors

:::
::::::

\vfill

\tiny

[@Yang2023-id]

## Example: Coronary artery disorder $\to$ prevention

\centerline{\includegraphics[width=.8\linewidth]{img/PGS_Khera_Kathiresan_NatGen_2018_CAD.pdf}}

* Top .5% of PGS values $\to$ five fold increase of CAD

\vfill

\tiny

[@Khera2018-cad]

## Example: PGS stratifies individuals' disease onset

\centerline{\includegraphics[width=.8\linewidth]{img/PGS_Mars_Ripatti_Nature_Medicine_2020_1.pdf}}

\vfill

* PRS $\approx$ PGS.

* We can partition cohorts based on PGS profiles.

\tiny

[@Mars2020-vw]


## GWAS for "Obsessive `ggplot` Disorder"^[just for illustration purposes]




 \scriptsize



 \normalsize


 \scriptsize



 \normalsize

A genotype matrix $X$ ($X_{ij} \in \{0,1,2\}$):


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-17-1} \end{center}

 \normalsize

## In GWAS, we can have case vs. control phenotype

\large

$Y_{i} = 1$ if case vs. $Y_{i} = 0$ if control

For our *OGD* GWAS:


 \scriptsize



 \normalsize

:::::: {.columns}
::: {.column width=.2}

:::
::: {.column width=.8}



 \large


``` r
table(y)
```

```
## y
##   0   1 
## 731 769
```

 \normalsize

:::
::::::


- E.g., Cancer vs. non-cancer, schizophrenia vs. no mental disorder

- The fundamental question is, "Which individuals/subjects can be truly labelled control/wild type?"

## Under the hood, there are quantitative "risk" scores


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-20-1} \end{center}

 \normalsize

- E.g., height, BMI, parents' age of death, how many packs of cigarettes, etc.

## Side note: Almost all complex traits are polygenic

Polygenic effects of common SNPs $\to$ polymorphism in phenotypes

\vfill

> A polygenic trait is a characteristic, such as height or skin color, that is influenced by two or more genes. **Because multiple genes are involved, polygenic traits do not follow the patterns of Mendelian inheritance.** Many polygenic traits are also influenced by the environment and are called multifactorial.

\vfill

https://www.genome.gov/genetics-glossary/Polygenic-Trait


## A polygenic score to predict disease prevalence


 \scriptsize



 \normalsize


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-22-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-23-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-24-1} \end{center}


}

 \normalsize

## Polygenic score $\propto$ the odds of the case vs. control

Log-odds ratio:

$$g(y) = \log\frac{p(\textsf{disease}|Y > y)}{p(\textsf{no disease}|Y > y)}$$

\vfill


 \scriptsize

\onslide<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-25-1} \end{center}


}

 \normalsize

## A PGS model to explain disease prevalence

Log-odds ratio:

$$g(y) = \log\frac{p(\textsf{disease}|Y > y)}{p(\textsf{no disease}|Y > y)}$$

\vfill

**Goal**: Estimate a function to predict this log-odds ratio.

\only<1>{$$g(y_{i}) \sim \beta_{1} X_{i1} + \beta_{2} X_{i2} \cdots + \epsilon$$}

\only<2>{$$g(y_{i}) \sim \underbrace{\beta_{1} X_{i1} + \beta_{2} X_{i2} \cdots }_{\textsf{\color{teal}\bf genetic effects}} + \epsilon$$}


## PGS estimation is a linear modelling with ${p \gg n}$

\only<1>{$$Y_{i} \sim X_{i1} \beta_{1} + X_{i2} \beta_{2} + \cdots$$}

\only<2>{$$\left(
\begin{array}{l}
Y_{1}\\
Y_{2}\\
\vdots\\
Y_{n}
\end{array}
\right)
\sim
\left(
\begin{array}{l}
X_{11}\\
X_{21}\\
\vdots\\
X_{n1}
\end{array}
\right)
\beta_{1} +
\left(
\begin{array}{l}
X_{12}\\
X_{22}\\
\vdots\\
X_{n2}
\end{array}
\right)
\beta_{2} +
\cdots$$}

\onslide<3-5>{$$\left(
\begin{array}{l}
Y_{1}\\
Y_{2}\\
\vdots\\
Y_{n}
\end{array}
\right)
\sim
\left(
\begin{array}{l}
X_{11}\\
X_{21}\\
\vdots\\
X_{n1}
\end{array}
\right)
\beta_{1} +
\left(
\begin{array}{l}
X_{12}\\
X_{22}\\
\vdots\\
X_{n2}
\end{array}
\right)
\beta_{2} +
\cdots +
\left(
\begin{array}{l}
X_{1p}\\
X_{2p}\\
\vdots\\
X_{np}
\end{array}
\right)
\beta_{p}
$$
}

\onslide<4-5>{
$$\mathbf{y} \sim \mathbf{x}_{1} \beta_{1} + \mathbf{x}_{2} \beta_{2} + \cdots + \mathbf{x}_{p} \beta_{p}$$

$$p \gg n$$
}

\only<5>{\large
$n \approx 10^{4}$
but
$p \approx 10^{6}$ for many large-scale GWAS}

## $p \gg n$: Why can't we just fit a model?


 \large


``` r
## lm.out <- lm(Y ~ X)
```

 \normalsize

\large

\begin{itemize}[<+-|alert@+>]
\item $p \gg n$: need to estimate $\beta_{1},\ldots,\beta_{p}$
\item How many samples? How many unknowns?
\item Is it computationally feasible?
\end{itemize}

## $p \gg n$: Variant-by-variant GWAS is not a bad idea


 \large


``` r
.gwas <- col_t_welch(X[Y == 0, ], X[Y == 1, ]) # univar.
```

 \normalsize


 \scriptsize

\onslide<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-28-1} \end{center}


}

 \normalsize

## A linear model including only these "GWAS" variants

:::::: {.columns}
::: {.column width=.45}


 \scriptsize


``` r
X.gwas <- X[, p.adjust(.gwas$pvalue) < .05]
head(X.gwas, 3)
```

```
##         41270163 41270937 41272533 41274422
## HG00159        0        0        0        0
## HG01495        1        1        1        1
## HG01519        2        2        2        2
```

 \normalsize

:::
::: {.column width=.45}


 \scriptsize



 \normalsize


 \scriptsize


``` r
lm.out <- lm(Y ~ X.gwas)
```

 \normalsize

:::
::::::

:::::: {.columns}
::: {.column width=.2}

:::
::: {.column width=.8}


 \scriptsize


``` r
anova(lm.out)
```

```
## Analysis of Variance Table
## 
## Response: Y
##             Df  Sum Sq Mean Sq F value    Pr(>F)    
## X.gwas       4  212.77  53.194  63.965 < 2.2e-16 ***
## Residuals 1495 1243.25   0.832                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

 \normalsize

:::
::::::

##

:::::: {.columns}
::: {.column width=.5}


 \scriptsize

\onslide<1-2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-33-1} \end{center}


}

 \normalsize
:::
::: {.column width=.5}


 \scriptsize

\onslide<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-34-1} \end{center}


}

 \normalsize

:::
::::::


## Do we need all of these GWAS variants?


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-35-1} \end{center}

 \normalsize

## Just one of the GWAS variants works fine


 \scriptsize



 \normalsize

:::::: {.columns}
::: {.column width=.5}


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-37-1} \end{center}

 \normalsize
:::
::: {.column width=.5}


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-38-1} \end{center}

 \normalsize

:::
::::::

## Simply add more variables?


 \scriptsize



 \normalsize

Let's do some experiments:

\large

1. Pick top *K* variants (ordering by GWAS, $p_{(1)} < p_{(2)} < \ldots$)

2. Take GWAS effect sizes $\beta$ (the mean difference between case vs. control)

3. Predict by a linear combination of these effects:

$$\sum_{j=1}^{K} X_{i(j)} \beta_{(j)}$$

## PGS is a variable selection problem

:::::: {.columns}
::: {.column width=.5}


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-40-1} \end{center}


}

 \normalsize

 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-41-1} \end{center}


}

 \normalsize

 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-42-1} \end{center}


}

 \normalsize

 \scriptsize

\only<4>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-43-1} \end{center}


}

 \normalsize

 \scriptsize

\only<5>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-44-1} \end{center}


}

 \normalsize

:::
::: {.column width=.5}


 \scriptsize

\onslide<5->{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-45-1} \end{center}


}

 \normalsize

:::
::::::


# Fundamental challenges in PGS (regression)

## GWAS statistics only roughly tag causal variants

\centerline{\includegraphics[width=.8\linewidth]{img/direct_indirect_associations.pdf}}

\vfill

\tiny

[@Balding2006-bn]

## Covariance between variants obfuscates GWAS interpretation

:::::: {.columns}
::: {.column width=.3}

* Linkage disequilibrium (LD) block

* The result of recombination events throughout generations

:::
::: {.column width=.7}

\includegraphics[height=.7\textheight]{img/altshuler_daly_lander_LD.pdf}

:::
::::::


## A typical setting of PGS training

::: {.block}
### PGS training

* Input:

    1. GWAS summary statistics (p-values, univariate effect sizes)

	2. $n \times p$ genotype matrix $X$, where $X_{ij}\in\{0,1,2\}$

	3. $n \times 1$ phenotype vector $\mathbf{y}$

* Output:

    - *Polygenic* effect sizes $\hat{\beta}_{1},\ldots,\hat{\beta}_{p}$

* Objective:

    - We want to predict unseen $Y^{\star}_{i}$ by $f(X^{\star}; \hat{\beta})$ as accurately as possible

	- We want $\beta_{j} = 0$ as many as possible; only some $\beta_{j}{\neq}0$.

:::


## Common strategies to deal with LD structures

:::::: {.columns}
::: {.column width=.3}

Strategy 1. Just use them all

:::
::: {.column width=.3}

Strategy 2. Pruning/clumping

:::
::: {.column width=.3}

Strategy 3. **Fine-mapping**

:::
::::::

\vfill

:::::: {.columns}
::: {.column width=.3}
\includegraphics[width=\linewidth]{img/GWAS_LD_strategy_useall.pdf}
:::
::: {.column width=.3}
\includegraphics[width=\linewidth]{img/GWAS_LD_strategy_pruning.pdf}
:::
::: {.column width=.3}
\includegraphics[width=\linewidth]{img/GWAS_LD_strategy_fm.pdf}
:::
::::::


## Causal variables explain a large fraction of var.


 \scriptsize


``` r
.lm <- lm(.sim$y.q ~ .sim$x[, .sim$causal, drop = FALSE] - 1)
```

 \normalsize


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-47-1} \end{center}

 \normalsize

## A general workflow of PGS estimation

\only<1>{\centerline{\includegraphics[height=.7\textheight]{img/PGS_Choi_Protocol_1.pdf}}}

\only<2>{\centerline{\includegraphics[height=.7\textheight]{img/PGS_Choi_Protocol_2.pdf}}}

\vfill

\small

[@Choi2020-rn]


## A proposed workflow of PGS estimation

\large

1. Take GWAS summary statistics and target population genotype $X$

2. Run SuSiE^[SuSiE can run with summary statistics only] chromosome by chromosome to select variables

3. Predict PGS per chromosome and aggregate them

## PGS Catalog: No need to train from scratch

\large

`https://www.pgscatalog.org/`

\vfill

\centerline{\includegraphics[height=.6\textheight]{img/PGS_catalog.png}}


# Expression Quantitative Trait Loci (G $\to$ Expr.)


## Transcriptome-wide association study (TWAS)

Using eQTL, we can understand GWAS mechanisms... [@Gamazon2015-pc]

\centerline{\includegraphics[width=\textwidth]{img/GWAS_blackbox_TWAS.pdf}}

\vfill

$$\mathbf{m}_{g} \sim X \boldsymbol{\alpha}_{g} + \boldsymbol{\epsilon}_{g} \quad \implies \quad \mathbf{y}_{\textsf{GWAS}} = \sum_{g \in \textsf{causal genes}} \mathbf{m}_{g} \beta_{g} + \boldsymbol{\epsilon}_{y}$$

## eQTL mapping as a linear regression problem

A variant-by-variant association test for a variant $j$:

\begin{eqnarray*}
H_{0}:
&\,&
\mathbb{E}\!\left[Y_{i}|X_{ij}\right] = \mu \\
& & \textsf{vs.} \\
\, H_{1}:&\,& \mathbb{E}\!\left[Y_{i}|X_{ij}\right] = \mu + X_{ij}\beta_{j}
\end{eqnarray*}

* For gene expression $Y$, test association with each variant $X_j$

* Model: $Y_{i} = X_{ij}\beta_{j} + \epsilon$

* Repeat for all gene-variant pairs in cis window (e.g., $\pm$ 1Mb)

## The high-dimensional challenge: $p \gg n$

Predicting expression from genotypes is a high-dimensional problem:

$$Y_{i} = \sum_{j=1}^{p} X_{ij} \theta_{j} + \epsilon$$

\small
where $\theta_{j}$ are the **multivariate/joint effects** (cf. $\beta_{j}$ for univariate GWAS effects)

\large

\begin{itemize}[<+-|alert@+>]
\item $p \gg n$: need to estimate $\theta_{1},\ldots,\theta_{p}$
\item How many samples? How many unknowns?
\item Is it computationally feasible?
\item Standard linear regression fails!
\end{itemize}





## Penalized regression: Motivation

\large

* Many predictors (variants), few samples

* Need to select relevant variants

* Avoid overfitting

* Two main approaches:
  - Ridge (L2): shrinkage without variable selection
  - LASSO (L1): shrinkage + sparse selection
  - Elastic Net: combination of both

## Application: PrediXcan framework

\Large

**Goal**: Predict gene expression from genotypes to understand GWAS

\normalsize

**Two-stage approach**:

1. **Training**: Build expression prediction model
   $$\hat{E}_g = \sum_{j} X_{j} \beta_{j}^{(g)}$$
   using LASSO/Elastic Net on eQTL reference panel

2. **Prediction**: Test predicted expression with trait
   $$Y \sim \hat{E}_g + \textsf{covariates}$$

\vfill

\tiny
Gamazon et al. (2015) Nature Genetics

## Lasso: a linear regression with L1 [@Tibshirani1996-hc]

Prior distribution
\centerline{$p(\boldsymbol{\theta}) = \textsf{Laplace}(\boldsymbol{\theta}| \lambda) \propto \exp\left(-\lambda\|\boldsymbol{\theta}\|_{1}\right)$}
where
\large
$$\|\boldsymbol{\theta}\|_{1} = \sum_{j=1}^{p} |\theta_{j}|,\,\textsf{\color{blue}L1-norm}.$$

\normalsize
Maximize
$$
\ln p(\mathbf{y}|X,\boldsymbol{\theta}) + \ln p(\boldsymbol{\theta}|\lambda)
= - \frac{1}{2\sigma^{2}} \sum_{i=1}^{n} (y_{i} - \mathbf{x}_{i} \boldsymbol{\theta})^{2}
- \lambda \|\boldsymbol{\theta}\|_{1}
$$

Minimize $L_{1}$-regularized error
$$
\sum_{i=1}^{n} (y_{i} - \mathbf{x}_{i} \boldsymbol{\theta})^{2}
+ \lambda \|\boldsymbol{\theta}\|_{1}
$$


## `glmnet` [@Friedman2010-kx] solves this opt. problem 

Goal (by variable-by-variable updates):\centerline{$\min_{\boldsymbol{\theta}} \quad
\overbrace{(\mathbf{y} - X\boldsymbol{\theta})^{\top}(\mathbf{y} - X\boldsymbol{\theta})}^{\textsf{\color{blue} RSS}} + \underbrace{\lambda \alpha \|\boldsymbol{\theta}\|_{1}}_{\textsf{\color{red} variable selection}} + \underbrace{\lambda (1 - \alpha) \|\boldsymbol{\theta}\|_{2}}_{\textsf{\color{magenta} shrinkage}}$}

\onslide<2->{
For each $\theta_{j}$,
}

\only<2>{$$
\hat{\theta}_{j}^{\textsf{glmnet}} \gets
\frac{S\left(
\sum_{i=1}^{n} X_{ij} (y_{i} - \hat{y}_{i}^{(-j)}),
\lambda\alpha
\right)}
{ \sum_{i=1}^{n} X_{ij}^{2} +\lambda (1- \alpha) }
$$
\tiny
Friedman {\it et al.}, Regularization Paths for Generalized Linear Models via Coordinate Descent (2010)
}

\only<3>{$$
\hat{\theta}_{j} \gets
\frac{\overset{\textsf{\color{red} threshold}}{S}
\left(
\sum_{i=1}^{n} X_{ij} \overbrace{(y_{i} - y_{i}^{(-j)})}^{\textsf{\color{red} residual w/o the variable } \theta_{j} },
\lambda\alpha
\right)}
{ \sum_{i=1}^{n} X_{ij}^{2} + \underbrace{\lambda (1- \alpha)}_{\textsf{\color{magenta} shrinkage}}}
$$
where $S(z, \tau)$ will set it to zero if $|z| < \tau$.
}

## Sum of Single Effect (SuSiE) regression


 \large


``` r
library(susieR)
Y.train <- .sim$y.q[.train]
susie.out <- susie(X, Y.train)
```

 \normalsize

\vfill

\large
$$\mathbf{y} = \sum_{l=1}^{L} \underbrace{ \sum_{j} \mathbf{x}_{j}
    \overset{\textsf{\color{magenta}probabilistic selection}}{\alpha_{j}^{(l)}}
	\overset{\textsf{\color{teal} single variant effect}}{\beta_{j}^{(l)}} }_{\textsf{layer-by-layer}}  + \boldsymbol{\epsilon}$$

\normalsize
where $\sum_{j=1}^{p} \alpha_{j}^{(l)} = 1$ for each layer $l$.

\vfill

\small

[@Wang2020-qu]


## Ideas behind SuSiE

\large

1. Estimating univariate effects is easy (GWAS)

2. Many variants can show similar effects (LD)

3. Let's weight variants probabilistically

4. Regress out the probabilistically reweighted effects

5. Repeat 1-4.

## SuSiE identifies top causal variants


 \scriptsize



 \normalsize


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-51-1} \end{center}

 \normalsize

## SuSiE provides credible sets


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-52-1} \end{center}

 \normalsize

**Credible sets**: Groups of variants that contain the causal variant with high probability (e.g., 95%)

Each color represents one independent signal (one "layer" in the model)


## Geometric intuition: L2, L1, and SuSiE

\centerline{$\mathbf{y} = \mathbf{x}_{1} \theta_{1} + \mathbf{x}_{2} \theta_{2} + \epsilon$}



:::::: {.columns}
::: {.column width=.32}


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-54-1} \end{center}


}

 \normalsize

:::
::: {.column width=.32}


 \scriptsize

\onslide<2->{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-55-1} \end{center}


}

 \normalsize

:::
::: {.column width=.32}


 \scriptsize

\onslide<3>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-56-1} \end{center}


}

 \normalsize

:::
::::::

\scriptsize

Blue cross (×): OLS solution; Red dot: Constrained solution; Blue contours: RSS levels


## Transcriptome-wide association study (TWAS)

Using eQTL, we can understand GWAS mechanisms... [@Gamazon2015-pc]

\centerline{\includegraphics[width=\textwidth]{img/GWAS_blackbox_TWAS.pdf}}

\vfill

$$\mathbf{m}_{g} \sim X \boldsymbol{\alpha}_{g} + \boldsymbol{\epsilon}_{g} \quad \implies \quad \mathbf{y}_{\textsf{GWAS}} = \sum_{g \in \textsf{causal genes}} \mathbf{m}_{g} \beta_{g} + \boldsymbol{\epsilon}_{y}$$

## Mendelian Randomization


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics[height=.5\textheight]{img/MR_DAG} \end{center}


}

 \normalsize

\vfill

\only<1>{
We will use a genetic variable $G$ to mimic an RCT
}

\only<2>{
If a genetic variable $G$ is a valid instrumental variable...
}


## MR in modern epidemiology studies

:::::: {.columns}
::: {.column width=.5}

\Large
\begin{itemize}
\item G: genotype
\item X: APOE protein
\item Y: cancer
\end{itemize}

$$G \to X \to Y$$


 \scriptsize

\onslide<2>{


\begin{center}\includegraphics[height=.2\textheight]{img/MR_GDS} \end{center}


}

 \normalsize

:::
::: {.column width=.5}


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics[height=.5\textheight]{img/MR_APOE} \end{center}


}

 \normalsize

:::
::::::


## MR in action: *FTO* $\to$ fat mass $\to$ obesity, diabetes

:::::: {.columns}
::: {.column width=.25}


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics[width=\linewidth]{img/MR_FTO_BMI_T2D} \end{center}


}

 \normalsize

:::
::: {.column width=.25}


 \scriptsize

\onslide<2->{


\begin{center}\includegraphics[width=\linewidth]{img/MR_FTO_BMI_T2D_2} \end{center}


}

 \normalsize

:::
::: {.column width=.45}

$${}$$

Using genotype as an instrumental variable, we test causality between
other exposure variables and downstream phenotypes.

$${}$$

\begin{itemize}[<+-| alert@+>]
\item Genotype in \textit{FTO} locus $\to$ T2D
\item \textit{FTO} locus $\to$ fat mass
\item Using \textit{FTO} as "instrumental variable", we can ask other MR questions
\end{itemize}

:::
::::::

\tiny

[@Frayling2007-gc]

## Why use Mendelian Randomization?

:::::: {.columns}
::: {.column width=.6}


 \scriptsize


\begin{center}\includegraphics[width=\linewidth]{img/MR_DAG} \end{center}

 \normalsize

:::
::: {.column width=.4}

\Large

1. We do not have enough knowledge about $U$

2. We do not have a way to intervene on $X$, i.e., $do(X=x)$

:::
::::::

\vfill

> "Genotypes are beautifully randomized" [@Fisher1952-fn]

## How can we learn causality between genetic disorders?

Three disorders (just for demonstration)

* OGD: Obsessive `ggplot` disorder

* PTRS: Post-traumatic `R` session stress syndrome

* CIS: Chronic indentation syndrome (as a result of `Python` coding)

## How can we learn causality between genetic traits?

In previous longitudinal studies, we found the following relationships.


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-63-1} \end{center}

 \normalsize

## Can we recover such relationships from GWAS data?


 \scriptsize

\only<1>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-64-1} \end{center}


}

 \normalsize


 \scriptsize

\only<2>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-65-1} \end{center}


}

 \normalsize


 \scriptsize

\only<3>{


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-66-1} \end{center}


}

 \normalsize

## MR measures mediation effects of $X$ to $Y$

<!-- \textsf{\color{teal} Given that $G$ is a valid instrumental variable...} -->

:::::: {.columns}
::: {.column width=.48}

\small

*Example*: 

- $X$: gene expression

- $Y$: disease phenotype

Suppose we have estimated

\onslide<1->{
\begin{eqnarray*}
G &\overset{\alpha}{\to}& X \\
X &=& G \hat{\color{teal}\boldsymbol{\alpha}} + \epsilon_{X}
\end{eqnarray*}
}
\onslide<2->{
and
\begin{eqnarray*}
G &\overset{\gamma}{\to}& Y \\
Y &=& G \hat{\color{magenta}\boldsymbol{\gamma}} + \epsilon_{Y}
\end{eqnarray*}
}

:::
::: {.column width=.48}

\onslide<3->{
\textbf{Goal}: What is the causal effect $\beta$ in $X \overset{\beta}{\to} Y$?
}

\begin{eqnarray*}
\onslide<3->{
Y &=& X \beta + \epsilon' \\
}
\onslide<4->{
G \hat{\color{magenta}\boldsymbol{\gamma}} + \epsilon_{Y} 
&=& 
(G \hat{\color{teal}\boldsymbol{\alpha}} + \epsilon_{X}) \beta + \epsilon' \\
}
\onslide<5->{
G \hat{\color{magenta}\boldsymbol{\gamma}}
&=& 
G \hat{\color{teal}\boldsymbol{\alpha}} {\color{red} \boldsymbol{\beta}} + \cdots
}
\end{eqnarray*}

\onslide<6>{
The answer is as simple as
\begin{eqnarray*}
\mathbb{E}\!\left[{\color{red} \beta}\right] &=& 
\frac{{\color{magenta}\boldsymbol{\gamma}}}
{{\color{teal}\boldsymbol{\alpha}}}
\end{eqnarray*}
}
:::
::::::

## Will any genetic variable work for MR analysis?

Goal: Is PTRS $\to$ OGD or PTRS $\to$ CIS?


 \scriptsize



 \normalsize


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-68-1} \end{center}

 \normalsize

## Select causal variants of exposure (PTRS)?

Goal: Is PTRS $\to$ OGD or PTRS $\to$ CIS?


 \scriptsize



 \normalsize


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-70-1} \end{center}

 \normalsize

## Definition of instrumental variable in MR study

:::::: {.columns}
::: {.column width=.55}


 \scriptsize

\onslide<1->{


\begin{center}\includegraphics[width=.8\linewidth]{img/MR_DAG_IV} \end{center}


}

 \normalsize

:::
::: {.column width=.45}

\begin{itemize}
\item<1-> IV1: The genetic variant $G$ is independent of the potential confounder variable $U$
\item<2-> IV2: The genetic variant is associated with the exposure $X$
\item<3-> IV3: The genetic variant is independent of the outcome $Y$ conditioning on $X$
\end{itemize}

:::
::::::

\vfill

[@Bowden2015-ex]

# Co-localization in genomics

## What is co-localization?

\Large

**Co-localization**: Two traits share the same causal variant(s) in a genomic region

\vfill

* If two GWAS signals co-localize, they likely share a causal variant

* Important for: MR validation, functional follow-up, understanding pleiotropy

\vfill

**Our simulation**:

* PTRS and OGD share 2 causal variants (co-localized)
* PTRS and CIS have independent causal variants (not co-localized)


## Fine-mapping with SuSiE




\large

* SuSiE (Sum of Single Effects) identifies causal variants

* Outputs PIP (Posterior Inclusion Probability) for each variant

* High PIP = likely to be causal

## Co-localization: PTRS and OGD (mediated)


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-73-1} \end{center}

 \normalsize

## No co-localization: PTRS and CIS (independent)


 \scriptsize


\begin{center}\includegraphics{../Fig/regression-classification/unnamed-chunk-74-1} \end{center}

 \normalsize

## `coloc.susie` compare two SuSiE results

\large

**Posterior probabilities (PP)** [@Wallace2021-xd]:

* **PP.H0**: No association in either trait
* **PP.H1**: Association with trait 1 only
* **PP.H2**: Association with trait 2 only
* **PP.H3**: Both associated, different causal variants
* **PP.H4**: Both associated, **shared causal variant** (co-localization)

\vfill

**Rule of thumb**: PP.H4 > 0.8 suggests strong evidence for co-localization

\vfill

\small
coloc.susie allows for multiple causal variants per trait


## Statistical test: coloc with SuSiE

Using coloc.susie [@Wallace2021-xd] to test for co-localization:


 \scriptsize


``` r
## PTRS vs OGD (should co-localize)
coloc.ptrs.ogd$summary
```

```
##    nsnps     hit1    hit2    PP.H0.abf    PP.H1.abf    PP.H2.abf PP.H3.abf
##    <int>   <char>  <char>        <num>        <num>        <num>     <num>
## 1:   841 10907286 9078190 1.688239e-21 0.0211344301 7.817507e-20 0.9786441
## 2:   841 49705291 9078190 2.045804e-06 0.0211225288 9.473239e-05 0.9780921
## 3:   841  5273571 9078190 1.436645e-04 0.0209848186 6.652486e-03 0.9717157
## 4:   841  4115684 9078190 7.781789e-08 0.0211264781 3.603412e-06 0.9782751
## 5:   841  9079464 9078190 9.825459e-22 0.0006701868 4.549747e-20 0.0290930
##       PP.H4.abf  idx1  idx2
##           <num> <int> <int>
## 1: 0.0002214617     1     1
## 2: 0.0006886163     3     1
## 3: 0.0005033485     4     1
## 4: 0.0005947005     5     1
## 5: 0.9702368166     2     1
```

 \normalsize

## Statistical test: coloc with SuSiE

Using coloc.susie [@Wallace2021-xd] to test for co-localization:


 \scriptsize


``` r
## PTRS vs CIS (should NOT co-localize)
coloc.ptrs.cis$summary
```

```
## NULL
```

 \normalsize

## Summary: Visual + Statistical evidence

\Large

**Co-localized (PTRS-OGD)**:

* High PIPs at same positions (mirror image)
* High PP.H4 from coloc.susie
* Suggests shared causal variant

\vfill

**Not co-localized (PTRS-CIS)**:

* High PIPs at different positions
* Low PP.H4, higher PP.H3
* Independent causal variants

## Working with GWAS summary statistics

Often we only have access to GWAS summary statistics, not individual-level data:

* $\hat{\beta}_{j}$: effect size for variant $j$

* $\textsf{SE}(\hat{\beta}_{j})$: standard error

* $Z_{j} = \hat{\beta}_{j} / \textsf{SE}(\hat{\beta}_{j})$: z-score

* $\chi^2_{j} = Z_{j}^2$: chi-square statistic

\vfill

**Question**: Can we infer biological mechanisms from summary statistics alone?

## From univariate to multivariate effects

Underlying assumption: traits are polygenic

\Large

$$Y_{i} = \sum_{j=1}^{p} X_{ij} \theta_{j} + \epsilon$$

\normalsize

**Problem**: Univariate GWAS gives us $\hat{\beta}_j$ for each variant $j$ separately

**Goal**: Understand the true multivariate effects $\theta_j$ from summary statistics

**Challenge**: LD structure confounds the relationship between $\hat{\beta}_j$ and $\theta_j$

## LD-score: Relating $\chi^2$ statistics to heritability

What is a generative model for $\chi_{j}^{2}$ statistics?

\only<1>{
We can show that the univariate z-score relates to multivariate effects:

$$Z_{j} = \frac{\sqrt{n}}{\sigma} \sum_{k} R_{jk} \theta_{k}  + \epsilon_{j}$$

where $R$ is the LD matrix and $\epsilon \sim \mathcal{N}\!\left(0,1\right)$.
}

\onslide<2->{
$$\mathbb{E}\!\left[\chi_{j}^{2}\right] = \mathbb{E}\!\left[Z_{j}^{2}\right] = \mathbb{E} \left( \sqrt{n} \sum_{k} R_{jk} \theta_{k} + \epsilon_{j} \right)^{2}$$
}

\only<3>{
If the effects are independent, i.e., $\mathbb{E}\!\left[\theta_{k}\theta_{j}\right] = 0$ for all $k \neq j$,

$$\mathbb{E}\!\left[\chi^{2}_{j}\right] = n \underbrace{\sum_{k} R_{jk}^{2}}_{\textsf{\color{blue} LD-score}} \mathbb{E}\!\left[ \theta_{k}^{2} \right] + 1$$
}

\vfill
\tiny

Bulik-Sullivan \emph{et al.}, \emph{Nature Genetics} (2014)

## Baseline LD-score to measure polygenic heritability

:::::: {.columns}
::: {.column width=.48}

Assuming that all the variants equally contribute,

$$\mathbb{E}\!\left[\theta_{k}^{2}\right] = \tau / p,$$

where $p$ is the total number of SNPs.

:::
::: {.column width=.48}

Defining an LD score for a variant/SNP $j$ as
$$l_{j} \overset{\textsf{def}}{=} \sum_{k} R^{2}_{jk},$$

:::
::::::

\vfill

We get

$$\mathbb{E}\!\left[\chi^{2}_{j}\right] = \underset{\textsf{\color{teal} sample size}}{n} \underset{\textsf{\color{blue} LD score}}{l_{j}} \underset{\textsf{\color{magenta}per SNP heritability}}{\frac{\tau}{p}} + 1$$

**Regression model**: Regress observed $\chi^2_j$ on LD scores $l_j$ to estimate $\tau$ (heritability)

\vfill
\tiny
Bulik-Sullivan *et al.*, \emph{Nature Genetics} (2014)

## Stratified LD-score regression (S-LDSC)

**Question**: Which genomic annotations explain GWAS heritability?

\vfill

\only<1>{
\centerline{\includegraphics[width=.85\textwidth]{img/GWAS_herit_partition.pdf}}
}

\only<2>{
Partition heritability across $C$ genomic annotations:

\Large
$$\mathbb{E}\!\left[\chi^{2}_{j}\right] = n \sum_{c=1}^{C} \ell_{j,c} \tau_{c} + 1$$

\normalsize

where:

- $\ell_{j,c}$: LD score for variant $j$ with respect to annotation $c$
- $\tau_{c}$: per-SNP heritability in annotation $c$
}

\vfill
\scriptsize
 [@Finucane2015-qr]

## S-LDSC $\to$ tissue/cell type contexts of GWAS

\centerline{\includegraphics[height=.75\textheight]{img/sLDSC.pdf}}

\vfill
\scriptsize

[@Finucane2015-bn]

## Take-home messages

\Large

1. Methods like SuSiE, coloc, and TWAS integrate multiple data types

2. LD score regression connects GWAS signals to heritability

3. Regression is fundamental to dissecting regulatory genomic mechanisms


## Reference {.allowframebreaks}
