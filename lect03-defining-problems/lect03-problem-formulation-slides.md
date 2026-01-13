---
title: "How can we formulate our research questions effectively?"
author: "Yongjin Park"
date: "13 January 2026"
classoption: "aspectratio=169"
fontsize: 12pt
output:
    beamer_presentation:
        theme: Madrid
        keep_md: true
        keep_tex: true
        latex_engine: xelatex
        slide_level: 2
header-includes:
  - \usepackage{cancel}
  - \usepackage{tikz}
  - \usetikzlibrary{positioning,calc,arrows.meta}
  - \tikzset{>={Stealth[length=3mm,width=2mm]}}
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




## Course Logistics

:::::: {.columns}
::: {.column width=.2}

:::
::: {.column width=.7}

\Large

**Please fill out the survey**

of your background

:::
::::::


## Advertisement

\Large

* [ASDA](https://asda.stat.ubc.ca/)

\normalsize

* https://asda.stat.ubc.ca/2025/11/10/free-stats-help.html

* If you and your team are interested in the following ASDA projects, let me know. 

* I will send you one-page project proposal document. 

## ASDA \#1: Blood RNAs for Febrile Illness Diagnosis

* Scoping review identified 286 RNA signatures from literature

* Signatures diagnose or predict outcomes of febrile illnesses (sepsis, TB, infections)

* Goal: Characterize and understand potential sources of biases

* Meta-analysis across many different studies

## ASDA \#2: `Breath print` data for Lung Cancer Detection

* Breath samples from 500+ participants (lung cancer vs. healthy controls)

* Volatile organic compounds (VOCs) analyzed via mass spectrometry

* High-dimensional data: ~500-1500 compounds per sample

* Goal: Build analysis pipeline for biomarker discovery and classification

## Roadmap of STAT/BIOF/GSAT540 2026

_data generation_ followed by _inference_ or validation

\vfill

\large

:::::: {.columns}
::: {.column width=.45}

2. _Statistics and matrix algebra_

3. \textbf{\color{magenta}Formulating biological questions} (Today)

4. \textbf{\color{teal} Hypothesis testing} (Thursday)

5. RNA-seq data analysis

6. Statistical inference

:::
::: {.column width=.45}

8. Single-cell genomics

9. Matrix factorization

10. Genome-wide association studies

11. Multiple hypothesis testing

:::
::::::

\normalsize

##

\huge

A biological problem $\implies$ 

\flushright 

a computational problem

\normalsize

## Guide 1: Crystallize abstract ideas into a concrete computational problem by writing a paper

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Research Model \#1

1. Idea
2. Do research
3. Write a paper

:::

:::
::: {.column width=.45}



::: {.block}
### Research Model \#2

* Idea
* **Write a paper**
* Do research

:::

:::
::::::

Which model would work better for you? Why?

\only<2>{
${}$
\centerline{\textbf{\color{teal} Discussion with your peers around you (5 min).}}
${}$
}

\vfill

\normalsize

\begin{enumerate}
\item<3-> Model \#2 forces us to be clear, focused
\item<3-> Crystallizes what we don’t understand yet (and what we don't need to know)
\item<3-> Opens the way to dialogue with others: reality check and collaboration
\end{enumerate}

\tiny

Simon Peyton Jones, [How to write a great research paper](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/07/How-to-write-a-great-research-paper.pdf)

\normalsize

##

\huge
\begin{center}
${}$ Model \#1 or \#2?
\end{center}

## Guide 2: Ten simple rules ...

:::::: {.columns}
::: {.column width=.3}

:::
::: {.column width=.66}

\large

William Noble,

[**Ten simple rules for defining a computational biology project.**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010786)

*PLoS Comput Biol* (2023)

\normalsize

:::
::::::


## Summary: Ten Simple Rules

\large

:::::: {.columns}
::: {.column width=.43}

1. \textbf{\color{teal} Write it down}
2. \textbf{\color{teal} Define the problem}
3. Define potential impact
4. Summarize related work
5. \textbf{\color{teal} State central hypothesis}

:::
::: {.column width=.47}

6. Sketch your approach
7. \textbf{\color{teal} Identify available data}
8. Choose validation measures
9. Set up version control
10. Share with colleagues

:::
::::::

\normalsize

$${}$$

\textbf{\color{blue}\large Discussion: How many items have you followed in the past?}

## Rule 1: Write it down

\Large

Document your idea in "writing"

\normalsize

\vfill

* Clarifies your thinking and the underlying problem

* Exposes half-baked concepts early

* Forces you to articulate assumptions

> "The act of writing your idea down will inevitably clarify it for you."

\vfill

**Even a rough draft is better than no draft**

## Rule 2: Define the problem

\Large

Clearly delineate the scope and nature of the problem

\normalsize

\vfill

* What exactly are you trying to solve?

* What is **in scope** vs. **out of scope**?

* Be specific, not vague


  * ma\vfill

Which sounds more specific?

:::::: {.columns}
::: {.column width=.45}

* Understand gene expression programs in myeloid cancer

:::
::: {.column width=.45}

* Predict enhancer-promoter interaction in K562 cells

:::
::::::


<!-------------------->
<!-- skip this item -->
<!-------------------->

<!------------------------------------------------------>
<!-- ## Rule 3: Define potential impact				  -->
<!-- 												  -->
<!-- \Large											  -->
<!-- 												  -->
<!-- Identify who benefits and how					  -->
<!-- 												  -->
<!-- \normalsize									  -->
<!-- 												  -->
<!-- \vfill											  -->
<!-- 												  -->
<!-- * Who is your audience?						  -->
<!-- * What concrete benefits will your work provide? -->
<!-- * How does it advance the field?				  -->
<!-- 												  -->
<!-- \vfill											  -->
<!-- 												  -->
<!-- This helps:									  -->
<!-- 												  -->
<!-- 1. Motivate stakeholders						  -->
<!-- 2. Write compelling grant applications			  -->
<!-- 3. Frame your paper's contribution				  -->
<!------------------------------------------------------>

## Rule 4: Summarize related work

\Large

Conduct literature review early

\normalsize

\vfill

> "A terrible trap to fall into is carrying out an entire research
> project only to find that the same work was already published before
> you began."

* Learn from others' approaches

* Identify gaps in current knowledge

* **Summarize** methods in a chronologically ordered table

## Rule 5: State central hypothesis

\Large

Reduce your idea to a single sentence

\large

\vfill

* Focus on the one, single, **overarching** theme

* Make it **testable**

* Keep it **simple and clear**

\vfill

\normalsize

* It takes practice and time.

<!---------------------------------------------------------------->
<!-- ## Rule 6: Sketch your approach						    -->
<!-- 														    -->
<!-- \Large													    -->
<!-- 														    -->
<!-- Describe the computational method in plain English		    -->
<!-- 														    -->
<!-- \normalsize											    -->
<!-- 														    -->
<!-- \vfill													    -->
<!-- 														    -->
<!-- * What is your strategy?								    -->
<!-- * Highlight strengths and weaknesses					    -->
<!-- * Compare to alternative approaches					    -->
<!-- * Don't dive into code yet								    -->
<!-- 														    -->
<!-- \vfill													    -->
<!-- 														    -->
<!-- **Goal:** Communicate the *idea* before the implementation -->
<!---------------------------------------------------------------->

## Rule 7: Identify available data

\Large

Ensure you have both inputs and validation data

\normalsize

\vfill

* What are the input?

* How do we know that our discovery is correct?

* Simulation studies are okay.

<!------------------------------------------------------------->
<!-- ## Rule 8: Choose validation measures					 -->
<!-- 														 -->
<!-- \Large													 -->
<!-- 														 -->
<!-- Design quantitative metrics aligned with end-user needs -->
<!-- 														 -->
<!-- \normalsize											 -->
<!-- 														 -->
<!-- \vfill													 -->
<!-- 														 -->
<!-- * Don't rely only on textbook metrics					 -->
<!-- * Consider what matters to biologists					 -->
<!-- * Define success criteria upfront						 -->
<!-- 														 -->
<!-- \vfill													 -->
<!-- 														 -->
<!-- Examples:												 -->
<!-- 														 -->
<!-- * Classification: AUC, precision, recall				 -->
<!-- * Prediction: correlation, RMSE						 -->
<!-- * Discovery: enrichment, biological plausibility		 -->
<!------------------------------------------------------------->

## Rule 9: Set up version control

\Large

Use `git` or a similar tool from day one

\normalsize

\vfill

* Track code evolution from project inception

    * Keep the dates

* Enable collaboration

* Maintain reproducibility

* Document your progress


## Rule 10: Share with colleagues

\Large

Seek feedback on feasibility and scientific importance

\normalsize

\vfill

* Prepare to "Think Again" (don't be embarrassed or offended by being wrong)

* Refine your ideas through discussion

* Find potential collaborators

## Ten Simple Rules (items) $\to$ A paper

:::::: {.columns}
::: {.column width=.43}

1. \textbf{\color{teal} Write it down}
2. \textbf{\color{teal} Define the problem}
3. Define potential impact
4. Summarize related work
5. \textbf{\color{teal} State central hypothesis}
6. Sketch your approach
7. \textbf{\color{teal} Identify available data}
8. Choose validation measures
9. Set up version control
10. Share with colleagues

:::
::: {.column width=.47}

* Abstract
  * Rule \#1, \#2, \#5

* Introduction
  * Rule \#2, \#3, \#4, \#5

* Results
  * Rule \#5, \#7, \#10

* Discussion
  * Rule \#10

* Methods
  * Rule \#6, \#7, \#8, \#9, \#10

:::
::::::


## Bottom line

\Large

${}$ \textbf{It's so critical to describe your problem clearly.}

\begin{itemize}
\item <2-> Recognize the underlying computational problem(s)
\item <3-> Make an actionable plan to solve them.
\item <4> Evaluate what/how it worked and didn't.
\end{itemize}

## Is defining a problem worth my time?

\centerline{\includegraphics[height=.72\textheight]{img/automation_2x.png}}

\vfill

\tiny
[https://xkcd.com/1319/](https://xkcd.com/1319/)
\normalsize

## Two flavors

:::::: {.columns}
::: {.column width=.45}

* \#1: Optimization approach

	* \textbf{\color{magenta} The art of abstraction}

	* Define **input**, **output**, and **metric** (objective function)

	* Find the same or similar existing problems in the literature

:::
::: {.column width=.45}

* \#2: Probabilistic Modelling

	* \textbf{ Think of a potential {\color{teal} generative model}}

	* Describe the problem in terms of random variables

	* Structural Equation Model (in epidemiology)

	* Use an existing inference algorithm

:::
::::::

# Define a problem as an optimization problem

## The art of abstraction: Same data, different questions

\Large

Single-cell RNA-seq data: a matrix of **genes $\times$ cells**

\normalsize

\vfill

Biological questions $\rightarrow$ computational problems


## Exercise: Eleven Grand Challenges in Single-Cell

\footnotesize

Lähnemann et al., *Genome Biology* (2020)

\normalsize

:::::: {.columns}
::: {.column width=.48}

**Transcriptomics**

1. Handling sparsity (dropout)
2. \textbf{\color{magenta} Defining cell identity}
3. Tracking cellular dynamics
4. Mapping to reference atlases

**Genomics**

5. Finding variants (scDNA-seq)
6. Identifying copy-number changes
7. Inferring regulatory networks

:::
::: {.column width=.48}

**Phylogenetics**

8. Understanding tumor evolution
9. Calling structural variants

**Integration**

10. Multi-modal integration
11. Scaling to millions of cells

:::
::::::

Let's work on one of these examples

## Example: Defining cell identity (Challenge 2)

\Large

**Biological question:** What cell types are in my sample?

\vfill


* What are the variables of interest?

* Is it important?

\normalsize


## How should we think of cells in this context?

\Large

Data: $m\times{n}$ expression matrix $\mathbf{\color{magenta} X}$ ($m$ genes $\times$ $n$ cells)

\vfill

\only<1>{

$$
\mathbf{X} = \begin{pmatrix}
x_{11} & x_{12} & \cdots & x_{1n} \\
x_{21} & x_{22} & \cdots & x_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
x_{m1} & x_{m2} & \cdots & x_{mn}
\end{pmatrix}
=
\left(
\mathbf{x}_{1}, \mathbf{x}_{2}, \ldots, \mathbf{x}_{n}
\right)$$
where
$m$ corresponds to 19k genes.
}

\vfill

\begin{itemize}
\item<2> Is a cell a vector of gene expressions over 19k genes?
\item<2> Is a cell a vector of gene expressions over several genes?
\item<2> Is a cell a point in a "cell type landscape?"
\item<2> ...
\end{itemize}

\normalsize

## Let's start writing out the problem

\Large

**Biological question:** What cell types are in my sample?

\normalsize

\vfill

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Input
Gene expression $X$ ($m$ genes $\times$ $n$ cells)
:::

:::
::: {.column width=.45}

::: {.block}
### Output
Assign each cell vector $\mathbf{x}_{j}$ to label $c_i \in \{1, \ldots, K\}$
:::

:::
::::::

## Brainstorming ...

* How did the cell types emerge in multicellular organisms?

* How did the cell express genes?

* What do we measure with a cell-level gene expression matrix?

* How do we know cells belong to a cell type A, not B?


## Example: Let's write it down with several sentences

1. We can consider each cell as a point along a differentiation trajectory.

2. We can represent the relationships between cells as a graph (connecting points).

3. Cells that share similar expression profiles are closely located in a graph space.

\vfill

Can we do better? Is it specific enough?


## How do we measure the distance between cells?

Multiple levels of distance measures between two cell vectors:

1. **Use all genes**: Compute distance using all $m$ features (e.g., Euclidean, cosine)

2. **Focus on selected genes**: Use only marker genes or highly variable genes

3. **Learn new coordinates**: Project cells into a lower-dimensional space (e.g., PCA) and compute distance there


## A toy example: 100 genes and 200 cells/samples




\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-1-1} \end{center}

## Compute cell-cell distance in Top 5 PCA space



:::::: {.columns}
::: {.column width=.2}

* row = cells; column = cells

${}$

Note: _Don't worry about what PCA is yet (`lect-09`)_

${}$

Instead of dealing with 100 (genes) dimensions, we only need to deal with 5 dimensions.

:::
::: {.column width=.4}

\onslide<1->{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-3-1} \end{center}


}

\tiny
blue (small); red (large)
\normalsize

:::
::: {.column width=.3}

\onslide<2>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-4-1} \end{center}


}

:::
::::::


## Visualize cells in PC space






\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-7-1} \end{center}

## Reminder: Our goal is to build a graph of cells

\large

* Graph = connecting dots (dot = cell)

* Tree = a special case of graph

* Minimum Spanning Tree = one way to connect the dots

* What is our assumption?

\normalsize


## Building Minimum Spanning Tree





\only<1>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-1} \end{center}

}
\only<2>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-2} \end{center}

}
\only<3>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-3} \end{center}

}
\only<4>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-4} \end{center}

}
\only<5>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-5} \end{center}

}
\only<6>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-6} \end{center}

}
\only<7>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-7} \end{center}

}
\only<8>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-8} \end{center}

}
\only<9>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-9} \end{center}

}
\only<10>{


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/mst_animation-10} \end{center}

}

## Build minimum spanning tree (MST)


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-10-1} \end{center}

## Graph-based 2D layout show the MST more explicitly



:::::: {.columns}
::: {.column width=.48}

**PCA layout**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-12-1} \end{center}

:::
::: {.column width=.48}

**Graph-based layout**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-13-1} \end{center}

:::
::::::

## What if we build Maximum Spanning Tree?



:::::: {.columns}
::: {.column width=.48}

**Minimum Spanning Tree**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-15-1} \end{center}

:::
::: {.column width=.48}

**Maximum Spanning Tree**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-16-1} \end{center}

:::
::::::


## Alternative: k-Nearest Neighbor (kNN) Graph



:::::: {.columns}
::: {.column width=.48}

**PCA layout**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-18-1} \end{center}

:::
::: {.column width=.48}

**Force-directed layout**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-19-1} \end{center}

:::
::::::



## Let's refine our paragraph on the problem definition.

1. We consider each cell as a point along a differentiation trajectory. We assume that \textbf{\color{teal} most genes are affected by the differentiation process}, so cell-to-cell distances are best captured in "principal" component analysis.

2. We represent the relationships between cells as a graph, \textbf{\color{magenta} more precisely, a tree structure reflecting their developmental lineage.}

3. Since gene expressions are \textbf{\color{olive} gradually changed} over the differentiation process, each \textbf{\color{olive} branch of the underlying tree structure} represents distinctive cell types.

4. We annotate cell type labels based on prior knowledge (more work needed).

## More formal definition of the problem

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Input
* Gene expression $X$ ($m$ genes $\times$ $n$ cells)
* Assumption: a differentiation process (minimum spanning tree)
:::

:::
::: {.column width=.45}

::: {.block}
### Output
* Cell-to-cell MST graph
* Assign each cell vector $\mathbf{x}_{j}$ to label $c_i \in \{1, \ldots, K\}$
:::

:::
::::::

## More specific problem definition of MST

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Input
Distance graph $G = (V, E, D)$
where

* $V$ = vertices (cells)

* $E$ = edges between cells

* $D$ = distance matrix

:::

:::
::: {.column width=.45}

::: {.block}
### Output
Tree $T$ connecting all cells
:::

:::
::::::

\vfill

:::::: {.columns}
::: {.column width=.45}

::: {.block}
### Objective
Minimize total edge weight (MST)
:::

:::
::: {.column width=.45}

$$
\min_{T} \sum_{(i,j) \in E(T)} D_{ij}
$$

:::
::::::

## Optional: more on identifying branches in the tree

:::::: {.columns}
::: {.column width=.5}


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-20-1} \end{center}

:::
::: {.column width=.4}

* True "state" is hidden (need to estimate)

* We will revisit how we can do clustering on graphs

* But the following three slides are there for the sake of completeness

:::
::::::


## Optional: Spectral clustering on MST



:::::: {.columns}
::: {.column width=.48}

**Eigenvalues of Laplacian**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-22-1} \end{center}

\tiny
Gap at k=3 suggests 3 clusters
\normalsize

:::
::: {.column width=.48}

**Spectral embedding**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-23-1} \end{center}

\tiny
Fiedler vector separates branches
\normalsize

:::
::::::

## Optional: Spectral cut results



:::::: {.columns}
::: {.column width=.48}

**True cell states**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-25-1} \end{center}

:::
::: {.column width=.48}

**Spectral clustering (k=3)**


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-26-1} \end{center}

:::
::::::


# Use probabilities to describe data generation process

## Why Probabilistic (Graphical) Model?

PGM forces us to be explicit about:

\vfill

1. **What we observe** (data)

2. **What we want to learn** (parameters, latent variables)

3. **What assumptions we make** (model structure)

4. **What we don't know** (uncertainty)

\vfill

> "All models are wrong, but some are useful." -- George Box

## Describe your science by stitching together random variables

\Large

Represent joint probability distributions using graphs

\normalsize

\vfill

:::::: {.columns}
::: {.column width=.45}

**Nodes** = Random variables

* Observed (shaded)
* Latent/hidden (open)
* Parameters (dots or small nodes)

:::
::: {.column width=.45}

**Edges** = Dependencies

* Arrow $\rightarrow$ = causal/generative direction
* No arrow = correlation/association

:::
::::::

## `lect-02`: What is statistics? 


\centerline{\includegraphics[width=.55\textwidth]{img/probability-inference.pdf}}

\tiny 

Wasserman, _All of Statistics_ (2003)

\vfill

\normalsize

As long as you can write out the probability...

## 


## Central question: How did we get our data?

\Large

1. **Write down your data generating model**

2. Inference: Fit observed data to the model

3. Validation: How well? How significant?

## Example: the number of mutations in HIV (`lect02`)

:::::: {.columns}
::: {.column width=.55}

Let's graphically represent

$$X \mid \lambda \sim \text{Poisson}(\lambda)$$

where

* $\lambda$: unknown mutation rate parameter (open)

* $X$: observed number of mutations (shaded)

:::
::: {.column width=.35}

\begin{center}
\begin{tikzpicture}[
    node distance=1.5cm,
    observed/.style={circle, draw, fill=gray!30, minimum size=1cm},
    latent/.style={circle, draw, minimum size=1cm},
    param/.style={circle, draw, fill=black, minimum size=0.2cm}
]
\node[latent] (lambda) {$\lambda$};
\node[observed, below=of lambda] (x) {$X$};

\draw[->, thick] (lambda) -- (x);
\end{tikzpicture}
\end{center}

:::
::::::


## Revisit the same question: cell type identification

**Data**: RNA-seq of **genes $\times$ cells**

\large

Let $X_{gi}$ be the expression value of a gene $g$ in a cell $i$.

$$
\mathbf{X} = \begin{pmatrix}
X_{11} & X_{12} & \cdots & X_{1n} \\
X_{21} & X_{22} & \cdots & X_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
X_{m1} & X_{m2} & \cdots & X_{mn}
\end{pmatrix}
=
\left(
\mathbf{x}_{1}, \mathbf{x}_{2}, \ldots, \mathbf{x}_{n}
\right)$$
where
$m$ corresponds to 19k genes.


## The number of mRNA counts in the same way?

:::::: {.columns}
::: {.column width=.55}

Let's graphically represent

$$X_{gi} \mid \lambda \sim \text{Poisson}(\lambda_{gi})$$

where

* $\lambda_{gi}$: unknown transcription rate parameter of a gene $g$ within a cell $i$ (open; latent; parameter)

* $X_{gi}$: the number of mRNA reads of a gene $g$ in a cell $i$ (shaded; observed; data)

:::
::: {.column width=.35}

\begin{center}
\begin{tikzpicture}[
    node distance=1.5cm,
    observed/.style={circle, draw, fill=gray!30, minimum size=1cm},
    latent/.style={circle, draw, minimum size=1cm},
    param/.style={circle, draw, fill=black, minimum size=0.2cm}
]
\node[latent] (lambda) {$\lambda_{gi}$};
\node[observed, below=of lambda] (x) {$X_{gi}$};

\draw[->, thick] (lambda) -- (x);
\end{tikzpicture}
\end{center}

:::
::::::

## Can we elaborate more on the $\lambda_{gi}$ (gene $\times$ cell)?




\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-28-1} \end{center}

\vfill

What's wrong with having $\lambda_{gi}$ for each $X_{gi}$?

## Let's think of the "generation" of $X$ (or $\lambda$)

HINT (just sorted rows and columns):




\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-30-1} \end{center}

## Write out data generating scheme

:::::: {.columns}
::: {.column width=.45}

\large

\onslide<1->{Of the $K$ cell types:

$$Z_{ik} = 1 \iff \text{cell $i$ in type $k$}$$

}

\onslide<2->{
For each cell type $k$, a gene $g$ is expressed with the rate

$$\lambda_{gk}$$

*Note*: do we have more $\lambda_{gi}$'s or $\lambda_{gk}$'s?
}

\normalsize

:::
::: {.column width=.45}

\large

\onslide<3->{
Given the cell type assignment $Z_{ik}$ and type specific expression rate parameters $\lambda$, we have

$$X_{gi} | Z_{ik}=1 \sim \mathsf{Poisson}(\lambda_{gk})$$
}

\normalsize

:::
::::::

## Or we could simply draw it out

\large

\begin{center}
\begin{tikzpicture}[
    node distance=1.2cm,
    observed/.style={circle, draw, fill=gray!30, minimum size=0.9cm},
    latent/.style={circle, draw, minimum size=0.9cm}
]
\node[latent] (z) {$Z_{ik}$};
\node[latent, right=1.5cm of z] (lambda) {$\lambda_{gk}$};
\node[observed, below=1.2cm of $(z)!0.5!(lambda)$] (x) {$X_{gi}$};

\draw[->, thick] (z) -- (x);
\draw[->, thick] (lambda) -- (x);
\end{tikzpicture}
\end{center}

What's next?

## Where do these $Z$ and $\lambda$ variables correspond to?



:::::: {.columns}
::: {.column width=.1}

:::
::: {.column width=.85}

\onslide<2>{\large $Z_{ik} \to \to \to$}

:::
::::::

:::::: {.columns}
::: {.column width=.1}

\onslide<2>{\large
\begin{align*}
\lambda_{gk} \\
\downarrow \\
\downarrow
\end{align*}
}

:::
::: {.column width=.85}


\begin{center}\includegraphics{../Fig/lect03-problem-formulation/unnamed-chunk-31-1} \end{center}

:::
::::::



## In many cases we have assignment $Z$ observed/given

\large

\begin{center}
\begin{tikzpicture}[
    node distance=1.2cm,
    observed/.style={circle, draw, fill=gray!30, minimum size=0.9cm},
    latent/.style={circle, draw, minimum size=0.9cm}
]
\node[observed] (z) {$Z_{ik}$};
\node[latent, right=1.5cm of z] (lambda) {$\lambda_{gk}$};
\node[observed, below=1.2cm of $(z)!0.5!(lambda)$] (x) {$X_{gi}$};

\draw[->, thick] (z) -- (x);
\draw[->, thick] (lambda) -- (x);
\end{tikzpicture}
\end{center}

What's next (`lect04`)? 

$$\lambda_{g1} = \lambda_{g2} \quad 
\text{or} \quad \lambda_{g1} \neq \lambda_{g2}$$

## As long as we have PGM...

\large

* Estimate the unknown parameter, such as $\lambda$

* Check how "well" generated data, such as $X$, can explain observed $X$

## Summary

\large

**Two approaches to formulate computational problems:**

\normalsize

:::::: {.columns}
::: {.column width=.45}

1. **Optimization approach**
   * Define input, output, metric
   * The art of abstraction
   * Example: MST for cell trajectories

:::
::: {.column width=.45}

2. **Probabilistic modeling**
   * Write down generative model
   * Use graphical models (PGM)
   * Example: Poisson model for gene expression

:::
::::::

\vfill

\large

**Key takeaway:** Write down your problem clearly before coding!


