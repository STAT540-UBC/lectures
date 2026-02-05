options(stringsAsFactors = FALSE)

`%&%` <- function(a,b) paste0(a, b)
`%r%` <- function(mat,rr) mat[rr, , drop = FALSE]
`%c%` <- function(mat,cc) mat[, cc, drop = FALSE]

library(tidyverse)
library(data.table)
library(ggrepel)

num.int <- function(...) format(..., justify="none", big.mark=",", drop0trailing = TRUE)

num.sci <- function(...) format(..., justify="none", digits=2, scientific = TRUE)

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

col.order <- function(pair.tab, .ro, ret.tab = FALSE) {

    M = pair.tab %>%
        dplyr::select(row, col, weight) %>%
        dplyr::mutate(row = factor(row, .ro)) %>%
        tidyr::pivot_wider(names_from = col, values_from = weight, values_fill = 0, values_fn = mean)

    co = order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co = colnames(M)[-1][co]
    if(ret.tab) {
        ret = pair.tab %>%
            dplyr::mutate(row = factor(row, .ro)) %>%
            dplyr::mutate(col = factor(col, .co))
    } else {
        ret = .co
    }
    return(ret)
}

order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)
    require(data.table)

    .tab <- pair.tab %>% dplyr::select(row, col, weight) %>% as.data.table

    M <- dcast(.tab, row ~ col, value.var = "weight", fill = 0, fun.aggregate = mean)
    rr <- M[[1]]  # Get first column (row names)
    cc <- colnames(M)[-1]

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    # Convert to matrix
    M.mat <- M[, -1] %>% as.matrix
    ro <- row.order(M.mat)

    ## log.msg('Sort the rows: %d', length(ro))
    co <- order(apply(M.mat[ro, , drop=FALSE], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret <- pair.tab %>%
            dplyr::mutate(row = factor(row, rr[ro])) %>%
            dplyr::mutate(col = factor(col, cc[co]))
    } else {
        ret <- list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}


.sort.matrix <- function(.X) {
    as.matrix(.X) %>%
        reshape2::melt() %>%
        rename(row = Var1, col = Var2, weight = value) %>%
        order.pair(ret.tab=TRUE) %>%
        as.data.table %>%
        dcast(row ~ col, value.var = "weight") %>%
        dplyr::select(-row) %>%
        as.matrix
}

.rnorm <- function(d1, d2) {
    matrix(rnorm(d1 * d2), d1, d2)
}

###############################################################
.matshow <- function(.mat, .lab=0, .axis.lab=0, .lw=0, .scale=TRUE) {

    library(ggrastr)

    .mat <- as.matrix(.mat)
    .cols <- colnames(.mat)
    .rows <- rownames(.mat)
    if(length(.cols) < ncol(.mat)){
        colnames(.mat) <- str_c("c", 1:ncol(.mat))
    }
    if(length(.rows) < nrow(.mat)){
        rownames(.mat) <- str_c("r", 1:nrow(.mat))
    }
    .cols <- colnames(.mat)
    .rows <- rownames(.mat)

    # Convert to data frame and ensure all columns are named
    .df <- as.data.frame(.mat, check.names = FALSE)
    # Fix any empty column names
    empty_names <- which(colnames(.df) == "" | is.na(colnames(.df)))
    if(length(empty_names) > 0) {
        colnames(.df)[empty_names] <- paste0("V", empty_names)
    }

    .dt <-
        .df %>%
        tibble::rownames_to_column("row") %>%
        as.data.table %>%
        melt(id.vars = "row", variable.name = "col") %>%
        dplyr::mutate(row = factor(as.character(row), rev(.rows))) %>%
        dplyr::mutate(col = factor(as.character(col), .cols))

    ret <-
        ggplot(.dt, aes(y = row, x = col, fill = value)) +
        theme(legend.position = "none")
    if(.axis.lab > 0){
        ret <-
            ret +
            theme(axis.text.x = element_text(size=.axis.lab, angle=90, vjust=1, hjust=1)) +
            theme(axis.text.y = element_text(size=.axis.lab))
    } else {
        ret <-
            ret +
            theme(axis.text = element_blank()) +
            theme(axis.ticks = element_blank())
    }
    ret <- ret +
        theme(axis.title = element_blank())

    if(.lw > 0){
        ret <- ret +
            ggrastr::rasterise(geom_tile(linewidth = .lw, colour = "black"), dpi=300)
    } else {
        ret <- ret + ggrastr::rasterise(geom_tile(), dpi=300)
    }

    if(.scale){
        ret <- ret +
            scale_fill_distiller(palette = "RdBu", direction = -1)
    } else {
        ret <- ret +
            scale_fill_distiller(palette = "Greys", direction = 1)
    }

    if(.lab > 0) {
        ret <- ret +
            geom_text(aes(label = round(value,1)), size = .lab)
    }

    return(ret)
}

################################################################
.gg.plot <- function(...) {
    ggplot(...) +
        theme_classic() +
        theme(axis.title = element_text(size=8)) +
        theme(axis.text = element_text(size=6)) +
        theme(legend.spacing = unit(.1, "lines"),
              legend.key.size = unit(.5, "lines"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=6),
              panel.background = element_rect(fill='transparent'),
              plot.background = element_rect(fill='transparent', color=NA),
              legend.background = element_rect(fill='transparent', size=0.05),
              legend.box.background = element_rect(fill='transparent'))
}

as.dt <- function(x, col.names=NULL) {
    .mat <- as.matrix(x)
    if(is.null(col.names)) col.names <- str_c("V",1:ncol(.mat))
    colnames(.mat) <- col.names
    as.data.table(.mat)
}

################################################################
if.needed <- function(.file, .code) {
    if(!all(file.exists(unlist(.file)))){ .code }
    stopifnot(all(file.exists(unlist(.file))))
}

################################################################
setup.env <- function(fig.dir) {

    ## save figures here
    dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)
    knitr::opts_chunk$set(warning = FALSE, message = FALSE, fig.path = fig.dir)
    knitr::opts_chunk$set(echo=FALSE, fig.align="center")

    ## allow the code to chunk set size="tiny" ##
    hook.chunk  <- knitr::knit_hooks$get("chunk")

    ## Default: normalsize -> scriptsize
    ## This will redefine normalsize
    knitr::knit_hooks$set(chunk = function(x, options) {
        x <- hook.chunk(x, options)
        ifelse(options$size != "normalsize",
               paste0("\n \\", options$size, "\n\n", x, "\n\n \\normalsize"),
               paste0("\n \\", "scriptsize", "\n\n", x, "\n\n \\normalsize"))
    })

    ## show plot one by one in beamer ##
    hook.plot <- knitr::knit_hooks$get("plot")

    knitr::knit_hooks$set(plot = function(x, options) {
        if (!is.null(options$onslide.plot)) {
            bf <- paste0("\\onslide<", options$onslide.plot, ">{")
            ret <- paste(c(bf, knitr::hook_plot_tex(x, options), "}"),
                         collapse = "\n")
            return(ret)
        } else if (!is.null(options$only.plot)) {
            bf <- paste0("\\only<", options$only.plot, ">{")
            ret <- paste(c(bf, knitr::hook_plot_tex(x, options), "}"),
                         collapse = "\n")
            return(ret)
        }
        return(hook.plot(x, options))
    })

    knitr::opts_chunk$set(dev = "cairo_pdf")
}

################################################################
make.gs.lol <- function(.dt) {
    .dt <- as.data.table(.dt) %>% unique()
    .list <-
        .dt[, .(gene = .(gene_symbol)), by = .(gs_name)] %>%
        as.list()
    .names <- .list$gs_name
    .ret <- .list$gene
    names(.ret) <- .names
    return(.ret)
}

read.panglao <- function(.organs = NULL, .file = "../data/PanglaoDB_markers_27_Mar_2020.tsv.gz"){
    if(is.null(.organs)){
        .dt <- fread(.file)
    } else {
        .dt <- fread(.file)[organ %in% .organs]
    }
    ret <- .dt %>%
        rename(gene_symbol = `official gene symbol`, gs_name = `cell type`) %>%
        mutate(gs_name = gsub(pattern="[ ]*cells", replacement = "", gs_name)) %>%
        make.gs.lol()
}

run.beta.gsea <- function(.dt, .db){
    vars <- unique(.dt$Var2)
    ret <- data.table()
    for(v in vars){
        .scores <- as.vector(scale(.dt[Var2 == v]$value))
        names(.scores) <- .dt[Var2 == v]$Var1
        .gsea <- fgsea::fgsea(pathways = .db,
                              stats = .scores,
                              scoreType = "pos") %>%
            as.data.table()
        .gsea[, Var2 := v]
        ret <- rbind(ret, .gsea)
    }
    return(ret)
}

sort.beta <- function(.dt, genes.selected){
    .genes <- rownames(.dt)
    .mat <- as.matrix(.dt)
    colnames(.mat) <- 1:ncol(.mat)
    rownames(.mat) <- .genes
    .mat <- .mat[rownames(.mat) %in% genes.selected, , drop = F]
    .melt <- reshape2::melt(.mat) %>% as.data.table()
    .melt <- .melt[order(`value`, decreasing = T), head(.SD, 1),
                   by = .(Var1, Var2)]
    .mat <- dcast(.melt, Var1 ~ Var2, value.var="value", fun.aggregate = max)
    .rows <- unlist(.mat[,1], use.names = F)
    .mat <- .mat[,-1]
    .ro <- apply(t(.mat), 2, which.max) %>%
        order(decreasing = F) %>%
        (function(x) .rows[x])
    .to <- colnames(.mat)
    .melt[, row := factor(`Var1`, .ro)]
    .melt[, variable := factor(`Var2`, .to)]
    return(.melt)
}

select.gsea.beta <- function(.gsea, .beta, ntop = 3, nmax = 30, sd.cutoff = 1e-4, padj.cutoff = .1){
    .cts <- unlist(.gsea[order(pval),
                         head(.SD, ntop),
                         by = .(Var2)][padj < padj.cutoff,
                                       .(pathway)])
    .gsea.show <- .gsea[pathway %in% .cts]
    .gsea.show[, variable := factor(`Var2`, 0:100)]
    .temp <- dcast(.gsea.show, variable ~ pathway,
                   value.var = "pval",
                   fun.aggregate = min)
    .po <- order(apply(-log10(1e-20 + .temp[,-1]), 2, which.max))
    .po <- colnames(.temp)[-1][.po]
    .gsea.show[, ct := factor(pathway, .po)]
    .sd <- .beta[, .(s = sd(`value`)), by = .(Var1)]
    .genes.show <- intersect(unlist(.gsea.show$leadingEdge),
                             .sd[s > sd.cutoff]$Var1)
    .beta.show <- .beta[Var1 %in% .genes.show]
    .temp <- copy(.gsea.show)
    .temp[, r := as.character(.I)]
    .leading.edges <- .temp[, leadingEdge[[1]], by=r] %>%
        left_join(.temp[, .(r, ct)]) %>%
        mutate(row = factor(V1, sort(unique(.beta.show$row))))
    .leading.genes <- .beta.show[row %in% .leading.edges$row][order(value, decreasing = T), head(.SD, 10), by = .(variable)][order(value, decreasing = T)]$row
    .leading.genes <- head(unique(.leading.genes), nmax)
    list(gsea = .gsea.show, beta = .beta.show, leading.edges = .leading.edges,
         leading.genes = .leading.genes)
}

