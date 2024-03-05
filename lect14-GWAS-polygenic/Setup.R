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
        tidyr::spread(key = col, value = weight, fill = 0)

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

    .tab = pair.tab %>% dplyr::select(row, col, weight)

    M = .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr = M[, 1] %>% unlist(use.names = FALSE)
    cc = colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro = row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## log.msg('Sort the rows: %d', length(ro))
    co = order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret = pair.tab %>%
            dplyr::mutate(row = factor(row, rr[ro])) %>%
            dplyr::mutate(col = factor(col, cc[co]))
    } else {
        ret = list(rows = rr[ro], cols = cc[co], M = M)
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

    .dt <-
        as.data.frame(.mat) %>%
        rownames_to_column("row") %>% 
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
            scale_fill_distiller(palette = "Paired", direction = -1)
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

    knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = TRUE)
    knitr::opts_chunk$set(dev = "cairo_pdf")
}

