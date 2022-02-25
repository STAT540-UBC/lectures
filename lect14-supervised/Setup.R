setup.env <- function(fig.dir) {

    ## save figures here
    dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)
    knitr::opts_chunk$set(warning = FALSE, message = FALSE, fig.path = fig.dir)
    knitr::opts_chunk$set(results = FALSE, fig.align="center")

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
}
