sort.theta <- function(.dt){

    .co <- order(apply(as.matrix(.dt[, -1]), 1, which.max))
    .to <- colnames(.dt)[-1]

    .melt <- .dt %>%
        mutate(col = factor(1:n(), rev(.co))) %>%
        melt(id.vars = c("col", "celltype")) %>%
        mutate(variable = factor(variable, .to))

    .cto <-
        .melt[,
              .(value = mean(value)),
              by = .(celltype, variable)] %>%
        dcast(celltype ~ variable, value.var = "value") %>%
        (function(x){
            .mat <- as.matrix(x[,-1])
            .o <- order(apply(.mat, 1, which.max), decreasing = T)
            rev(x$celltype[.o])
        })

    .melt[, celltype := factor(celltype, .cto)]
    return(.melt)
}

read.theta <- function(.file){
    return(sort.theta(fread(.file)))
}

plot.theta <- function(.file, with.ct = F){
    .dt <- read.theta(.file)
    .thm <-
        theme(axis.text.y = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        theme(legend.title = element_blank()) +
        theme(legend.position = c(1,1)) +
        theme(legend.justification = c(1,1))

    ret <-
        ggplot(.dt, aes(variable, col, fill = value)) +
        ggrastr::rasterise(geom_tile(), dpi=300) +
        scale_fill_distiller(direction=1) +
        .thm
    
    if(with.ct){
        ret <- ret +
            theme(strip.text.y = element_text(angle=0, hjust=0)) +
            facet_grid(celltype ~ ., scales="free", space="free")
    }
    return(ret)
}

################################################################

sort.beta <- function(.dt, genes.selected){

    .genes <- rownames(.dt)
    .mat <- as.matrix(.dt)
    colnames(.mat) <- 1:ncol(.mat)
    rownames(.mat) <- .genes

    .mat <- .mat[rownames(.mat) %in% genes.selected, , drop = F]
    .melt <- reshape2::melt(.mat) %>% as.data.table()
    .melt <- .melt[order(`value`, decreasing = T), head(.SD, 1),
                   by = .(Var1, Var2)] # remove redundancy

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

read.beta <- function(.file, ...){
    return(sort.beta(fread(.file), ...))
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

select.gsea.beta <- function(.gsea, .beta, ntop = 3, nmax = 30, sd.cutoff = 1e-4, padj.cutoff = .1){

    ## cell types (pathways) to show
    .cts <- unlist(.gsea[order(pval),
                         head(.SD, ntop),
                         by = .(Var2)][padj < padj.cutoff,
                                       .(pathway)])

    .gsea.show <- .gsea[pathway %in% .cts]
    .gsea.show[, variable := factor(`Var2`, 0:100)] # 100 should be enough

    .temp <- dcast(.gsea.show, variable ~ pathway,
                   value.var = "pval",
                   fun.aggregate = min)

    .po <- order(apply(-log10(1e-20 + .temp[,-1]), 2, which.max))
    .po <- colnames(.temp)[-1][.po]
    .gsea.show[, ct := factor(pathway, .po)]

    ## genes to show
    .sd <- .beta[, .(s = sd(`value`)), by = .(Var1)]

    .genes.show <- intersect(unlist(.gsea.show$leadingEdge),
                             .sd[s > sd.cutoff]$Var1)
    .beta.show <- .beta[Var1 %in% .genes.show]

    ## leading edges to show
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
