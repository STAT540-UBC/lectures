if(!require(torch)) install.package("torch")
if(!require(mmutilR)) remotes::install_github("causalpathlab/mmutilR")

scdata.loader <- dataset(
    name = "single-cell data loader",
    #' @param .hdr single-cell fileset header
    #' @param MY.DEV torch device
    initialize = function(.hdr, MY.DEV = torch_device("cpu")) {
        self$sc.data <- sc.data <- self$.add.ext(.hdr)
        self$DEV <- MY.DEV
        self$info <- mmutilR::rcpp_mmutil_info(self$sc.data$mtx)
        ncells <- self$info$max.col
        self$mem.idx <- mmutilR::rcpp_mmutil_read_index(sc.data$idx)
        self$read.full()
    },
    read.full = function() {
        .mtx <- self$sc.data$mtx          # matrix market data file
        .idx <- self$mem.idx              # disk memory index file
        .loc <- 1:self$info$max.col       # simply read
        .out <- mmutilR::rcpp_mmutil_read_columns_sparse(.mtx, .idx, .loc)
        .ind <- rbind(.out$col, .out$row) # column is row & vice versa
        .sz <- c(.out$max.col, .out$max.row)
        xx <- torch_sparse_coo_tensor(.ind, .out$val, .sz,
                                      dtype = torch_float())
        xx <- xx$to_dense()
        self$x.list <- lapply(1:nrow(xx), function(r) {
            cat(r,"\r",file=stderr());flush(stderr())
            xx[r, , drop=FALSE]$to_sparse()$to(dtype = torch_float16(),
                                               device=self$DEV)
        })
        message("Load the full data")
    },
    .getitem = function(.loc) {
        torch_cat(self$x.list[.loc], dim=1)$to_dense()$to(dtype=torch_float32(), device=self$DEV)
    },
    .length = function() {
        self$info$max.col
    },
    .dim = function() {
        self$info$max.row
    },
    dim = function() {
        self$.dim()
    },
    .add.ext = function(.hdr){
        .names <- c("mtx","row","col","idx")
        paste0(.hdr, c(".mtx.gz", ".rows.gz", ".cols.gz", ".mtx.gz.index")) %>%
            as.list() %>%
            (function(x) { names(x) <- .names; x })
    },
    .files.exist = function(...){
        lapply(..., file.exists) %>%
            unlist() %>%
            all()
    })
