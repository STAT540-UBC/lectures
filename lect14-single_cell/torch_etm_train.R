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

#' @param .mean mean vector
#' @param .lnvar log-variance (independent)
normal.stoch <- function(.mean, .lnvar) {
    .eps <- torch_randn_like(.lnvar)
    .sig <- torch_exp(.lnvar / 2.)
    .mean + .eps * .sig
}

kl.loss <- function(.mean, .lnvar) {
    -0.5 * torch_sum(1. + .lnvar - torch_pow(.mean, 2.) - torch_exp(.lnvar), dim = -1);
}

#' FC ReLU layer generator
#' @param in_ input dimension
#' @param out_ output dimension
#' @return a generator
#'
build.fc.relu <- function(in_, out_) {
    nn_sequential(nn_linear(in_, out_),
                  nn_relu())
}

#' Build a stack of layers
#' @param in_ input dimension (first layer)
#' @param layers dimensions for the subsequent layers
#' @param generator a shared generator function
#' @param name a header for the names
#'
build.stack <- function(in_, layers, generator = nn_linear, name = "", stdizer = NULL) {
    ret <-
        nn_module(
            classname = "stack.layers",
            initialize = function(in_,
                                  .layers,
                                  .generator,
                                  .name) {
                d.prev <- in_
                for(l in 1:length(.layers)) {
                    d.this <- .layers[l]
                    .name.l <- paste0(name, ".a", (l - 1))
                    self$add_module(.name.l, module = .generator(d.prev, d.this))
                    if(!is.null(stdizer)){
                        .name.l <- paste0(name, ".d", (l - 1))
                        self$add_module(.name.l, module = stdizer(d.this))
                    }
                    d.prev <- d.this
                }
            },
            forward = function(input) {
                for(module in private$modules_) {
                    input <- module(input)
                }
                input
            })
    return(ret(in_, layers, generator, name))
}

#' @param n.in
#' @param K
#' @param d
build.etm.encoder <-
    nn_module(
        classname = "Encoder",
        initialize = function(n.in, K, layers) {
            self$fc <- build.stack(n.in, layers,
                                   generator = build.fc.relu,
                                   name = "enc.fc")
            d <- layers[length(layers)]
            self$K <- K
            self$z.mean <- nn_linear(d, K)
            self$z.lnvar <- nn_linear(d, K)
        },
        forward = function(xx) {
            xx <- nnf_normalize(xx, dim=2)
            ss <- self$fc(xx)
            mm <- self$z.mean(ss)
            lv <- torch_clamp(self$z.lnvar(ss), -4.0, 4.0)
            z <- normal.stoch(mm, lv)
            list(z = z, z.mean = mm, z.lnvar = lv)
        })

#' @param n.out
#' @param K
#' @param jitter
build.etm.decoder <-
    nn_module(
        classname = "Decoder",
        initialize = function(n.out, K, jitter = .1) {
            self$lbeta <- nn_parameter(torch_randn(K, n.out) * jitter)
            self$beta <- nn_log_softmax(2) # topic x variant (softmax for each topic)
            self$hid <- nn_log_softmax(2)  # sample x topic (softmax for each sample)
        },
        forward = function(zz) {
            .beta <- self$beta(self$lbeta)
            .hh <- self$hid(zz)
            torch_mm(torch_exp(.hh), torch_exp(.beta))
        })

#' @param n.in input data dimension
#' @param n.out output data dimension
#' @param K latent state dimension
build.etm <-
    nn_module(
        classname = "ETM",
        initialize = function(n.in, n.out, K, d = 16) { # Register model parametesr
            self$enc <- build.etm.encoder(n.in, K, layers=d)   # encoder model
            self$dec <- build.etm.decoder(n.out, K)     # decoder model
        },
        forward = function(xx) {                        # Define the forward pass
            .enc <- self$enc(xx)                        # Converts data to the latent
            .pr <- self$dec(.enc$z)                     # Converts the latent to recon
            list(recon = .pr, z.mean = .enc$z.mean, z.lnvar = .enc$z.lnvar)
        })


#' @param xx observed data (vocabulary frequency)
#' @param pr reconstructed probability of observing vocabulary
etm.llik <- function(xx, pr, eps = 1e-8){
    torch_sum(xx * torch_log(pr + eps), dim=-1)
}
