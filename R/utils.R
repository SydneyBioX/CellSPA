
#' Subset Spatial experiment object by column
#'
#'
#' @param spe A Spatial Experiment object
#' @param col_idx The column index to keep.
#'
#' @return A SpatialExperiment object
#' @export
subset <- function(spe, col_idx) {
    spe <- spe[, col_idx]
    if (!is.null(spe@metadata$CellSegOutput)) {
        keep <- spe@metadata$CellSegOutput$cell_id %in% spe$cell_id
        spe@metadata$CellSegOutput <- spe@metadata$CellSegOutput[keep, ]
    }

    return(spe)
}

.add_metrics <- function(spe, metrics_name, type) {
    spe@metadata$CellSPA$metrics[[type]] <- sort(unique(c(spe@metadata$CellSPA$metrics[[type]],
                                                          metrics_name)))
    return(spe)
}


.add_dataset_metrics <- function(spe) {
    spe@metadata$CellSPA$dataset_level_metrics <- list(num_cells = ncol(spe),
                                                       num_genes = nrow(spe))
    spe <- .add_metrics(spe, c("num_cells", "num_genes"), "dataset_level")
    return(spe)
}


.initialise_CellSPA_list <- function(spe) {
    spe@metadata$CellSPA <- list()
    spe@metadata$CellSPA$metrics <- vector("list", 4)
    names(spe@metadata$CellSPA$metrics) <- c("dataset_level",
                                             "gene_level",
                                             "cell_level",
                                             "bin_level")
    return(spe)
}


# From Matrix.utils::aggregate.Matrix() function, which has been archived on CRAN

#' @export
aggregate.Matrix <- function(x, groupings = NULL, form = NULL, fun = 'sum', ...) {
    if (!methods::is(x,'Matrix')) {
        x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
    }
    if (fun=='count') {
        x <- x != 0
    }
    groupings2 <- groupings
    if (!methods::is(groupings2,'data.frame')) {
        groupings2 <- as(groupings2,'data.frame')
    }

    groupings2 <- data.frame(lapply(groupings2, as.factor))
    groupings2 <- data.frame(interaction(groupings2,sep = '_'))
    colnames(groupings2) <- 'A'
    if(is.null(form)) {
        form <- stats::as.formula('~0+.')
    }

    form <- stats::as.formula(form)
    mapping <- dMcast(groupings2,form)
    colnames(mapping) <- substring(colnames(mapping),2)
    result <- Matrix::t(mapping) %*% x


    attr(result,'crosswalk') <- grr::extract(groupings,match(rownames(result),groupings2$A))
    return(result)
}

# From Matrix.utils::dMcast() function, which has been archived on CRAN

dMcast <- function(data, formula, fun.aggregate='sum',
                   value.var=NULL, as.factors=FALSE,
                   factor.nas=TRUE, drop.unused.levels=TRUE) {
    values <- 1
    if( !is.null(value.var)) {
        values<-data[,value.var]
    }

    alltms <- stats::terms(formula, data = data)
    response <- rownames(attr(alltms,'factors'))[attr(alltms,'response')]
    tm <- attr(alltms,"term.labels")
    interactionsIndex <- grep(':',tm)
    interactions <- tm[interactionsIndex]
    simple <- setdiff(tm,interactions)
    i2 <- strsplit(interactions,':')
    newterms <- unlist(lapply(i2, function (x) {
        paste("paste(", paste(x,collapse=','),",","sep='_'",")")
    }))
    newterms <- c(simple,newterms)
    newformula <- stats::as.formula(paste('~0+',paste(newterms,collapse='+')))
    allvars <- all.vars(alltms)
    data <- data[,c(allvars),drop=FALSE]
    if (as.factors) {
        data <- data.frame(lapply(data,as.factor))
    }

    characters <- unlist(lapply(data,is.character))
    data[,characters] <- lapply(data[,characters,drop = FALSE], as.factor)
    factors <- unlist(lapply(data, is.factor))
    #Prevents errors with 1 or fewer distinct levels
    data[, factors] <- lapply(data[, factors, drop = FALSE], function (x) {
        if(factor.nas) {
            if(any(is.na(x))) {
                levels(x) <- c(levels(x),'NA')
                x[is.na(x)] <- 'NA'
            }
        }

        if (drop.unused.levels) {
            if(nlevels(x) != length(stats::na.omit(unique(x)))) {
                x <- factor(as.character(x))
            }
        }


        y <- stats::contrasts(x, contrasts = FALSE, sparse = TRUE)
        attr(x,'contrasts') <- y
        return(x)
    })
    #Allows NAs to pass
    attr(data,'na.action') <- stats::na.pass
    result <- Matrix::sparse.model.matrix(newformula, data, drop.unused.levels = FALSE, row.names=FALSE)
    brokenNames <- grep('paste(',colnames(result),fixed = TRUE)
    colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames], function (x) {
        x <- gsub('paste(',replacement = '', x = x, fixed = TRUE)
        x <- gsub(pattern = ', ',replacement='_', x = x, fixed = TRUE)
        x <- gsub(pattern = '_sep = \"_\")', replacement = '', x=x, fixed = TRUE)
        return(x)
    })

    result <- result*values
    # if (isTRUE(response>0)) {
    #     responses = all.vars(terms(as.formula(paste(response, '~0'))))
    #     result <- aggregate.Matrix(result, data[,responses,drop=FALSE],
    #                                fun = fun.aggregate)
    # }
    return(result)
}

