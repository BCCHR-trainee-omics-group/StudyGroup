annotateModification2 <- function (pvalue_annotation, mod_data, weight_by = NULL, weight_by_method = "Stouffer",
          mod_included = NULL, mod_corr = TRUE, mod_type = "methylation",
          verbose = FALSE)
{ 
  if (mod_type %in% names(slot(slot(pvalue_annotation, "modifications"),
                               "metadata")$elements)) {
    stop("Provided data set is already loaded as mod_type")
  }
  unique_feature_names <- unique(unlist(slot(pvalue_annotation,
                                             "annotation"))$feature)
  if (missing(weight_by)) {
    weight_by <- rep("pval", length(unique_feature_names[!unique_feature_names %in%
                                                           c("original", "tss")]))
  }
  if (is.null(names(weight_by))) {
    if (is.null(mod_included)) {
      mod_included <- unique_feature_names[!unique_feature_names %in%
                                             c("original", "tss")]
    }
    names(weight_by) <- mod_included
  }
  if (!is.null(names(weight_by))) {
    mod_included <- names(weight_by)
    if (!all(mod_included %in% unique_feature_names)) {
      stop("Provided weight names must match those in\n                     unique(GenomicRanges::mcols(unlist(pvalue_annotation@annotation))$feature)")
    }
  }
  if (any(!c(-1, 1) %in% unique(sign(mod_data[, 4])))) {
    message("WARNING: Effects should provide a direction,\n                    but these effects are all in the same direction.")
  }
  if (any(mod_data[, 5] < 0, mod_data[, 5] > 1)) {
    stop("P-values must be between 0 and 1")
  }
  mod_grange <- GenomicRanges::GRanges(seqnames = mod_data[,
                                                           1], ranges = IRanges::IRanges(start = mod_data[, 2],
                                                                                         end = mod_data[, 3]), effect = mod_data[, 4], pval = mod_data[,
                                                                                                                                                       5], type = mod_type)
  temp_annotation <- unlist(slot(pvalue_annotation, "annotation"))
  overlap_mods <- GenomicRanges::findOverlaps(temp_annotation,
                                              mod_grange)
  mod_grange_overlaps <- mod_grange[S4Vectors::subjectHits(overlap_mods)]
  GenomicRanges::mcols(mod_grange_overlaps) <- cbind(GenomicRanges::mcols(temp_annotation[as.numeric(S4Vectors::queryHits(overlap_mods))]),
                                                     GenomicRanges::mcols(mod_grange_overlaps))
  mod_grange_overlaps <- split(mod_grange_overlaps, mod_grange_overlaps$name)
  temp_annotation <- split(temp_annotation, temp_annotation$name)
  if (mod_corr == TRUE) {
    if (verbose == TRUE) {
      message("Computing correlation matrices")
    }
    temp_split_mod_grange <- split(mod_grange, GenomicRanges::seqnames(mod_grange))
    precede_follow_each_element <- lapply(temp_split_mod_grange,
                                          function(chr) {
                                            temp_chr <- IRanges(start(chr), end(chr))
                                            temp_precede <- precede(temp_chr)
                                            temp_follow <- follow(temp_chr)
                                            temp_precede[which(is.na(temp_precede))] <- which(is.na(temp_precede))
                                            temp_follow[which(is.na(temp_follow))] <- which(is.na(temp_follow))
                                            chr[c(temp_follow, temp_precede)]
                                          })
    mod_grange_corr <- unlist(GRangesList(precede_follow_each_element))
    duplicate_each_chr <- lapply(temp_split_mod_grange, function(chr) {
      c(chr, chr)
    })
    duplicate_each_chr <- unlist(GRangesList(duplicate_each_chr))
    mod_grange_corr$distance <- IRanges::distance(duplicate_each_chr,
                                                  mod_grange_corr)
    mod_grange_corr$pval2 <- duplicate_each_chr$pval
    mod_grange_corr <- mod_grange_corr[which(mod_grange_corr$pval2 <
                                               0.05)]
    quantile_distances_mod_corr <- Hmisc::cut2(mod_grange_corr$distance,
                                               g = 500, onlycuts = TRUE)
    quantile_distances_mod_corr[length(quantile_distances_mod_corr)] <- 2.5e+08
    mod_grange_corr$cat <- cut(mod_grange_corr$distance,
                               breaks = quantile_distances_mod_corr)
    mod_grange_corr <- split(mod_grange_corr, mod_grange_corr$cat)
    mod_grange_corr2 <- lapply(mod_grange_corr, function(j) {
      mean((sapply(1:500, function(i) {
        index <- sample(1:length(j), replace = TRUE)
        cor(qnorm(1 - j$pval[index]), qnorm(1 - j$pval2[index]))
      })))
    })
    correlations <- as.data.frame(do.call(rbind, mod_grange_corr2))
    final_corr <- data.frame(correlations, as.character(names(mod_grange_corr2)),
                             stringsAsFactors = FALSE)
    final_corr <- rbind(c(0.9, paste("(-1, ", quantile_distances_mod_corr[1],
                                     "]", sep = "")), final_corr)
    rm(mod_grange_corr)
  }
  combined_pvalues_list <- sapply(mod_included, function(i) {
    if (verbose == TRUE) {
      message(paste("Combining p-values over:", i))
    }
    temp <- subset(unlist(mod_grange_overlaps), unlist(mod_grange_overlaps)$feature ==
                     i)
    ref_data <- unlist(slot(pvalue_annotation, "annotation"))
    ref_data <- subset(ref_data, ref_data$feature == "tss")
    ref_data <- ref_data[temp$name]
    suppressWarnings(temp$distance <- distance(ref_data,
                                               temp) + 2)
    temp <- split(temp, temp$name)
    forreturn <- lapply(temp, function(each) {
      each_length <- length(each)
      each_effect <- each$effect[order(each$pval)]
      each_pval <- each$pval[order(each$pval)]
      distances <- each$distance[order(each$pval)]
      if (length(each_pval) > 1) {
        if (mod_corr == TRUE) {
          corr_mat <- matrix(as.numeric(final_corr[match(cut(as.matrix(dist(start(each)[order(each$pval)])),
                                                             breaks = c(-1, quantile_distances_mod_corr)),
                                                         final_corr[, 2]), 1]), ncol = each_length)
          diag(corr_mat) <- 1
          corr_mat <- abs(corr_mat)
          chol_d <- try(chol(corr_mat), silent = TRUE)
          while (is(chol_d, "try-error")) {
            index <- as.numeric(strsplit(strsplit(chol_d[1],
                                                  "the leading minor of order ")[[1]][2],
                                         " is not positive")[[1]][1]) - 1
            chol_d <- try(chol(corr_mat[1:index, 1:index]),
                          silent = TRUE)
            each_pval <- each_pval[1:index]
            each_effect <- each_effect[1:index]
            distances <- distances[1:index]
            each_length <- index
          }
          each_pval <- 1 - pnorm(abs(solve(t(chol_d)) %*%
                                       qnorm(1 - each_pval/2)))
          each_pval <- replace(each_pval, each_pval ==
                                 0, 1e-09)
          each_pval <- each_pval * 2
        }
        if (weight_by_method == "Stouffer") {
          if (weight_by[i] == "distance") {
            out_mean <- weighted.mean(each_effect, w = (1/log(distances)))
            out_pval <- stoufferTest(each_pval, weights = (1/log(distances)))
          }
          else if (weight_by[i] %in% c("pval", "p.value",
                                       "pvalue", "p_val")) {
            if (length(unique(each_pval))==1){
            if (unique(each_pval)==1){
              each_pval <- as.numeric(gsub(1, 0.9999999999, each_pval))
            }}
            out_mean <- weighted.mean(each_effect, w = -log(each_pval))
            out_pval <- stoufferTest(each_pval, weights = NULL)
          }
          else {
            out_mean <- mean(each_effect, na.rm = TRUE)
            out_pval <- stoufferTest(each_pval, weights = NULL)
          }
        }
        else if (weight_by_method %in% c("minimum", "Sidak",
                                         "sidak")) {
          index <- which(each_pval == min(each_pval))
          if (length(index) > 1) {
            index <- index[which(abs(each_effect[index]) ==
                                   max(abs(each_effect[index])))][1]
          }
          out_mean <- each_effect[index]
          out_pval <- 1 - (1 - each_pval[index])^length(each_pval)
        }
        else if (weight_by_method == "binomial") {
          index <- which(each_pval == min(each_pval))
          if (length(index) > 1) {
            index <- index[which(abs(each_effect[index]) ==
                                   max(abs(each_effect[index])))][1]
          }
          out_mean <- each_effect[index]
          out_pval <- (1 - pbinom(q = length(which(each_pval <
                                                     0.05)), size = each_length, prob = 0.05))
        }
        else if (weight_by_method %in% c("Fisher", "fisher",
                                         "chisq", "chi")) {
          out_pval <- 1 - pchisq(-2 * sum(log(each_pval)),
                                 each_length * 2)
          out_mean <- mean(sign(each_effect), na.rm = TRUE)
        }
      }
      else {
        out_mean <- each_effect
        out_pval <- each_pval
      }
      c(out_mean, out_pval, each_length)
    })
    do.call(rbind, forreturn)
  })
  if (verbose == TRUE) {
    message("Quantile permuting scores")
  }
  combined_pvalues_list <- lapply(combined_pvalues_list, function(each_feature) {
    categories <- data.frame(categories = as.numeric(Hmisc::cut2(each_feature[,
                                                                              3], g = 100)))
    categories_table <- data.frame(table(categories))
    trans_p <- cbind(trans = qnorm(1 - each_feature[, 2]/2),
                     plyr::join(categories, categories_table, by = "categories"))
    trans_p[, 1] <- replace(trans_p[, 1], is.infinite(trans_p[,
                                                              1]), max(subset(trans_p, !is.infinite(trans_p[, 1])),
                                                                       na.rm = TRUE))
    num_list <- split(trans_p$trans, trans_p$categories)
    rand_list <- sapply(1:length(num_list), function(i) {
      as.matrix(sapply(1:500, function(j) {
        sample(num_list[[as.numeric(i)]], replace = TRUE)
      }))
    })
    new_pval <- apply(trans_p, 1, function(i) {
      length(which(rand_list[[as.numeric(i[2])]] > as.numeric(i[1])))/(500 *
                                                                         as.numeric(i[3]))
    })
    new_pval <- replace(new_pval, new_pval == 0, min(subset(new_pval,
                                                            new_pval != 0), na.rm = TRUE))
    each_feature[, 2] <- new_pval
    each_feature <- as.data.frame(each_feature)
    each_feature
  })
  if (verbose == TRUE) {
    message("Scores have been adjusted")
  }
  newmods <- 
               unlist(mod_grange_overlaps)
  names(newmods) <- NULL
  newmods <- split(newmods, newmods$name)
  output_m_summary <- suppressWarnings(as.data.frame(c(list(names = names(mod_grange_overlaps)),
                                                       lapply(combined_pvalues_list, function(x) {
                                                         x[match(names(mod_grange_overlaps), rownames(x)),
                                                           1:2]
                                                       }))))
  rownames(output_m_summary) <- output_m_summary[, 1]
  output_m_summary <- output_m_summary[, -1]
  colnames(output_m_summary) <- paste(mod_type, apply(expand.grid(c("effect",
                                                                    "pvalue"), mod_included), 1, function(i) {
                                                                      paste(i[2], i[1], sep = "_")
                                                                    }), sep = "_")
  newmetadata <- slot(slot(pvalue_annotation, "modifications"),
                      "metadata")
  if (is.null(newmetadata$m_summary)) {
    newmetadata$m_summary <- output_m_summary
  }
  else {
    newmetadata$m_summary <- merge(newmetadata$m_summary,
                                   output_m_summary, by = 0, all = TRUE)
    rownames(newmetadata$m_summary) <- newmetadata$m_summary[,
                                                             1]
    newmetadata$m_summary <- newmetadata$m_summary[, -1]
  }
  newmetadata[["elements"]][[mod_type]]$weight_by <- weight_by
  newmetadata[["elements"]][[mod_type]]$weight_by_method <- weight_by_method
  newmetadata$elementnames <- c(newmetadata$elementnames, paste(mod_type,
                                                                mod_included, sep = "_"))
  slot(newmods, "metadata") <- newmetadata
  slot(pvalue_annotation, "modifications") <- newmods
  pvalue_annotation
}  

  
                                         
    