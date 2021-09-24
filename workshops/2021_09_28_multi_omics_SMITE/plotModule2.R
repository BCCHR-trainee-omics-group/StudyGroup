plotModule2 <- function (pvalue_annotation, p_thresh = 0.05, which_network = 1, 
                         goseq = FALSE, layout = "fr", legend = TRUE, namestyle = "symbol", 
                         suppress_details = FALSE, meth_hi_col = "blue", meth_low_col = "yellow1", 
                         meth_mid_col = "gray90", exp_hi_col = "red1", exp_low_col = "chartreuse1", 
                         exp_mid_col = "gray90", label_scale = TRUE, label_shadow = FALSE, 
                         compare_plot = FALSE, pdf_out = NULL) 
{
  module_output <- slot(slot(pvalue_annotation, "score_data"), 
                        "module_output")
  if (!is.null(pdf_out)) {
    unlink(pdf_out)
    if (compare_plot == TRUE) {
      pdf(pdf_out, width = 24, height = 8.5)
    }
    else {
      pdf(pdf_out, width = 16, height = 8.5)
    }
  }
  if (compare_plot == FALSE) {
    par(mfrow = c(1, 1))
    if (is.null(pdf_out)) {
      if (goseq == FALSE) {
        if (legend == TRUE) {
          if (!names(dev.cur()) %in% c("RStudioGD", "pdf")) {
            dev.new(height = 10, width = 12)
          }
        }
        else {
          if (!names(dev.cur()) %in% c("RStudioGD", "pdf")) {
            dev.new(height = 8, width = 8)
          }
        }
      }
      else {
        if (!names(dev.cur()) %in% c("RStudioGD", "pdf")) {
          dev.new(height = 10, width = 16)
        }
      }
    }
  }
  else {
    legend <- FALSE
    goseq <- FALSE
    suppress_details <- TRUE
    if (!names(dev.cur()) %in% c("RStudioGD", "pdf")) {
      dev.new(height = 10, width = 20)
    }
    par(mfrow = c(1, 2))
  }
  for (n_plot in which_network) {
    if (goseq == TRUE) {
      if (length(module_output$goseqOut) == 0) {
        stop("Goseq analysis has not been performed.")
      }
    }
    name.eid <- names(module_output$modules)[n_plot]
    eid <- module_output$modules[[n_plot]]
    network <- module_output$network
    stat <- extractScores(pvalue_annotation)
    pval <- exp(stat/(-2))
    if (class(network) == "graphNEL") {
      network <- igraph::graph_from_graphnel(network)
      adj_mat_network <- igraph::as_adjacency_matrix(network)
      stat.v <- stat[which(names(stat) %in% rownames(adj_mat_network))]
      stat.v <- stat.v[order(names(stat.v))]
      adj_mat_network <- adj_mat_network[order(rownames(adj_mat_network)), 
                                         order(colnames(adj_mat_network))]
      temp_vstat <- apply(adj_mat_network, 1, function(v) return(v * 
                                                                   stat.v))
      W <- (temp_vstat + t(temp_vstat))/2
      Graph_adj_mat <- igraph::graph_from_adjacency_matrix(W, 
                                                           mode = "undirected", weighted = TRUE)
      igraph::V(Graph_adj_mat)$weight <- stat.v
      network <- Graph_adj_mat
    }
    vect2color <- function(v, palettev, breaks) {
      w <- v
      for (i in 1:length(palettev)) {
        w[which(v >= breaks[i] & v < breaks[i + 1])] = palettev[i]
      }
      return(w)
    }
    vertexPalette.v <- colorRampPalette(c("white", "gray85", 
                                          "gray65", "salmon"))(50)
    vertexBreaks.v <- rev(-2 * log(Hmisc::cut2(pval, g = 50, 
                                               onlycuts = TRUE)))
    vertexBreaks.v[51] <- vertexBreaks.v[51] + 0.001
    edgePalette.v <- colorRampPalette(c("white", "gray85", 
                                        "gray65", "salmon"))(50)
    edgeBreaks.v <- Hmisc::cut2(igraph::E(network)$weight, 
                                g = 50, onlycuts = TRUE)
    h <- igraph::induced_subgraph(network, eid)
    stat.v <- stat[igraph::V(h)$name]
    pval.v <- pval[igraph::V(h)$name]
    pval.v[which(pval.v == 0)] <- 1e-09
    par(mar = c(4, 0, 2, 0))
    igraph::E(h)$color <- vect2color(igraph::E(h)$weight, 
                                     edgePalette.v, edgeBreaks.v)
    igraph::V(h)$color <- vect2color(stat.v, vertexPalette.v, 
                                     vertexBreaks.v)
    vl <- igraph::V(h)$name
    igraph::V(h)$color[which(1 - pchisq(igraph::V(h)$weight, 
                                        2) < p_thresh)] <- "red"
    igraph::E(h)$color[which(1 - pchisq(igraph::E(h)$weight, 
                                        4) < p_thresh)] <- "red"
    if (layout == "circle") {
      layout1 <- igraph::layout_in_circle(h)
    }
    if (layout == "fr") {
      layout1 <- igraph::layout_with_fr(h)
    }
    if (layout == "dh") {
      layout1 <- igraph::layout_with_dh(h)
    }
    if (layout == "kk") {
      layout1 <- igraph::layout_with_kk(h)
    }
    layout1_scaled <- cbind(scales::rescale(layout1[, 1], 
                                            to = c(-1, 1)), scales::rescale(layout1[, 2], to = c(-1, 
                                                                                                 1)))
    counter <- 0
    while (counter < 2) {
      plot(h, layout = layout1, vertex.label = "", vertex.frame.color = "black", 
           vertex.label.dist = 0.1, vertex.label.font = 3, 
           vertex.label.color = "black", vertex.size = if (length(igraph::V(h)) < 
                                                           50) {
             15
           }
           else {
             15 * 13/length(igraph::V(h))
           }, edge.width = if (length(igraph::V(h)) < 50) {
             2
           }
           else {
             1
           }, ylim = c(-1, 1.5), xlim = c(-1, 1))
      halfCircle <- function(x, y, r, r2 = 0.75, quarter = FALSE, 
                             start = 0, end = pi, nsteps = 30, col = NULL, 
                             lwd = 1, border = NULL) {
        if (isTRUE(quarter)) {
          rs <- seq(start, end, len = nsteps)
          xc <- x + r * cos(rs)
          yc <- y + r * sin(rs)
          xc2 <- x + r * r2 * cos(rs)
          yc2 <- y + r * r2 * sin(rs)
          polygon(c(xc, rev(xc2)), c(yc, rev(yc2)), col = col, 
                  lwd = lwd, border = border)
        }
        if (!isTRUE(quarter)) {
          rs <- seq(start, end, len = nsteps)
          xc <- x + r * cos(rs)
          yc <- y + r * sin(rs)
          polygon(xc, yc, col = col, lwd = lwd, border = border)
        }
      }
      arctext <- function(x, y, r, start, end, what, cex = 1) {
        delta = (end - start)/4
        rs <- seq(start, end, len = 3)
        xc <- x + r * cos(rs[2])
        yc <- y + r * sin(rs[2])
        text(xc, yc, what, srt = 180 + atan2(((y + r * 
                                                 sin(rs[3])) - (y + r * sin(rs[1]))), ((x + 
                                                                                          r * cos(rs[3])) - (x + r * cos(rs[1])))) * 
               180/pi, cex = cex)
      }
      methcol <- c(meth_low_col, meth_mid_col, meth_hi_col, "white")
      names(methcol) <- c("Dc", "NS", "In", "NoData")
      expcol <- c(exp_low_col, exp_mid_col, exp_hi_col, "white")
      names(expcol) <- c("Down", "NS", "Up", "NoData")
      pval_data <- slot(slot(pvalue_annotation, "score_data"), 
                        "pval_data")
      genes_score <- slot(pvalue_annotation, "score_data")@genes
      effect_data <- slot(slot(pvalue_annotation, "score_data"), 
                          "effect_data")
      signs_idx <- slot(slot(pvalue_annotation, "score_data"), 
                        "signs_index")
      if (any(suppress_details == FALSE, counter == 1)) {
        for (i in 1:nrow(layout1_scaled)) {
          halfCircle(x = layout1_scaled[i, 1], y = layout1_scaled[i, 
                                                                  2], r = ifelse(length(igraph::V(h)) < 50, 
                                                                                 0.075, 0.025), start = pi/2, end = 2 * pi/2, 
                     quarter = TRUE, lwd = 1, col = ifelse(!is.na(pval_data$expression_pvalue[which(genes_score %in% 
                                                                                                      igraph::V(h)$name[i])]), expcol[ifelse(abs(pval_data$expression_pvalue[which(genes_score %in% 
                                                                                                                                                                                     igraph::V(h)$name[i])]) < p_thresh, ifelse(sign(effect_data$expression_effect[which(genes_score %in% 
                                                                                                                                                                                                                                                                           igraph::V(h)$name[i])]) == 1, 3, 1), 2)], 
                                                           expcol[4]))
          start <- pi
          delta <- (3 * pi/2)/nrow(signs_idx)
          for (j in signs_idx[, 3]) {
            score_graph_col <- returnPvalueCol(slot(pvalue_annotation, 
                                                    "score_data"), j)[which(genes_score %in% 
                                                                              igraph::V(h)$name[i])]
            halfCircle(x = layout1_scaled[i, 1], y = layout1_scaled[i, 
                                                                    2], r = ifelse(length(igraph::V(h)) < 50, 
                                                                                   0.075, 0.025), start = start, end = start + 
                         delta, quarter = TRUE, col = ifelse(!is.na(score_graph_col), 
                                                             methcol[ifelse(abs(score_graph_col) < p_thresh, 
                                                                            ifelse(sign(effect_data[, grep(j, colnames(effect_data))][which(genes_score %in% 
                                                                                                                                              igraph::V(h)$name[i])]) == 1, 3, 1), 
                                                                            2)], methcol[4]))
            start <- start + delta
          }
        }
      }
      if (any(legend == TRUE, counter == 1)) {
        if (any(suppress_details == FALSE, counter == 
                1)) {
          num_factors <- nrow(signs_idx)
          halfCircle(x = -1.25, y = 1.25, r = 0.4, start = pi/2, 
                     end = pi, quarter = TRUE)
          start <- pi/2
          delta <- pi/8
          for (g in 1:4) {
            halfCircle(x = -1.25, y = 1.25, r = 0.37, 
                       r2 = 0.89, start = start + delta * (g - 
                                                             1), end = start + delta * g, col = expcol[g], 
                       quarter = TRUE)
            arctext(x = -1.25, y = 1.25, r = 0.35, start = start + 
                      delta * (g - 1), end = start + delta * 
                      g, names(expcol)[g], cex = 0.5)
          }
          text(-1.25 + 1.2 * cos(seq(start, start + delta * 
                                       2, len = 30)[15]), 1.25 + 0.4 * 0.75 * sin(seq(start, 
                                                                                      start + delta * 2, len = 30)[15]), "expression", 
               cex = 0.6)
          start <- pi
          delta <- (3 * pi/2)/num_factors
          for (j in signs_idx[, 3]) {
            halfCircle(x = -1.25, y = 1.25, r = 0.4, 
                       start = start, end = start + delta, quarter = TRUE)
            for (g in 1:4) {
              start <- start
              end <- start + delta
              delta2 <- (end - start)/4
              halfCircle(x = -1.25, y = 1.25, r = 0.37, 
                         r2 = 0.89, start = start + delta2 * (g - 
                                                                1), end = start + delta2 * g, col = methcol[g], 
                         quarter = TRUE)
              arctext(x = -1.25, y = 1.25, r = 0.35, 
                      start = start + delta2 * (g - 1), end = start + 
                        delta2 * g, ifelse(num_factors <= 4, 
                                           names(methcol)[g], substring(names(methcol)[g], 
                                                                        1, 1)), cex = 0.5)
            }
            text(-1.25 + 0.57 * cos(seq(start, start + 
                                          delta, len = 30)[15]), 1.25 + 0.57 * 0.75 * 
                   sin(seq(start, start + delta, len = 30)[15]), 
                 paste(strsplit(j, "_")[[1]], collapse = "\n"), 
                 cex = 0.6)
            start <- start + delta
          }
        }
        halfCircle(x = -1.25, y = 1.25, r = 0.3, end = 2 * 
                     pi, col = "white")
        vertexPalette.v[which(vertexBreaks.v >= min(igraph::V(h)$weight[which(igraph::V(h)$weight >= 
                                                                                qchisq(1 - p_thresh, 2))])) - 1] <- "red"
        points(seq(-1.45, -1.05, length.out = 50), rep(1.35, 
                                                       50), col = vertexPalette.v, pch = 15, cex = 2.5)
        text(-1.2, 1.44, expression("node             ", 
                                    Chi[2]^2))
        text(-1.45, 1.295, round(vertexBreaks.v[1], 2))
        text(-1.25, 1.295, round(vertexBreaks.v[26], 
                                 2))
        text(-1.07, 1.295, paste(">", round(qchisq(1 - 
                                                     p_thresh, 2), 2), sep = ""))
        if (length(which(igraph::E(h)$weight >= qchisq(1 - 
                                                       p_thresh, 4))) > 0) {
          edgePalette.v[which(edgeBreaks.v >= min(igraph::E(h)$weight[which(igraph::E(h)$weight >= 
                                                                              qchisq(1 - p_thresh, 4))])) - 1] <- "red"
        }
        points(seq(-1.45, -1.05, length.out = 50), rep(1.15, 
                                                       50), col = edgePalette.v, pch = 15, cex = 2.5)
        text(-1.2, 1.24, expression("edge             ", 
                                    Chi[4]^2))
        text(-1.45, 1.095, round(edgeBreaks.v[1], 2))
        text(-1.25, 1.095, round(edgeBreaks.v[26], 2))
        text(-1.07, 1.095, paste(">", round(qchisq(1 - 
                                                     p_thresh, 4), 2), sep = ""))
      }
      if (label_shadow == T) {
        addShadoA30Pext(layout1_scaled[, 1], layout1_scaled[, 
                                                            2], vl, font = 2, cex = if (label_scale == 
                                                                                        TRUE) {
                                                              scales::rescale(stat.v, to = (c(0.5, 2)))
                                                            }
                        else {
                          0.5
                        }, bg = "white", col = "black")
      }
      if (label_shadow == F) {
        text(layout1_scaled[, 1], layout1_scaled[, 2], 
             vl, font = 2, cex = if (label_scale == TRUE) {
               scales::rescale(stat.v, to = (c(0.5, 2)))
             }
             else {
               0.5
             }, bg = "white", col = "black")
      }
      if (namestyle == "refseq") {
        ref2eg <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egREFSEQ2EG)
        eg2sym <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
        text(layout1_scaled[, 1], layout1_scaled[, 2] - 
               0.05, sapply(vl, function(k) {
                 ifelse(is.null(ref2eg[[k]]), return(NA), return(eg2sym[[ref2eg[[k]]]]))
               }))
        text(0, 1.65, paste("Network built around", name.eid, 
                            ifelse(is.null(ref2eg[[name.eid]]), NA, eg2sym[[ref2eg[[name.eid]]]])))
      }
      if (namestyle == "symbol") {
        text(0, 1.7, paste("Network built around", name.eid, 
                           "\nChi-square P-value=", round(module_output$moduleStats[[n_plot]][2], 
                                                          4)))
      }
      if (goseq == TRUE) {
        if (nrow(module_output$goseqOut[[n_plot]]) > 
            0) {
          text("Num\nGenes", x = 1.35, y = 1.6, font = 2)
          text("Enriched\nPathway/Term", x = 2, y = 1.6, 
               font = 2)
          for (i in 1:nrow(module_output$goseqOut[[n_plot]])) {
            text(module_output$goseqOut[[n_plot]][i, 
                                                  4], x = 1.3, y = seq(1.4, -1, length.out = nrow(module_output$goseqOut[[n_plot]]))[i], 
                 adj = c(0, 0))
            text(module_output$goseqOut[[n_plot]][i, 
                                                  6], x = 1.5, y = seq(1.4, -1, length.out = nrow(module_output$goseqOut[[n_plot]]))[i], 
                 adj = c(0, 0))
          }
        }
        else {
          text("No enriched terms from\nGoseq", x = 1.5, 
               y = 1.6, font = 2)
        }
      }
      if (counter == 1) {
        break
      }
      counter <- 2
    }
    if (is.null(pdf_out)) {
      if (n_plot != which_network[length(which_network)]) {
        message("Press key to go to next plot")
        readline()
      }
    }
  }
  if (!is.null(pdf_out)) {
    dev.off()
  }
}


returnPvalueCol <- function(pval_object, col_name){
  
  return(slot(pval_object,
              "pval_data")[, grep(col_name,
                                  colnames(slot(pval_object,
                                                "pval_data")))])
}