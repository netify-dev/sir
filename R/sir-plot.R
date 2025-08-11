# Avoid R CMD check warnings about non-standard evaluation in ggplot2
utils::globalVariables(c("Sender", "Receiver", "Value", "density", 
                         "Iteration", "Deviance", "Estimate", "Coefficient",
                         "Lower", "Upper", "Type", "original_weight", "name"))

#' Plot Methods for SIR Objects using ggplot2
#' 
#' @description
#' Creates diagnostic and visualization plots for SIR models using ggplot2,
#' including influence matrix heatmaps, distributions, convergence diagnostics,
#' and network visualizations.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_point geom_hline 
#' @importFrom ggplot2 geom_smooth scale_fill_gradient2 scale_fill_viridis_c
#' @importFrom ggplot2 theme_minimal theme labs element_text element_blank
#' @importFrom ggplot2 stat_qq stat_qq_line facet_wrap coord_fixed
#' @importFrom ggplot2 geom_density geom_vline geom_line geom_errorbarh 
#' @importFrom ggplot2 scale_color_manual after_stat
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom tidyr pivot_longer
#' @importFrom viridis scale_fill_viridis
#' @importFrom stats median density reorder sd
#' @importFrom grDevices dev.off pdf

#' Plot Diagnostics for SIR Model using ggplot2
#' 
#' @param x A sir object
#' @param which Which plots to produce (1-6):
#'   1 = Influence Matrix A heatmap
#'   2 = Influence Matrix B heatmap  
#'   3 = Distribution of A matrix values
#'   4 = Distribution of B matrix values
#'   5 = Network statistics over time (if available)
#'   6 = Coefficient plot with confidence intervals
#' @param combine Logical, whether to combine plots using patchwork
#' @param title Main title for combined plot
#' @param theme_base Base ggplot2 theme to use
#' @param ... Additional arguments (unused)
#' @return A ggplot2 object or patchwork combination
#' @export
plot.sir <- function(x, 
                     which = 1:4,
                     combine = TRUE,
                     title = NULL,
                     theme_base = theme_minimal(),
                     ...) {
  
  plots <- list()
  
  # Plot 1: Influence Matrix A Heatmap
  if (1 %in% which && !is.null(x$A)) {
    # Convert matrix to long format for ggplot
    A_df <- expand.grid(
      Sender = factor(1:nrow(x$A)),
      Receiver = factor(1:ncol(x$A))
    )
    A_df$Value <- as.vector(x$A)
    # Set diagonal to NA for better visualization
    A_df$Value[A_df$Sender == A_df$Receiver] <- NA
    
    p1 <- ggplot(A_df, aes(x = Sender, y = Receiver, fill = Value)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "#2166AC", 
        mid = "white", 
        high = "#B2182B",
        midpoint = 0,
        name = "Influence",
        na.value = "grey90"
      ) +
      coord_fixed() +
      labs(
        title = "Influence Matrix A (Sender Effects)",
        subtitle = "How nodes influence others",
        x = "Sender Node",
        y = "Receiver Node"
      ) +
      theme_base +
      theme(
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    plots$A_heatmap <- p1
  }
  
  # Plot 2: Influence Matrix B Heatmap
  if (2 %in% which && !is.null(x$B)) {
    B_df <- expand.grid(
      Sender = factor(1:nrow(x$B)),
      Receiver = factor(1:ncol(x$B))
    )
    B_df$Value <- as.vector(x$B)
    B_df$Value[B_df$Sender == B_df$Receiver] <- NA
    
    p2 <- ggplot(B_df, aes(x = Sender, y = Receiver, fill = Value)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "#5E4FA2", 
        mid = "white", 
        high = "#F46D43",
        midpoint = 0,
        name = "Influence",
        na.value = "grey90"
      ) +
      coord_fixed() +
      labs(
        title = "Influence Matrix B (Receiver Effects)",
        subtitle = "How nodes are influenced by others",
        x = "Sender Node",
        y = "Receiver Node"
      ) +
      theme_base +
      theme(
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    plots$B_heatmap <- p2
  }
  
  # Plot 3: Distribution of A matrix (off-diagonal)
  if (3 %in% which && !is.null(x$A)) {
    A_offdiag <- x$A[row(x$A) != col(x$A)]
    A_dist_df <- data.frame(Value = A_offdiag)
    
    p3 <- ggplot(A_dist_df, aes(x = Value)) +
      geom_histogram(
        aes(y = after_stat(density)), 
        bins = 30, 
        fill = "#2166AC", 
        alpha = 0.7,
        color = "white"
      ) +
      geom_density(color = "#B2182B", linewidth = 1) +
      geom_vline(
        aes(xintercept = mean(Value)), 
        color = "#B2182B", 
        linetype = "dashed", 
        linewidth = 1
      ) +
      geom_vline(
        aes(xintercept = median(Value)), 
        color = "#2166AC", 
        linetype = "dashed", 
        linewidth = 1
      ) +
      labs(
        title = "Distribution of A Matrix Values",
        subtitle = paste0("Mean: ", round(mean(A_offdiag), 4), 
                         " | Median: ", round(median(A_offdiag), 4)),
        x = "Influence Value",
        y = "Density"
      ) +
      theme_base +
      theme(
        panel.grid.minor = element_blank()
      )
    
    plots$A_dist <- p3
  }
  
  # Plot 4: Distribution of B matrix (off-diagonal)
  if (4 %in% which && !is.null(x$B)) {
    B_offdiag <- x$B[row(x$B) != col(x$B)]
    B_dist_df <- data.frame(Value = B_offdiag)
    
    p4 <- ggplot(B_dist_df, aes(x = Value)) +
      geom_histogram(
        aes(y = after_stat(density)), 
        bins = 30, 
        fill = "#5E4FA2", 
        alpha = 0.7,
        color = "white"
      ) +
      geom_density(color = "#F46D43", linewidth = 1) +
      geom_vline(
        aes(xintercept = mean(Value)), 
        color = "#F46D43", 
        linetype = "dashed", 
        linewidth = 1
      ) +
      geom_vline(
        aes(xintercept = median(Value)), 
        color = "#5E4FA2", 
        linetype = "dashed", 
        linewidth = 1
      ) +
      labs(
        title = "Distribution of B Matrix Values",
        subtitle = paste0("Mean: ", round(mean(B_offdiag), 4), 
                         " | Median: ", round(median(B_offdiag), 4)),
        x = "Influence Value",
        y = "Density"
      ) +
      theme_base +
      theme(
        panel.grid.minor = element_blank()
      )
    
    plots$B_dist <- p4
  }
  
  # Plot 5: Deviance/convergence history (if available)
  if (5 %in% which && !is.null(x$history) && !is.null(x$history$DEV)) {
    dev_history <- x$history$DEV
    if (nrow(dev_history) > 1) {
      dev_df <- data.frame(
        Iteration = 1:nrow(dev_history),
        Deviance = dev_history[, 2]
      )
      
      p5 <- ggplot(dev_df, aes(x = Iteration, y = Deviance)) +
        geom_line(color = "#2166AC", linewidth = 1) +
        geom_point(color = "#B2182B", size = 2) +
        labs(
          title = "Model Convergence",
          subtitle = "Deviance over iterations",
          x = "Iteration",
          y = "Deviance"
        ) +
        theme_base +
        theme(
          panel.grid.minor = element_blank()
        )
      
      plots$convergence <- p5
    }
  }
  
  # Plot 6: Coefficient plot with confidence intervals
  if (6 %in% which && !is.null(x$summ) && "se" %in% colnames(x$summ)) {
    coef_df <- data.frame(
      Coefficient = rownames(x$summ),
      Estimate = x$summ$coef,
      SE = x$summ$se,
      stringsAsFactors = FALSE
    )
    
    # Calculate confidence intervals
    coef_df$Lower <- coef_df$Estimate - 1.96 * coef_df$SE
    coef_df$Upper <- coef_df$Estimate + 1.96 * coef_df$SE
    
    # Remove rows with NA standard errors
    coef_df <- coef_df[!is.na(coef_df$SE), ]
    
    if (nrow(coef_df) > 0) {
      # Add coefficient type for coloring
      coef_df$Type <- ifelse(grepl("^\\(Z\\)", coef_df$Coefficient), "Exogenous",
                             ifelse(grepl("^\\(alphaW\\)", coef_df$Coefficient), "Alpha",
                                   ifelse(grepl("^\\(betaW\\)", coef_df$Coefficient), "Beta", "Other")))
      
      p6 <- ggplot(coef_df, aes(x = Estimate, y = reorder(Coefficient, Estimate))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_errorbarh(
          aes(xmin = Lower, xmax = Upper, color = Type),
          height = 0.2,
          linewidth = 0.8
        ) +
        geom_point(
          aes(color = Type),
          size = 3
        ) +
        scale_color_manual(
          values = c(
            "Exogenous" = "#2166AC",
            "Alpha" = "#B2182B",
            "Beta" = "#5E4FA2",
            "Other" = "#333333"
          )
        ) +
        labs(
          title = "Coefficient Estimates",
          subtitle = "With 95% confidence intervals",
          x = "Estimate",
          y = "Coefficient",
          color = "Type"
        ) +
        theme_base +
        theme(
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9)
        )
      
      plots$coef_plot <- p6
    }
  }
  
  # Return plots
  if (length(plots) == 0) {
    cli::cli_alert_warning("No plots to display. Check that the model object contains the necessary components.")
    return(NULL)
  }
  
  # Combine plots if requested
  if (combine && length(plots) > 1) {
    # Use patchwork to combine plots
    combined <- patchwork::wrap_plots(plots, ncol = 2)
    
    if (!is.null(title)) {
      combined <- combined + 
        plot_annotation(
          title = title,
          theme = theme(plot.title = element_text(size = 16, face = "bold"))
        )
    }
    
    return(combined)
  } else if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}

#' Create a network visualization of influence matrices
#' 
#' @param x A sir object
#' @param matrix Which matrix to plot ("A" or "B")
#' @param threshold Threshold for displaying edges (absolute value)
#' @param node_labels Optional node labels
#' @param layout Layout algorithm (default "fr" for Fruchterman-Reingold)
#' @return A ggplot2 object
#' @export
plot_network.sir <- function(x, 
                             matrix = c("A", "B"),
                             threshold = 0.1,
                             node_labels = NULL,
                             layout = "fr") {
  
  matrix <- match.arg(matrix)
  
  # Check if igraph and ggraph are available
  if (!requireNamespace("igraph", quietly = TRUE) || 
      !requireNamespace("ggraph", quietly = TRUE)) {
    cli::cli_alert_warning("Network plots require 'igraph' and 'ggraph' packages")
    return(NULL)
  }
  
  # Get the appropriate matrix
  adj_matrix <- if (matrix == "A") x$A else x$B
  
  # Apply threshold
  adj_matrix[abs(adj_matrix) < threshold] <- 0
  
  # Create igraph object (use absolute values for weights since FR layout requires positive weights)
  g <- igraph::graph_from_adjacency_matrix(
    abs(adj_matrix),
    mode = "directed",
    weighted = TRUE,
    diag = FALSE
  )
  
  # Store original weights as an edge attribute for coloring
  igraph::E(g)$original_weight <- as.vector(adj_matrix[adj_matrix != 0])
  
  # Add node labels
  if (!is.null(node_labels)) {
    igraph::V(g)$name <- node_labels
  } else {
    # Default node labels
    igraph::V(g)$name <- as.character(1:igraph::vcount(g))
  }
  
  # Create network plot
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        alpha = abs(original_weight),
        color = original_weight,
        width = abs(original_weight)
      ),
      arrow = grid::arrow(length = grid::unit(2, "mm"))
    ) +
    ggraph::scale_edge_color_gradient2(
      low = "#2166AC",
      mid = "grey80",
      high = "#B2182B",
      midpoint = 0,
      name = "Influence"
    ) +
    ggraph::scale_edge_width(range = c(0.5, 2), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.3, 1), guide = "none") +
    ggraph::geom_node_point(size = 5, color = "#333333") +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      size = 3
    ) +
    ggplot2::labs(
      title = paste0("Network Visualization - Matrix ", matrix),
      subtitle = paste0("Edges with |influence| > ", threshold)
    ) +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11)
    )
  
  return(p)
}