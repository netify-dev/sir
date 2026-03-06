# avoid R CMD check warnings about non-standard evaluation in ggplot2
utils::globalVariables(c("Sender", "Receiver", "Value", "density",
						 "Iteration", "Deviance", "Estimate", "Coefficient",
						 "Lower", "Upper", "Type", "Significance",
						 "original_weight", "name"))

#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_point geom_hline
#' @importFrom ggplot2 scale_fill_distiller scale_color_brewer
#' @importFrom ggplot2 theme_bw theme labs element_text element_blank element_rect
#' @importFrom ggplot2 facet_wrap coord_fixed
#' @importFrom ggplot2 geom_density geom_vline geom_line geom_errorbarh
#' @importFrom ggplot2 after_stat
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom tidyr pivot_longer
#' @importFrom stats median density reorder sd
#' @keywords internal
NULL

#' Diagnostic Plots for a Fitted SIR Model
#'
#' Produces a selection of diagnostic plots for model assessment. By default,
#' plots 1-4 (influence matrix heatmaps and distributions) are shown. Use
#' the \code{which} argument to select specific plots. All plots use
#' \code{ggplot2} and are combined via \code{patchwork} when \code{combine = TRUE}.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param which Integer vector selecting which plots to produce. Options:
#'   \describe{
#'     \item{1}{Heatmap of sender influence matrix A. Shows how each node
#'       influences others' outgoing ties.}
#'     \item{2}{Heatmap of receiver influence matrix B. Shows how each node
#'       affects others' incoming ties.}
#'     \item{3}{Histogram and density of off-diagonal A values. Useful for
#'       assessing the overall strength and distribution of sender effects.}
#'     \item{4}{Histogram and density of off-diagonal B values. Same for
#'       receiver effects.}
#'     \item{5}{Convergence trace plot showing deviance across ALS iterations.
#'       Iteration history is always stored in the fitted model.}
#'     \item{6}{Coefficient plot with 95\% confidence intervals. Requires
#'       standard errors (\code{calc_se = TRUE}). Parameters are grouped by
#'       type (exogenous, alpha, beta).}
#'   }
#' @param combine Logical. If TRUE (default), combines selected plots into
#'   a single patchwork layout. If FALSE, returns a list of individual plots.
#' @param title Optional character string for the combined plot title.
#' @param theme_base A ggplot2 theme applied to all plots. Default is
#'   \code{theme_bw()}.
#' @param ... Additional arguments (unused).
#' @return When \code{combine = TRUE} and multiple plots are requested, a
#'   \code{patchwork} object. When a single plot is requested, a \code{ggplot}
#'   object. When \code{combine = FALSE}, a named list of \code{ggplot} objects.
#'
#' @examples
#' \dontrun{
#' model <- sir(Y, W, X, family = "poisson")
#'
#' # default: influence heatmaps and distributions
#' plot(model)
#'
#' # all plots combined with a title
#' plot(model, which = 1:6, title = "SIR Diagnostics")
#'
#' # individual plots for custom arrangement
#' plots <- plot(model, which = c(1, 6), combine = FALSE)
#' plots$A_heatmap
#' plots$coef_plot
#' }
#' @export
plot.sir <- function(x, 
					 which = 1:4,
					 combine = TRUE,
					 title = NULL,
					 theme_base = theme_bw(),
					 ...) {
  
  plots <- list()
  
  # plot 1: influence matrix A heatmap
  if (1 %in% which && !is.null(x$A)) {
	# for dynamic W (3D A), average across time
	A_mat <- if (length(dim(x$A)) == 3) apply(x$A, c(1, 2), mean) else x$A
	A_df <- expand.grid(
	  Sender = factor(1:nrow(A_mat)),
	  Receiver = factor(1:ncol(A_mat))
	)
	A_df$Value <- as.vector(A_mat)
	A_df$Value[A_df$Sender == A_df$Receiver] <- NA

	p1 <- ggplot(A_df, aes(x = Sender, y = Receiver, fill = Value))
	p1 <- p1 + geom_tile()
	p1 <- p1 + scale_fill_distiller(palette = "RdBu", direction = -1,
									 name = "Influence", na.value = "grey90")
	p1 <- p1 + coord_fixed()
	p1 <- p1 + labs(title = "Influence Matrix A (Sender Effects)",
					 x = "Sender Node", y = "Receiver Node")
	p1 <- p1 + theme_base
	p1 <- p1 + theme(
	  panel.border = element_blank(),
	  axis.text = element_text(size = 8),
	  axis.text.x = element_text(angle = 45, hjust = 1),
	  axis.ticks = element_blank(),
	  legend.position = "top"
	)

	plots$A_heatmap <- p1
  }
  
  # plot 2: influence matrix B heatmap
  if (2 %in% which && !is.null(x$B)) {
	B_mat <- if (length(dim(x$B)) == 3) apply(x$B, c(1, 2), mean) else x$B
	B_df <- expand.grid(
	  Sender = factor(1:nrow(B_mat)),
	  Receiver = factor(1:ncol(B_mat))
	)
	B_df$Value <- as.vector(B_mat)
	B_df$Value[B_df$Sender == B_df$Receiver] <- NA

	p2 <- ggplot(B_df, aes(x = Sender, y = Receiver, fill = Value))
	p2 <- p2 + geom_tile()
	p2 <- p2 + scale_fill_distiller(palette = "PuOr", direction = -1,
									 name = "Influence", na.value = "grey90")
	p2 <- p2 + coord_fixed()
	p2 <- p2 + labs(title = "Influence Matrix B (Receiver Effects)",
					 x = "Sender Node", y = "Receiver Node")
	p2 <- p2 + theme_base
	p2 <- p2 + theme(
	  panel.border = element_blank(),
	  axis.text = element_text(size = 8),
	  axis.text.x = element_text(angle = 45, hjust = 1),
	  axis.ticks = element_blank(),
	  legend.position = "top"
	)

	plots$B_heatmap <- p2
  }
  
  # plot 3: distribution of A matrix (off-diagonal)
  if (3 %in% which && !is.null(x$A)) {
	A_use <- if (length(dim(x$A)) == 3) apply(x$A, c(1, 2), mean) else x$A
	A_offdiag <- A_use[row(A_use) != col(A_use)]
	A_dist_df <- data.frame(Value = A_offdiag)

	p3 <- ggplot(A_dist_df, aes(x = Value))
	p3 <- p3 + geom_histogram(aes(y = after_stat(density)), bins = 30)
	p3 <- p3 + geom_density(linewidth = 1)
	p3 <- p3 + geom_vline(aes(xintercept = mean(Value)),
						   linetype = "dashed", linewidth = 0.8)
	p3 <- p3 + geom_vline(aes(xintercept = median(Value)),
						   linetype = "dotted", linewidth = 0.8)
	p3 <- p3 + labs(title = "Distribution of A Matrix Values",
					 subtitle = paste0("Mean: ", round(mean(A_offdiag), 4),
									  " | Median: ", round(median(A_offdiag), 4)),
					 x = "Influence Value", y = "Density")
	p3 <- p3 + theme_base
	p3 <- p3 + theme(panel.border = element_blank(),
					  panel.grid.minor = element_blank(),
					  axis.ticks = element_blank(),
					  legend.position = "top")

	plots$A_dist <- p3
  }
  
  # plot 4: distribution of B matrix (off-diagonal)
  if (4 %in% which && !is.null(x$B)) {
	B_use <- if (length(dim(x$B)) == 3) apply(x$B, c(1, 2), mean) else x$B
	B_offdiag <- B_use[row(B_use) != col(B_use)]
	B_dist_df <- data.frame(Value = B_offdiag)

	p4 <- ggplot(B_dist_df, aes(x = Value))
	p4 <- p4 + geom_histogram(aes(y = after_stat(density)), bins = 30)
	p4 <- p4 + geom_density(linewidth = 1)
	p4 <- p4 + geom_vline(aes(xintercept = mean(Value)),
						   linetype = "dashed", linewidth = 0.8)
	p4 <- p4 + geom_vline(aes(xintercept = median(Value)),
						   linetype = "dotted", linewidth = 0.8)
	p4 <- p4 + labs(title = "Distribution of B Matrix Values",
					 subtitle = paste0("Mean: ", round(mean(B_offdiag), 4),
									  " | Median: ", round(median(B_offdiag), 4)),
					 x = "Influence Value", y = "Density")
	p4 <- p4 + theme_base
	p4 <- p4 + theme(panel.border = element_blank(),
					  panel.grid.minor = element_blank(),
					  axis.ticks = element_blank(),
					  legend.position = "top")

	plots$B_dist <- p4
  }
  
  # plot 5: deviance/convergence history
  if (5 %in% which && !is.null(x$history) && !is.null(x$history$DEV)) {
	dev_history <- x$history$DEV
	if (nrow(dev_history) > 1) {
	  dev_df <- data.frame(
		Iteration = 1:nrow(dev_history),
		Deviance = dev_history[, 2]
	  )

	  p5 <- ggplot(dev_df, aes(x = Iteration, y = Deviance))
	  p5 <- p5 + geom_line(linewidth = 1)
	  p5 <- p5 + geom_point(size = 2)
	  p5 <- p5 + labs(title = "Model Convergence",
					   subtitle = "Deviance over iterations",
					   x = "Iteration", y = "Deviance")
	  p5 <- p5 + theme_base
	  p5 <- p5 + theme(panel.border = element_blank(),
						panel.grid.minor = element_blank(),
						axis.ticks = element_blank(),
						legend.position = "top")

	  plots$convergence <- p5
	}
  }
  
  # plot 6: coefficient plot with significance-colored confidence intervals
  if (6 %in% which && !is.null(x$summ) && "se" %in% colnames(x$summ)) {
	coef_df <- data.frame(
	  Coefficient = rownames(x$summ),
	  Estimate = x$summ$coef,
	  SE = x$summ$se,
	  stringsAsFactors = FALSE
	)
	coef_df$Lower <- coef_df$Estimate - 1.96 * coef_df$SE
	coef_df$Upper <- coef_df$Estimate + 1.96 * coef_df$SE
	coef_df <- coef_df[!is.na(coef_df$SE), ]

	if (nrow(coef_df) > 0) {
	  # compute p-values and significance bands
	  coef_df$z <- coef_df$Estimate / coef_df$SE
	  coef_df$pval <- 2 * (1 - pnorm(abs(coef_df$z)))
	  coef_df$Significance <- cut(coef_df$pval,
								  breaks = c(0, 0.01, 0.05, 0.10, 1),
								  labels = c("p < 0.01", "p < 0.05", "p < 0.10", "n.s."),
								  include.lowest = TRUE, right = TRUE)

	  p6 <- ggplot(coef_df, aes(x = Estimate, y = reorder(Coefficient, Estimate)))
	  p6 <- p6 + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
	  p6 <- p6 + geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Significance),
								 height = 0.2, linewidth = 0.8)
	  p6 <- p6 + geom_point(aes(color = Significance), size = 3)
	  p6 <- p6 + ggplot2::scale_color_brewer(palette = "Dark2", drop = FALSE)
	  p6 <- p6 + labs(title = "Coefficient Estimates",
					   subtitle = "With 95% confidence intervals (colored by significance)",
					   x = "Estimate", y = "Coefficient", color = "Significance")
	  p6 <- p6 + theme_base
	  p6 <- p6 + theme(panel.border = element_blank(),
						panel.grid.major.y = element_blank(),
						axis.text.y = element_text(size = 9),
						axis.ticks = element_blank(),
						legend.position = "top")

	  plots$coef_plot <- p6
	}
  }
  
  # return plots
  if (length(plots) == 0) {
	cli::cli_alert_warning("No plots to display. Check that the model object contains the necessary components.")
	return(NULL)
  }
  
  # combine plots if requested
  if (combine && length(plots) > 1) {
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

#' Network Graph Visualization of Influence Matrices
#'
#' Draws the estimated influence matrix as a directed network graph, where
#' edges represent influence weights between nodes. Edge color and width
#' encode the sign and magnitude of influence. Requires the \code{igraph}
#' and \code{ggraph} packages to be installed.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param matrix Character string: \code{"A"} (default) for sender effects
#'   or \code{"B"} for receiver effects.
#' @param threshold Numeric. Edges with absolute influence below this value
#'   are hidden. Default is 0.1. Increase for cleaner plots with dense
#'   networks.
#' @param node_labels Optional character vector of node names. If NULL,
#'   nodes are labeled 1 through m.
#' @param layout Character string specifying the graph layout algorithm.
#'   Default is \code{"fr"} (Fruchterman-Reingold). Other options include
#'   \code{"kk"} (Kamada-Kawai) and \code{"circle"}.
#' @return A \code{ggplot} object produced by \code{ggraph}. Returns NULL
#'   if \code{igraph} or \code{ggraph} are not installed.
#' @export
plot_sir_network <- function(x,
							 matrix = c("A", "B"),
							 threshold = 0.1,
							 node_labels = NULL,
							 layout = "fr") {
  
  matrix <- match.arg(matrix)
  
  # check if igraph and ggraph are available
  if (!requireNamespace("igraph", quietly = TRUE) || 
	  !requireNamespace("ggraph", quietly = TRUE)) {
	cli::cli_alert_warning("Network plots require 'igraph' and 'ggraph' packages")
	return(NULL)
  }
  
  # get the appropriate matrix (average across time for dynamic W)
  adj_raw <- if (matrix == "A") x$A else x$B
  adj_matrix <- if (length(dim(adj_raw)) == 3) apply(adj_raw, c(1, 2), mean) else adj_raw
  
  # apply threshold
  adj_matrix[abs(adj_matrix) < threshold] <- 0
  
  # create igraph object (absolute values for FR layout)
  g <- igraph::graph_from_adjacency_matrix(
	abs(adj_matrix),
	mode = "directed",
	weighted = TRUE,
	diag = FALSE
  )
  
  # store original weights for coloring
  igraph::E(g)$original_weight <- as.vector(adj_matrix[adj_matrix != 0])
  
  # add node labels
  if (!is.null(node_labels)) {
	igraph::V(g)$name <- node_labels
  } else {
	# default node labels
	igraph::V(g)$name <- as.character(1:igraph::vcount(g))
  }
  
  # create network plot
  p <- ggraph::ggraph(g, layout = layout) +
	ggraph::geom_edge_link(
	  ggplot2::aes(
		alpha = abs(original_weight),
		color = original_weight,
		width = abs(original_weight)
	  ),
	  arrow = grid::arrow(length = grid::unit(2, "mm"))
	) +
	ggraph::scale_edge_color_distiller(
	  palette = "RdBu",
	  direction = -1,
	  name = "Influence"
	) +
	ggraph::scale_edge_width(range = c(0.5, 2), guide = "none") +
	ggraph::scale_edge_alpha(range = c(0.3, 1), guide = "none") +
	ggraph::geom_node_point(size = 5) +
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