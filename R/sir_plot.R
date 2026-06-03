# avoid R CMD check warnings about non-standard evaluation in ggplot2
utils::globalVariables(c("Sender", "Receiver", "Influenced", "Source",
						 "Value", "density",
						 "Iteration", "Deviance", "Estimate", "Coefficient",
						 "Lower", "Upper", "Type", "Significance",
						 "original_weight", "name"))

#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_point geom_hline
#' @importFrom ggplot2 theme_bw theme labs element_text element_blank element_rect
#' @importFrom ggplot2 facet_wrap coord_fixed
#' @importFrom ggplot2 geom_density geom_vline geom_line geom_errorbar
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
#' @param node_labels Optional character vector of node names used for the
#'   heatmap row/column tick labels (plots 1-2). If NULL (default), labels are
#'   taken from the \code{dimnames} of the influence matrix \code{x$A} when
#'   present, otherwise \code{1:m}.
#' @param period Optional integer selecting a single time slice for the
#'   heatmaps when the fit uses dynamic influence covariates (\code{x$A} is
#'   \code{m x m x T}). If NULL (default), dynamic heatmaps show the
#'   time-average across all periods (clearly labeled as such).
#' @param theme_base A ggplot2 theme applied to all plots. Default is
#'   \code{theme_bw()}.
#' @param ... Additional arguments (unused).
#' @return When \code{combine = TRUE} and multiple plots are requested, a
#'   \code{patchwork} object. When a single plot is requested, a \code{ggplot}
#'   object. When \code{combine = FALSE}, a named list of \code{ggplot} objects.
#'
#' @examples
#' \donttest{
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
#'            calc_se = TRUE, seed = 1)
#' plot(fit)                                  # influence heatmaps + distributions
#' plot(fit, which = 1:6, title = "SIR Diagnostics")
#' plots <- plot(fit, which = c(1, 6), combine = FALSE)
#' plots$A_heatmap
#' }
#' @export
plot.sir_fit <- function(x,
					 which = 1:4,
					 combine = TRUE,
					 title = NULL,
					 node_labels = NULL,
					 period = NULL,
					 theme_base = theme_bw(),
					 ...) {

	plots <- list()

	# collapse a dynamic influence matrix to one display matrix
	slice_influence <- function(M) {
	if (length(dim(M)) == 3) {
	  Tn <- dim(M)[3]
	  if (!is.null(period)) {
		if (period < 1 || period > Tn) {
		  cli::cli_abort("`period` must be between 1 and {Tn} (the number of time slices).")
		}
		list(mat = M[, , period], sub = paste0("Time slice t = ", period, " of ", Tn))
	  } else {
		list(mat = apply(M, c(1, 2), mean),
			 sub = paste0("Time-average over ", Tn, " periods (sign-flipping influence may cancel)"))
	  }
	} else {
	  list(mat = M, sub = NULL)
	}
	}

	# resolve heatmap tick labels
	resolve_labels <- function(mat) {
	m <- nrow(mat)
	if (!is.null(node_labels)) {
	  if (length(node_labels) != m) {
		cli::cli_abort("`node_labels` has length {length(node_labels)} but the influence matrix has {m} nodes.")
	  }
	  return(as.character(node_labels))
	}
	dn <- dimnames(mat)
	if (!is.null(dn) && !is.null(dn[[1]])) return(as.character(dn[[1]]))
	as.character(seq_len(m))
	}

	# build a centered diverging fill scale
	sym_fill <- function(values, low, high) {
	rng <- range(values, na.rm = TRUE, finite = TRUE)
	lim <- max(abs(rng))
	if (!is.finite(lim) || lim == 0) lim <- 1  # guard all-equal / degenerate
	ggplot2::scale_fill_gradient2(
	  low = low, mid = "#F7F7F7", high = high, midpoint = 0,
	  limits = c(-lim, lim), name = "Influence", na.value = "grey90"
	)
	}

	# plot 1: influence matrix A heatmap
	if (1 %in% which && !is.null(x$A)) {
	# prepare sender-side influence heatmap data
	A_sl <- slice_influence(x$A)
	A_mat <- A_sl$mat
	A_labs <- resolve_labels(A_mat)
	# rows are influenced nodes and columns are source nodes
	A_df <- expand.grid(
	  Influenced = factor(A_labs, levels = A_labs),
	  Source = factor(A_labs, levels = A_labs)
	)
	A_df$Value <- as.vector(A_mat)
	A_df$Value[A_df$Influenced == A_df$Source] <- NA

	# subtitle: keep the index legend, append the time note when dynamic
	A_sub <- "A[i,k]: source k -> influenced sender i"
	if (!is.null(A_sl$sub)) A_sub <- paste0(A_sub, "\n", A_sl$sub)

	p1 <- ggplot(A_df, aes(x = Source, y = Influenced, fill = Value))
	p1 <- p1 + geom_tile()
	if (nrow(A_mat) > 25) {
		keep <- A_labs[unique(round(seq(1, nrow(A_mat), length.out = 12)))]
		p1 <- p1 + ggplot2::scale_x_discrete(breaks = keep)
		p1 <- p1 + ggplot2::scale_y_discrete(breaks = keep)
	}
	p1 <- p1 + sym_fill(A_df$Value, low = "#2166AC", high = "#B2182B")
	p1 <- p1 + coord_fixed()
	p1 <- p1 + labs(title = "Sender Influence A",
					 subtitle = A_sub,
					 x = "Source node (k)", y = "Influenced sender (i)")
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
	B_sl <- slice_influence(x$B)
	B_mat <- B_sl$mat
	B_labs <- resolve_labels(B_mat)
	# B[j, l] = how much source receiver l (column) shapes influenced receiver j's
		# (row) incoming ties: in eta = A X B', (X B')_{.,j} = sum_l x_{.,l} B[j,l].
	B_df <- expand.grid(
	  Influenced = factor(B_labs, levels = B_labs),
	  Source = factor(B_labs, levels = B_labs)
	)
	B_df$Value <- as.vector(B_mat)
	B_df$Value[B_df$Influenced == B_df$Source] <- NA

	B_sub <- "B[j,l]: source receiver l -> influenced receiver j"
	if (!is.null(B_sl$sub)) B_sub <- paste0(B_sub, "\n", B_sl$sub)

	p2 <- ggplot(B_df, aes(x = Source, y = Influenced, fill = Value))
	p2 <- p2 + geom_tile()
	if (nrow(B_mat) > 25) {
		keep <- B_labs[unique(round(seq(1, nrow(B_mat), length.out = 12)))]
		p2 <- p2 + ggplot2::scale_x_discrete(breaks = keep)
		p2 <- p2 + ggplot2::scale_y_discrete(breaks = keep)
	}
	p2 <- p2 + sym_fill(B_df$Value, low = "#542788", high = "#E08214")
	p2 <- p2 + coord_fixed()
	p2 <- p2 + labs(title = "Receiver Influence B",
					 subtitle = B_sub,
					 x = "Source node (l)", y = "Influenced receiver (j)")
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
	A_sl3 <- slice_influence(x$A)
	A_use <- A_sl3$mat
	A_offdiag <- A_use[row(A_use) != col(A_use)]
	A_dist_df <- data.frame(Value = A_offdiag)

	A_sub3 <- paste0("Mean: ", round(mean(A_offdiag), 4),
					 " | Median: ", round(median(A_offdiag), 4))
	if (!is.null(A_sl3$sub)) A_sub3 <- paste0(A_sub3, "\n", A_sl3$sub)

	p3 <- ggplot(A_dist_df, aes(x = Value))
	p3 <- p3 + geom_histogram(aes(y = after_stat(density)), bins = 30)
	p3 <- p3 + geom_density(linewidth = 1)
	p3 <- p3 + geom_vline(aes(xintercept = mean(Value)),
						   linetype = "dashed", linewidth = 0.8)
	p3 <- p3 + geom_vline(aes(xintercept = median(Value)),
						   linetype = "dotted", linewidth = 0.8)
	p3 <- p3 + labs(title = "Distribution of A Matrix Values",
					 subtitle = A_sub3,
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
	B_sl4 <- slice_influence(x$B)
	B_use <- B_sl4$mat
	B_offdiag <- B_use[row(B_use) != col(B_use)]
	B_dist_df <- data.frame(Value = B_offdiag)

	B_sub4 <- paste0("Mean: ", round(mean(B_offdiag), 4),
					 " | Median: ", round(median(B_offdiag), 4))
	if (!is.null(B_sl4$sub)) B_sub4 <- paste0(B_sub4, "\n", B_sl4$sub)

	p4 <- ggplot(B_dist_df, aes(x = Value))
	p4 <- p4 + geom_histogram(aes(y = after_stat(density)), bins = 30)
	p4 <- p4 + geom_density(linewidth = 1)
	p4 <- p4 + geom_vline(aes(xintercept = mean(Value)),
						   linetype = "dashed", linewidth = 0.8)
	p4 <- p4 + geom_vline(aes(xintercept = median(Value)),
						   linetype = "dotted", linewidth = 0.8)
	p4 <- p4 + labs(title = "Distribution of B Matrix Values",
					 subtitle = B_sub4,
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

	  # early warm-up iterations can blow up (deviance ~1e26), which flattens
	  # the converged tail to a flat line on a linear axis. use a log10 y-axis
	  # so the informative tail stays readable. fall back to linear (and say so)
	  # only when some deviance is non-positive/non-finite (cannot be logged).
	  use_log <- all(is.finite(dev_df$Deviance)) && all(dev_df$Deviance > 0)
	  dev_sub <- if (use_log) "Deviance over iterations (log10 y-axis)" else
		"Deviance over iterations"

	  p5 <- ggplot(dev_df, aes(x = Iteration, y = Deviance))
	  p5 <- p5 + geom_line(linewidth = 1)
	  p5 <- p5 + geom_point(size = 2)
	  if (use_log) {
		p5 <- p5 + ggplot2::scale_y_log10()
	  }
	  p5 <- p5 + labs(title = "Model Convergence",
					   subtitle = dev_sub,
					   x = "Iteration",
					   y = if (use_log) "Deviance (log10)" else "Deviance")
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
	  # colour encodes direction x significance: blue significant positive,
	  # red significant negative, grey not significant (classical Wald p-value)
	  coef_df$z <- coef_df$Estimate / coef_df$SE
	  coef_df$pval <- 2 * (1 - pnorm(abs(coef_df$z)))
	  sig <- coef_df$pval < 0.05
	  coef_df$Significance <- factor(
		  ifelse(!sig, "Not significant",
			  ifelse(coef_df$Estimate > 0, "Significant positive", "Significant negative")),
		  levels = c("Significant positive", "Significant negative", "Not significant"))

	  p6 <- ggplot(coef_df, aes(x = Estimate, y = reorder(Coefficient, Estimate)))
	  p6 <- p6 + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
		  p6 <- p6 + ggplot2::geom_errorbar(
			  aes(xmin = Lower, xmax = Upper, color = Significance),
			  orientation = "y", width = 0.2, linewidth = 0.8)
	  p6 <- p6 + geom_point(aes(color = Significance), size = 3)
	  p6 <- p6 + ggplot2::scale_color_manual(
		  values = c("Significant positive" = "#2166AC",
					 "Significant negative" = "#B2182B",
					 "Not significant"      = "grey60"),
		  drop = FALSE)
	  p6 <- p6 + labs(title = "Coefficient Estimates",
					   subtitle = "95% Wald intervals from classical SEs; cross-check with boot_sir()",
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
#' @param node_labels Optional character vector of node names. If NULL, labels
#'   are taken from the influence-matrix dimnames when present, otherwise nodes
#'   are labeled 1 through m.
#' @param layout Character string specifying the graph layout algorithm.
#'   Default is \code{"fr"} (Fruchterman-Reingold). Other options include
#'   \code{"kk"} (Kamada-Kawai) and \code{"circle"}.
#' @return A \code{ggplot} object produced by \code{ggraph}. Arrows point from
#'   the source/influencer node to the influenced node. Returns NULL
#'   if \code{igraph} or \code{ggraph} are not installed.
#' @examples
#' \donttest{
#' if (requireNamespace("igraph", quietly = TRUE) &&
#'     requireNamespace("ggraph", quietly = TRUE)) {
#'   dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#'   fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#'   plot_sir_network(fit, matrix = "A", threshold = 0.1)
#' }
#' }
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
		if (is.null(node_labels)) {
			dn <- dimnames(adj_matrix)
			if (!is.null(dn) && !is.null(dn[[1]])) {
				node_labels <- as.character(dn[[1]])
			}
		}
	
	# apply threshold
	adj_matrix[abs(adj_matrix) < threshold] <- 0
	
	# orient rows as sources for igraph
	graph_matrix <- t(adj_matrix)

	# build graph with absolute weights for layout
	g <- igraph::graph_from_adjacency_matrix(
		abs(graph_matrix),
		mode = "directed",
		weighted = TRUE,
		diag = FALSE
	)
	
	# keep signed weights for edge color
	edge_ends <- igraph::ends(g, igraph::E(g), names = FALSE)
	igraph::E(g)$original_weight <- graph_matrix[edge_ends]

	# symmetric limit so the diverging edge color is anchored at 0: equal-magnitude
	# +/- influence get equally-saturated opposite colors, 0 maps to neutral.
	edge_lim <- max(abs(igraph::E(g)$original_weight))
	if (!is.finite(edge_lim) || edge_lim == 0) edge_lim <- 1

	# add node labels
		if (!is.null(node_labels)) {
		igraph::V(g)$name <- as.character(node_labels)
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
	ggraph::scale_edge_color_gradient2(
	  low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
	  limits = c(-edge_lim, edge_lim),
	  name = "Influence",
	  guide = ggraph::guide_edge_colourbar()
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
	ggplot2::theme_bw() +
	ggplot2::theme(
	  panel.border = ggplot2::element_blank(),
	  axis.ticks = ggplot2::element_blank(),
	  legend.position = "top",
	  plot.title = ggplot2::element_text(size = 14, face = "bold"),
	  plot.subtitle = ggplot2::element_text(size = 11)
	)
	
	return(p)
}
