# Network Graph Visualization of Influence Matrices

Draws the estimated influence matrix as a directed network graph, where
edges represent influence weights between nodes. Edge color and width
encode the sign and magnitude of influence. Requires the `igraph` and
`ggraph` packages to be installed.

## Usage

``` r
plot_sir_network(
  x,
  matrix = c("A", "B"),
  threshold = 0.1,
  node_labels = NULL,
  layout = "fr"
)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- matrix:

  Character string: `"A"` (default) for sender effects or `"B"` for
  receiver effects.

- threshold:

  Numeric. Edges with absolute influence below this value are hidden.
  Default is 0.1. Increase for cleaner plots with dense networks.

- node_labels:

  Optional character vector of node names. If NULL, nodes are labeled 1
  through m.

- layout:

  Character string specifying the graph layout algorithm. Default is
  `"fr"` (Fruchterman-Reingold). Other options include `"kk"`
  (Kamada-Kawai) and `"circle"`.

## Value

A `ggplot` object produced by `ggraph`. Returns NULL if `igraph` or
`ggraph` are not installed.
