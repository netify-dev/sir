# Cast Directed Dyadic Data into Array Format

Transforms long-format dyadic data (edge list) into a 3D array suitable
for SIR model analysis. Handles both dyadic and monadic variables with
proper placement in the network structure.

## Usage

``` r
cast_array(dyad_data, var, monadic = FALSE, row = FALSE)
```

## Arguments

- dyad_data:

  Data frame in long format with columns i, j, t, and the value variable
  specified by var parameter.

- var:

  Character string naming the column containing the values to be placed
  in the array.

- monadic:

  Logical indicating whether the variable is monadic (node-level) rather
  than dyadic (edge-level). Default is FALSE.

- row:

  Logical, only used when monadic=TRUE. If TRUE, uses sender attributes;
  if FALSE, uses receiver attributes. Default is FALSE.

## Value

Three-dimensional array (m x m x T) where:

- First dimension: Sender nodes (i)

- Second dimension: Receiver nodes (j)

- Third dimension: Time periods (t)

Missing edges are filled with zeros.

## Details

This function converts network data from "long" format (one row per
edge/time) to "array" format (3D array indexed by sender, receiver,
time).

**Expected Input Format:** Data frame with columns:

- i: Sender node identifier

- j: Receiver node identifier

- t: Time period

- var: The variable value for this edge

**Monadic Variables:** When monadic=TRUE, the function treats the
variable as node-level rather than edge-level. The values are placed on
the diagonal of each time slice:

- row=TRUE: Uses sender (i) attributes

- row=FALSE: Uses receiver (j) attributes

**Missing Data:** Missing edges in the input are filled with zeros in
the output array. This assumes that absence of an edge means zero
interaction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example dyadic data
dyad_data <- data.frame(
  i = c(1,1,2,2,3,3),
  j = c(2,3,1,3,1,2),
  t = c(1,1,1,1,1,1),
  trade = c(100,150,80,120,90,110)
)

# Convert to array
trade_array <- cast_array(dyad_data, "trade")
dim(trade_array)  # 3 x 3 x 1

# Monadic example (node GDP)
node_data <- data.frame(
  i = rep(1:3, each=3),
  j = rep(1:3, 3),
  t = 1,
  gdp = rep(c(1000,2000,1500), each=3)
)

gdp_array <- cast_array(node_data, "gdp", monadic=TRUE, row=TRUE)
diag(gdp_array[,,1])  # Shows node GDPs on diagonal
} # }
```
