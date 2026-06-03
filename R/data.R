#' ICEWS Inter-State Material Conflict (Monthly)
#'
#' A longitudinal directed network of monthly material-conflict event counts
#' between countries, derived from the Integrated Crisis Early Warning System
#' (ICEWS) and used to illustrate the Social Influence Regression model of
#' Minhas & Hoff (2025). The data are subset to the 50 countries most involved
#' in material conflict over the sample window (February 2005 through December
#' 2012) to keep the bundled object small. See the installed provenance note
#' \code{system.file("icews-provenance.md", package = "sir")} for derivation
#' details and repository-file checksums.
#'
#' This bundled 50-country by 95-month slice is provided purely to demonstrate
#' the package on realistic data. It is \strong{illustrative only} and is
#' \emph{not} the estimation sample or the published results of Minhas & Hoff
#' (2025); fits on this slice will not reproduce the numbers reported in the
#' paper, and should not be cited as such.
#'
#' @format A named \code{list} with the array inputs \code{\link{sir}} expects:
#' \describe{
#'   \item{Y}{Integer array \code{50 x 50 x 95}. \code{Y[i, j, t]} is the number
#'     of material-conflict events initiated by country \code{i} toward country
#'     \code{j} in month \code{t}. Rows are senders, columns receivers; the
#'     diagonal is stored as zero in the bundled object and ignored by
#'     \code{\link{sir}} for one-mode fits.}
#'   \item{X}{Numeric array \code{50 x 50 x 95}, the influence-carrying network
#'     state \code{log(Y_{ij,t-1} + 1)} (the logged, lagged outcome). This is
#'     what the estimated influence matrices act on. The first bundled month
#'     uses a retained January 2005 lag from the source archive, so it cannot be
#'     reconstructed from the bundled \code{Y} alone.}
#'   \item{W}{Numeric array \code{50 x 50 x 4} of static (time-averaged)
#'     influence covariates that parameterize the sender/receiver influence
#'     matrices A and B. Slices: \code{int} (intercept anchor for the
#'     \eqn{\alpha_1 = 1} identifiability constraint), \code{ally} (alliance
#'     ties), \code{verbCoop} (logged verbal cooperation), and \code{minDistLog}
#'     (log minimum geographic distance).}
#'   \item{Z}{Numeric array \code{50 x 50 x 5 x 95} of exogenous dyadic
#'     covariates with direct effects (theta). Slices: \code{mConf} (lagged
#'     material conflict \eqn{i \to j}), \code{mConf_ji} (lagged reciprocal
#'     conflict \eqn{j \to i}), \code{minDistLog}, \code{ally}, and
#'     \code{verbCoop} (logged).}
#'   \item{countries}{Character vector of the 50 country names (the row/column
#'     labels of \code{Y}).}
#'   \item{dates}{Character vector of 95 month labels from \code{"2005-02-01"}
#'     through \code{"2012-12-01"}.}
#'   \item{metadata}{List with source DOI, replication DOI, raw archive path,
#'     access date, country-selection rule, date window, and transformation notes.}
#' }
#'
#' @details
#' This is the canonical applied example for the package. A typical fit:
#' \preformatted{
#'   data(icews)
#'   fit <- sir(icews$Y, W = icews$W, X = icews$X, Z = icews$Z,
#'              family = "poisson", seed = 1)
#'   summary(fit)
#' }
#' The covariate construction (logging verbal cooperation, the \code{int}
#' intercept slice as the influence anchor) mirrors the replication archive for
#' the published analysis.
#'
#' @source Integrated Crisis Early Warning System (ICEWS), via
#'   \code{replArchive/data/socRegData.rda} from the replication archive for
#'   Minhas, S. & Hoff, P. D. (2025), "Decomposing Network Influence: Social
#'   Influence Regression", \emph{Political Analysis}. Replication archive DOI:
#'   \doi{10.7910/DVN/VTFDX6}. See also the bundled provenance note:
#'   \code{system.file("icews-provenance.md", package = "sir")}.
#'
#' @references
#' Minhas, S. & Hoff, P. D. (2025). Decomposing Network Influence: Social
#' Influence Regression. \emph{Political Analysis}.
#' \doi{10.1017/pan.2025.10013}.
#'
#' @examples
#' data(icews)
#' dim(icews$Y)            # 50 x 50 x 95
#' head(icews$countries)
#' \dontrun{
#' fit <- sir(icews$Y, W = icews$W, X = icews$X, Z = icews$Z,
#'            family = "poisson", seed = 1)
#' summary(fit)
#' }
"icews"
