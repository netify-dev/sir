## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = FALSE
)
library(sir)
library(ggplot2)

## ----simulate-----------------------------------------------------------------
dat = sim_sir(
	m = 15, T_len = 10, p = 2, q = 1,
	family = "poisson",
	alpha = c(1, 0.4),
	beta  = c(0.5, -0.3),
	theta = -0.2,
	seed  = 42
)

