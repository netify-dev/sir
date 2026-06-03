
test_that("get_scen_array sets the requested variable to the scenario value", {
	scen_vals = list(
		proximity = c(1, 2, 3),
		alliance = c(10, 20, 30)
	)

	arr = get_scen_array(
		"alliance",
		scen_vals = scen_vals,
		node_names = c("n1", "n2"),
		var_names = names(scen_vals),
		n_time = 2,
		value = 30
	)

	expect_equal(dim(arr), c(2, 2, 2, 2))
	off <- row(arr[, , 1, 1]) != col(arr[, , 1, 1])
	expect_true(all(arr[, , "proximity", ][off] == 2))
	expect_true(all(arr[, , "alliance", ][off] == 30))
	expect_true(all(diag(arr[, , "proximity", 1]) == 0))
	expect_true(all(diag(arr[, , "alliance", 1]) == 0))
	expect_equal(dimnames(arr)[[1]], c("n1", "n2"))
	expect_equal(dimnames(arr)[[3]], names(scen_vals))
})

test_that("get_scen_array uses the target scenario mean by default", {
	scen_vals = list(
		proximity = c(1, 2, 3),
		alliance = c(10, 20, 30)
	)

	arr = get_scen_array(
		"proximity",
		scen_vals = scen_vals,
		node_names = c("n1", "n2"),
		var_names = names(scen_vals)
	)

	off <- row(arr[, , 1]) != col(arr[, , 1])
	expect_true(all(arr[, , "proximity"][off] == 2))
	expect_true(all(arr[, , "alliance"][off] == 20))
	expect_true(all(diag(arr[, , "proximity"]) == 0))
})

test_that("get_scen_vals supports data frames, variable filters, and time slices", {
	df_vals = get_scen_vals(data.frame(a = 1:10, b = 11:20), vars = "b")
	expect_equal(names(df_vals), "b")
	expect_equal(unname(df_vals$b["50%"]), 15.5)

	arr = array(seq_len(3 * 3 * 2 * 4), dim = c(3, 3, 2, 4),
				dimnames = list(NULL, NULL, c("w1", "w2"), NULL))
	val_all = get_scen_vals(arr, vars = "w2")
	val_t2 = get_scen_vals(arr, vars = "w2", time = 2)
	expect_equal(names(val_all), "w2")
	expect_false(isTRUE(all.equal(val_all$w2, val_t2$w2)))

		val_undir = get_scen_vals(arr, vars = "w1", time = 1, directed = FALSE)
		expect_equal(length(val_undir$w1), 5)
		val_diag = get_scen_vals(arr, vars = "w1", time = 1, include_diag = TRUE)
		expect_equal(unname(val_diag$w1["10%"]), unname(stats::quantile(as.numeric(arr[, , "w1", 1]), 0.1)))
		expect_error(get_scen_vals(arr, vars = "w1", time = 1.5), "time")
	})

test_that("get_scen_array returns static 3D by default and can keep square bipartite diagonals", {
	scen_vals = list(w1 = c(1, 2, 3))
	arr = get_scen_array("w1", scen_vals, node_names = c("s1", "s2"),
						  var_names = "w1", value = 3)
	expect_equal(dim(arr), c(2, 2, 1))
	expect_true(all(diag(arr[, , "w1"]) == 0))

	arr_bp = get_scen_array("w1", scen_vals, node_names = c("s1", "s2"),
							 var_names = "w1", value = 3, zero_diag = FALSE)
	expect_true(all(diag(arr_bp[, , "w1"]) == 3))
	expect_error(
		get_scen_array("w1", scen_vals, node_names = c("s1", "s2"),
					   var_names = "w1", n_time = 1.5),
		"n_time"
	)
})
