test_that("enumeration has the multinomial size and reproduces the paper's exact toy distribution", {
  net <- ex_toy(); e <- ex_enumerate(net, "stacking")
  expect_equal(nrow(e$stats), 560 * 560)
  ex <- ex_exact(net, c(stacking = 0))
  expect_equal(unname(ex$mean_value), 1.75, tolerance = 1e-9)
  expect_equal(as.vector(round(ex$pmf[as.character(0:4)], 3)), c(0.129, 0.257, 0.450, 0.096, 0.043))
})

test_that("enumerated statistics agree with direct computation", {
  net <- ex_toy(); set.seed(3)
  e <- ex_enumerate(net, c("stacking", "hom_agency", "hom_club", "coloc"))
  PA <- ergm.exclusive:::.multiset_perms(sort(net$alpha)); PG <- ergm.exclusive:::.multiset_perms(sort(net$gamma))
  for (k in sample(nrow(e$stats), 25)) {
    ia <- (k - 1) %% nrow(PA) + 1; ig <- (k - 1) %/% nrow(PA) + 1
    x <- exclusive_net(PA[ia, ], PG[ig, ], net$nat, net$club_cty)
    expect_equal(unname(e$stats[k, ]), unname(ex_stats(x)))
  }
})
