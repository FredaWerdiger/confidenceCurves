# test
df <- confidenceCurves::makeConfidenceCurves(
  theta.estimator = -0.22,
  confidence.lower = -0.36,
  confidence.upper = -0.07
  )

testthat::expect_equal(df$conf.benefit,
                       0.999,
                       tolerance = 0.001)
