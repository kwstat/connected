
test_that("con_check works", {
  expect_equal(
    con_check(data_eccleston, ~ row + trt),
    structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                1L, 1L, 1L), levels = "1", class = "factor")
  )

  expect_equal(
    con_check(data_eccleston, ~ row + col + trt),
    structure(c(16L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L),
              levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"),
              class = "factor")
  )
})

test_that("con_view", {
  set.seed(42)
  data_fernando = transform(data_fernando,
                            y=stats::rnorm(9, mean=100))
  expect_warning(con_view(data_fernando, y ~ gen*herd, cluster=FALSE,
           main = "Fernando unsorted"))
})

test_that("con_filter", {
  tab <- data.frame(gen=c("G1","G1","G1","G1", "G2","G2","G2", "G3"),
                    state=c("S1","S2","S3","S4", "S1","S2","S4", "S1"))
  expect_warning( con_filter(tab, ~ 2 * state / gen) )
  })
