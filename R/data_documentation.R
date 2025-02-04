# data.R

## #' tit
## #'
## #' Des
## #'
## #' @details
## #' 
## #' @format
## #'
## #' @source
## #'
## #' @references
## #'
## #' @examples
## #'
## "data_weeks2"

## ---------------------------------------------------------------------------

#' Data from Eccleston & Russell
#'
#' @details
#' A dataframe with 3 treatment factors.
#' Each pair of factors is connected, but the 3 factors are disconnected.
#' The 'trt' column uses numbers to match Eccleston (Table 1, Design 1)
#' and letters to match Foulley (Table 13.3).
#' 
#' @source
#' Eccleston, J. and K. Russell (1975).
#' Connectedness and orthogonality in multi-factor designs.
#' Biometrika, 62, 341-345.
#' https://doi.org/10.1093/biomet/62.2.341
#' 
#' @references
#' Foulley, J. L., Bouix, J., Goffinet, B., & Elsen, J. M. (1990).
#' Connectedness in Genetic Evaluation.
#' Advanced Series in Agricultural Sciences, 277–308.
#' https://doi.org/10.1007/978-3-642-74487-7_13
#'
#' @examples
#' # Each pair of factors is connected
#' con_check(data_eccleston, ~ row + trt)
#' con_check(data_eccleston, ~ col + trt)
#' con_check(data_eccleston, ~ row + col)
#' # But all three factors are COMPLETELY disconnected
#' con_check(data_eccleston, ~ row + col + trt)
#'
#' set.seed(42)
#' data_eccleston <- transform(data_eccleston,
#'   y = rnorm(nrow(data_eccleston), mean=100))
#' con_view(data_eccleston, y ~ row*col, xlab="row", ylab="col")
#' con_view(data_eccleston, y ~ row*trt, xlab="row", ylab="trt")
#' con_view(data_eccleston, y ~ col*trt, xlab="col", ylab="trt")
#' 
"data_eccleston"


## ---------------------------------------------------------------------------

#' Data from Fernando et al.
#'
#' @details
#' A dataframe with 2 treatment factors.
#' The treatment combinations form 2 disconnected groups.
#' 
#' @source
#' Fernando et al. (1983).
#' Identifying All Connected Subsets In A Two-Way Classification Without Interaction.
#' J. of Dairy Science, 66, 1399-1402. Table 1.
#' https://doi.org/10.3168/jds.S0022-0302(83)81951-1
#' 
#' @examples
#' library(lfe)
#' cbind(data_fernando,
#'       .group=con_check(data_fernando, ~ gen + herd))
#' library(connected)
#' set.seed(42)
#' data_fernando = transform(data_fernando,
#'   y=stats::rnorm(9, mean=100))
#' con_view(data_fernando, y ~ gen*herd, cluster=FALSE,
#'   main = "Fernando unsorted")
#' con_view(data_fernando, y ~ gen*herd, main="Fernando clustered")
#'
"data_fernando"

## ---------------------------------------------------------------------------

#' Data from Searle
#'
#' @details
#' A dataframe with 2 treatment factors.
#' The treatment combinations form 3 disconnected groups.
#' 
#' @source
#' Searle (1971).
#' Linear Models. Page 324.
#' 
#' @examples
#' cbind(data_searle,
#'       .group=con_check(data_searle, ~ f1 + f2))
#' data_searle = transform(data_searle,
#'   y = rnorm(nrow(data_searle), mean=100))
#' con_view(data_searle, y ~ f1*f2, cluster=FALSE, main="Searle unsorted")
#' con_view(data_searle, y ~ f1*f2, main="Searle clustered")
#' 
"data_searle"

## ---------------------------------------------------------------------------

#' Data from Tosh
#'
#' @details
#' A dataframe with 3 treatment factors.
#' The treatment combinations form 2 disconnected groups.
#' 
#' @source
#' Tosh, J. J., and J. W. Wilton. (1990).
#' Degree of connectedness in mixed models.
#' Proceedings of the 4th World Congress on Genetics applied to Livestock Production, 480-483.
#' Page 481.
#' 
#' @examples
#' cbind(data_tosh,
#'       .group=con_check(data_tosh, ~ a + b + c))
#' data_tosh = transform(data_tosh,
#'   y = rnorm(nrow(data_tosh), mean=100))
#' library(connected)
#' con_view(data_tosh, y ~ b * c)
#' 
"data_tosh"

## ---------------------------------------------------------------------------

#' Data from Weeks & Williams example 1
#'
#' @details
#' A dataframe with 3 treatment factors.
#' The treatment combinations are connected.
#' 
#' Note: This data is based on Table 1 of Weeks & Williams.
#' Table 2 is missing treatment combination (1,2,4).
#' 
#' @source
#' Weeks, David L. & Donald R. Williams (1964).
#' A Note on the Determination of Connectedness in an N-Way Cross Classification.
#' Technometrics, 6:3, 319-324. Table 1.
#' http://dx.doi.org/10.1080/00401706.1964.10490188
#' 
#' @examples
#' library(lfe)
#' cbind(data_weeks1,
#'       .group=con_check(data_weeks1, ~ f1+f2+f3))
#'
"data_weeks1"

## ---------------------------------------------------------------------------

#' Data from Weeks & Williams example 2
#'
#' @details
#' A dataframe with 3 treatment factors.
#' The treatment combinations form 4 disconnected groups.
#' 
#' Note: This data is based on Table 3 of Weeks & Williams.
#' The groups defined in the text are missing some combinations.
#' 
#' @source
#' Weeks, David L. & Donald R. Williams (1964).
#' A Note on the Determination of Connectedness in an N-Way Cross Classification.
#' Technometrics, 6:3, 319-324. Table 3.
#' http://dx.doi.org/10.1080/00401706.1964.10490188
#' 
#' @examples
#' library(lfe)
#' cbind(data_weeks2,
#'       .group=con_check(data_weeks2, ~f1 + f2 + f3))
#'
"data_weeks2"





