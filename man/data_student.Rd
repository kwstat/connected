\name{data_student}
\alias{data_student}
\docType{data}
\title{
  Simulated Data of Student Test Scores
}
\description{
  Simulated Data of Student Test Scores
}
\usage{data("data_student")}
\format{
  A data frame with 41 observations on the following 4 variables.
  \describe{
    \item{\code{student}}{student ID}
    \item{\code{class}}{class/subject}
    \item{\code{test1}}{score on test 1}
    \item{\code{test2}}{score on test 2}
  }
}
\details{

  This simulated data imagines the following scenario:
  
  In a school, 12 students (A-M) were tested in 6 classes (math, chemistry,
  physics, art, horticulture, welding).
  
  Not all students were enrolled in all classes.
  
  In each class, students were given test 1, then later test 2.
  Test scores could be 0-100.
  
  (Using 2 tests is typical for learning assessment, measuring
  intervention effectiveness, comparison of test forms, etc.)
  
  Some students missed class on the day of a test and had no score for
  that test.
  
  Student B felt ill during art test 1 and was allowed to re-take test 1
  (both test 1 scores are included and the same test 2 score was used). 

}
\source{
  None
}
\references{
  None
}
\examples{
data(data_student)
}
\keyword{datasets}
