# data_definition.R

data_eccleston <- data.frame(
  row=factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)),
  col=factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)),
  trt=factor(c("A1","B2","E5","F6","C3","D4","G7","H8",
               "H8","F6","A1","C3","G7","E5","B2","D4"))      
)

## ---------------------------------------------------------------------------

data_fernando <- data.frame(
  gen=factor(c('G2','G1','G1','G6','G3','G4','G5','G5','G7')),
  herd=factor(c('H1','H1','H3','H3','H2','H2','H2','H4','H4'))
)

## ---------------------------------------------------------------------------

data_searle <- data.frame(
  f1 = factor(c(1,2,2,3,3,4,4,4,5,5,6,6,7,7)),
  f2 = factor(c(1,3,4,2,6,1,5,7,3,4,3,4,2,6))
)

## ---------------------------------------------------------------------------

data_tosh <- matrix(c(1,1,1, # group 1
                      1,1,6,
                      1,2,1,
                      1,2,2,
                      1,2,4,
                      2,1,1,
                      2,1,6,
                      2,2,2,
                      2,2,6,
                      1,3,3, # group 2
                      1,3,5,
                      1,4,3,
                      2,3,3,
                      2,3,5,
                      2,4,5), byrow=TRUE, ncol=3)
data_tosh <- as.data.frame(data_tosh)
colnames(data_tosh) <- c("a","b","c")
data_tosh <- transform(data_tosh,
                       a=factor(a), b=factor(b), c=factor(c))

## ---------------------------------------------------------------------------

# Note, 1,2,4 is missing in Table 2.
# This example is connected.
data_weeks1 <- matrix(c(1,1,2,
                        1,1,3,
                        1,1,4,
                        1,2,1,
                        1,2,2,
                        1,2,4,
                        1,3,1,
                        1,3,2,
                        1,3,4,
                        2,1,2,
                        2,1,4,
                        2,2,1,
                        2,2,3,
                        2,2,4,
                        2,3,1,
                        2,3,4), ncol=3, byrow=TRUE)
data_weeks1 <- as.data.frame(data_weeks1)
colnames(data_weeks1) <- c("f1","f2","f3")
set.seed(42)
data_weeks1 <- transform(data_weeks1,
                         f1=as.factor(f1),
                         f2=as.factor(f2),
                         f3=as.factor(f3)
                         )

## ---------------------------------------------------------------------------

# This is based on Weeks & Williams Table 3.
data_weeks2 <- matrix(c(2,2,2,
                   4,2,2,
                   2,2,4,
                   4,2,4,
                   2,4,2,
                   4,4,2,
                   2,4,4,
                   4,4,4,
                   2,1,1, # set 2
                   2,1,3,
                   2,1,5,
                   2,3,1,
                   2,3,3,
                   2,3,5,
                   2,5,1,
                   2,5,3,
                   2,5,5,
                   4,1,1,
                   4,1,3,
                   4,1,5,
                   4,3,1,
                   4,3,3,
                   4,3,5,
                   4,5,1,
                   4,5,3,
                   4,5,5,
                   1,2,1, # set 3
                   1,2,3,
                   1,2,5,
                   1,4,1,
                   1,4,3,
                   1,4,5,
                   3,2,1,
                   3,2,3,
                   3,2,5,
                   3,4,1,
                   3,4,3,
                   3,4,5,
                   5,2,1,
                   5,2,3,
                   5,2,5,
                   5,4,1,
                   5,4,3,
                   5,4,5,
                   1,1,2, # set 4
                   1,1,4,
                   1,3,2,
                   1,3,4,
                   1,5,2,
                   1,5,4,
                   3,1,2,
                   3,1,4,
                   3,3,2,
                   3,3,4,
                   3,5,2,
                   3,5,4,
                   5,1,2,
                   5,1,4,
                   5,3,2,
                   5,3,4,
                   5,5,2,
                   5,5,4
                   ), ncol=3, byrow=TRUE)
data_weeks2 <- as.data.frame(data_weeks2)
colnames(data_weeks2) <- c("f1","f2","f3")
data_weeks2 <- transform(data_weeks2,
                         f1=factor(f1),
                         f2=factor(f2),
                         f3=factor(f3)
                         )
