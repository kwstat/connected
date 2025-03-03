
dat <- read.csv("c:/drop/rpack/connected/data-raw/student_data.txt")
dat
head(dat)
libs(dplyr)
dat <- mutate(dat,
              test1=ifelse(test1<0, NA, test1),
              test2=ifelse(test2<0, NA, test2),
              class=ifelse(class=="m", "math", class),
              class=ifelse(class=="c", "chem", class),
              class=ifelse(class=="p", "phys", class),
              class=ifelse(class=="a", "art", class),
              class=ifelse(class=="h", "hort", class),
              class=ifelse(class=="w", "weld", class),
              student=toupper(student))
dat <- filter(dat, !(is.na(test1) & is.na(test2)))
head(dat)
tail(dat)

data_student <- dat
# Add to the connected package
write.table(data_student, file = "c:/drop/rpack/connected/data/data_student.txt",
            sep = "\t", row.names = FALSE)
promptData(data_student, filename = "c:/drop/rpack/connected/man/data_student.Rd")


## ---------------------------------------------------------------------------

if(FALSE) {
  library(lattice)
  xyplot(test2~test1, data_student)
  
  libs(reshape2)
  dat2 <- melt(dat)
  dat2 <- rename(dat2, test=variable)
  dotplot(data=dat2, value~test|class, group=student)

  
  con_check(data_student, ~ class + student)
  con_check(data_student, test1 ~ class + student)
  
  con_view(data_student, test1~student*class, main="test1", xlab="student", ylab="class")
  con_view(data_student, test2~student*class, main="test2", xlab="student", ylab="class")

  #d3 <- con_filter(data_student, ~ 7*student/class) # probably not what you want
  d3 <- con_filter(data_student, test1 ~ 7*student/class)
  con_view(d3, test2~student*class, main="test2", xlab="student", ylab="class")

  data_student |> tabyl(student,class)  
  con_filter(data_student, test1 ~ 7*student/class) |>
    tabyl(student,class)

  con_view(test2~student*class, main="test2", xlab="student", ylab="class")
  
  con_concur(data_student, ~student/ class)
  con_concur(data_student, test1~student/ class)
}
