r2 <- function(x, y) {
  ymean <- mean(y)
  ss_total <- sum((y - ymean) ^ 2)
  ss_reg <- sum((y - ymean) ^ 2)
  ss_residual <- sum((y - x) ^ 2)

  1 - ss_residual / ss_total
}


x = c(39.312, 78.304, 25.067, 8.039, 84.232, 24.578, 59.199, 87.419, 33.381, 8.919)
y = c(3.833, 66.236, 84.03, 92.198, 26.383, 31.929, 35.311, 82.228, 0.47, 74.187)


