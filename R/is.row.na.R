is.row.na<-function (data)
{
  n <- dim(data)
  n1 <- n[1]
  n2 <- n[2]
  ind <- rep(FALSE, n1)

  for (i in 1:n1) {
    if (sum(is.finite(data[i, ]))==n2)
      ind[i] <- TRUE
  }
  ind
}
