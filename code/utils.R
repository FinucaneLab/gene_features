# Not %in% function
"%ni%" <- Negate("%in%")

# t-test
DiffTTest <- function (object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE, 
          assay.type = "RNA") 
{
  data.test <- GetAssayData(object = object, assay.type = assay.type, 
                            slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  if (print.bar) {
    iterate.fxn = pblapply
  }
  else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(x = iterate.fxn(X = genes.use, FUN = function(x) {
    as.numeric(t.test(x = data.test[x, cells.1], y = data.test[x, cells.2])$statistic)
  }))
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

# T stat
my.t.test <- function(c){
  n <- sqrt(length(c))
  mean(c)*n/sd(c)
}