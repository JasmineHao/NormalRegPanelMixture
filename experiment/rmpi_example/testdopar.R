library(doParallel)
cl <- makeForkCluster(20)
registerDoParallel(cl)
n<-20000
m<-100000

t <- Sys.time()
l <- foreach(i =n:m) %do% {
  sqrt(i)
}

print(Sys.time() - t)
print(sum(sapply(l, function(x) x[1])))
