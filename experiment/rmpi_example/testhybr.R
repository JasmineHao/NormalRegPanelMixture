library(doParallel)
library("Rmpi")
options(echo=FALSE)
registerDoParallel(cores=20)
myrank<-mpi.comm.rank(0)
ranks <-mpi.comm.size(0)
n<-0
m<-0
if (myrank==0)
{
   # Master rank set number of iterations
   n<-20000
   m<-100000
}
# y = mpi.bcast (x, type, rank , comm ) , type=1 integer, type=2 double
n<-mpi.bcast(n,1,0,0);
m<-mpi.bcast(m,1,0,0);
 
# For 2 compute nodes
if (myrank==0) {
   istart<- 1
   iend  <- n/ranks
} else {
   istart<- n/ranks+1
   iend  <- n
}

ls<-0.0
lsum<-foreach(i=istart:iend,.combine='+') %dopar% {
   for(j in 1:m) {
       ls<-ls+sqrt(i+j)+cos(i+j)+sin(i+j)
   }
   ls
}
# y = mpi.allreduce(x, type, op="sum",comm), type=2 double
gsum<-mpi.allreduce(lsum,2,op="sum",0)
sprintf("Sum %f",gsum)
mpi.quit()