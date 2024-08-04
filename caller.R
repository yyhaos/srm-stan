library(rstan)
library(rstudioapi)
# current_file_path <- rstudioapi::getActiveDocumentContext()$path
# current_dir <- dirname(current_file_path)
# setwd(current_dir)
source('_functions.R') # only use the simulation part

srm_model_file = 'srm.stan'
#####################################################################
#                            Data simulation  (copy from Reference) #
#####################################################################
#Problem structure setting
stateSet = c('A','B','C','D','E','F','G')
transSet = c('A#','AB',      
             'BA','B#','BC', 
             'CA','CB','CD', 
             'DB','DC','D#','DE', 
             'EC','E#','EF', 
             'FC','FD','FG') 
u <- length(transSet) #number of state transitions
z <- length(stateSet) - 1 #number of states minus 1
n <- 1000 #number of respondents

#Generate true value
set.seed(1)
theta <- rnorm(n)
#One needs to change the following line according to where the file is stored on your compute.
lambda <- as.matrix(read.csv('lambda19.csv')) #data simulation
dataList <- SimData0(lambda, theta, transSet)
scv <- dataList$scv
stl <- dataList$stl
#save data
# write.csv(cbind(scv, theta), 'True Value of Ability and Response Sequence.csv')


#####################################################################
#                            Data process for stan model            #
#####################################################################
# Target Data Format(from stan model):
# int<lower=1> J;                     // number of students
# int<lower=1> Z;                     // number of states
# array[Z] int outDegree;             // out degree of state
# int<lower=1> U;                     // number of transition
# array[U] int I;                     // Index value of transition
# array[Z] int transitionOffset;      // transition offset of state
# 
# int<lower=1> N;                     // number of observations
# array[N] int<lower=1, upper=J> jj;  // student for observation n
# array[N] int<lower=1, upper=Z> ss;  // start state for observation n
# array[N] int<lower=1, upper=Z> tt;  // transition for observation n
J = length(scv) # number of students
Z = length(stateSet) - 1 # number of states
outDegree = c(2,3,3,4,3,3)
U = length(transSet)
correctTransSetIndex = c(2, 5, 8, 12, 15, 18) 
I = rep(-1, U)
I[correctTransSetIndex] = 1
transitionOffset = c(1,3,6,9,13,16)
N <- sum(sapply(stl, length))
ss <- character(N)
tt <- integer(N)
jj <- integer(N)
# handle ss, ee, tt
index <- 1
for (row in seq_along(stl)) {
  for (i in stl[[row]]) {
    trans <- transSet[i]
    ss[index] <- substr(trans, 1, 1)
    tt[index] <- i
    jj[index] <- row
    index <- index + 1
  }
}
ss_indices <- match(ss, stateSet)

data_list <- list(
  J = J, # number of students
  Z = Z, # number of states
  outDegree = outDegree, # out degree of state
  U = U, # number of transition
  I = I, # Index value of transition
  transitionOffset = transitionOffset, # transition offset of state
  N = N, # number of observations
  jj = jj, # student for observation n
  ss = ss_indices, # start state for observation n
  tt = tt # transition for observation n
)


#####################################################################
#                            call stan                              #
#####################################################################

# run
fit <- stan(
  file = srm_model_file,
  data = data_list,
  iter = 30,            # iteration times
  chains = 4            # chains number by default
)


#####################################################################
#                            results                                #
#####################################################################
lambda_samples = extract(fit, pars = "lambda")$lambda
lambda_mean <- apply(lambda_samples, c(2), mean)
print(lambda_mean)
# Test Result:
# [1] -3.01044487  3.01044487 -1.49688445 -1.44155881  2.93844325 -0.17876285 -0.10367147  0.28243433  0.16122394 -0.06455934
# [11] -0.06005281 -0.03661179  0.25986249 -0.06222694 -0.19763555  0.75176154  1.24608916 -1.99785070
# True lambda:
# -3.015774971
# 3.015774971
# -1.527008364
# -1.469300759
# 2.996309123
# -0.166030616
# -0.138085059
# 0.304115675
# 0.197084306
# -0.059904982
# -0.127730483
# -0.00944884
# 0.23502588
# -0.076307083
# -0.158718798
# 0.7036586
# 1.228928174
# -1.932586774
