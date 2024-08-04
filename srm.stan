data {
  int<lower=1> J;                     // number of students
  int<lower=1> Z;                     // number of states
  array[Z] int outDegree;             // out degree of state
  int<lower=1> U;                     // number of transition
  array[U] int I;                     // Index value of transition
  array[Z] int transitionOffset;      // transition offset of state
  
  int<lower=1> N;                     // number of observations
  array[N] int<lower=1, upper=J> jj;  // student for observation n
  array[N] int<lower=1, upper=Z> ss;  // start state for observation n
  array[N] int<lower=1, upper=U> tt;  // transition for observation n
}

parameters {
  array[J] real theta;                // ability of student
  array[U] real lambda_raw;           // unconstrained lambda
}

transformed parameters {
  array[U] real lambda;               // constrained lambda (sum of lambda is zero of same part)
  {
    int idx = 1;
    for (z in 1:Z) {
      real sum_lambda = 0;
      for (k in 1:outDegree[z]) {
        lambda[idx] = lambda_raw[idx];
        sum_lambda += lambda_raw[idx];
        idx += 1;
      }
      // Adjust the last element to ensure the sum is zero
      if (outDegree[z] > 0) {
        lambda[idx - 1] -= sum_lambda;
      }
    }
  }
}

model {
  // Prior distributions
  theta ~ normal(0, 1);
  
  // Likelihood
  for (n in 1:N) {
    int st = ss[n];
    int tr = tt[n];
    int j = jj[n];
    vector[outDegree[st]] E;          // each transition 
    int offset = transitionOffset[st];

    for (k in 1:outDegree[st]) {
      int index = k + offset - 1;
      E[k] = exp(lambda[index] + I[index] * theta[j]);
    }
    target += log(E[tr - transitionOffset[st] + 1] / sum(E));
  }
}
