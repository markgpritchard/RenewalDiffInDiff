

functions {
  array[, ] real generateRmat(
    real alpha,
    int ntimes, 
    int nlocations
  ) {
    array[ntimes, nlocations] real Rmat;
      for (t in 1:ntimes) { 
        for (j in 1:nlocations) { 
          Rmat[t, j] = exp(
            alpha
          );
        }
      }
    return Rmat;
  }
}

data {
  int<lower=2> ntimes;
  int<lower=2> nlocations;
  int<lower=1> nseedtimes;
  int<lower=1> glength;
  array[ntimes, nlocations] int incidence;
  array[ntimes, nlocations] int interventions;
  array[nlocations] int Ns;
  array[glength] real g;
}

parameters {
  real alpha;
}

transformed parameters {
  array[ntimes, nlocations] real Rmat = generateRmat(alpha, ntimes, nlocations);
}

model {
  alpha ~ normal(0, 1);
  for (t in 1:nlocations) { 
    for (j in 1:nlocations) {
      incidence[t, j] ~ normal(Rmat[t, j], 1);
    }
  }
}

