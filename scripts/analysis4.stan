functions {
  array[] real generatethetavector(
    real thetazero, 
    array[] real rawtheta, 
    int ntimes
  ) {
    array[ntimes] real theta;
    theta[1] = thetazero;
    for (t in 2:ntimes) {
      theta[t] = theta[(t - 1)] + rawtheta[(t - 1)];
    }
    return theta;
  }

  array[, ] real generateRmat(
    array[, ] int interventions,
    real alpha,
    array[] real gamma,
    real gammavariance,
    array[] real theta,
    real thetavariance,
    real tau,
    int ntimes, 
    int nlocations
  ) {
    array[ntimes, nlocations] real Rmat;
    for (t in 1:ntimes) { 
      for(j in 1:nlocations) { 
        Rmat[t, j] = exp(
          alpha
          + gamma[j] * sqrt(gammavariance)
          + theta[t] * sqrt(thetavariance)
          + tau * interventions[t, j]
        );
      }
    }
    return Rmat;
  }

  array[] real generateiotazeros(
    array[, ] real incidence, int nseedtimes, int nlocations
  ) {
    array[nlocations] real iota_zeros;
    for(j in 1:nlocations) {
      iota_zeros[j] = fmax(
        0.0,
        log(sum(incidence[1:nseedtimes, j]) / nseedtimes)
      );
    }
    return iota_zeros;
  }

  real calculateexpectedsusceptible() {
    return 1.0;
  }

  real calculateexpectedI_seedtimes(real iota_zero, int t, int nseedtimes) {
    return exp(iota_zero - 2 * (nseedtimes - t));
  }

  real calculatecumulativegenerationnumbers(
    array[] real predictedinfections, 
    array[] real g,
    int t,
    int glength
  ) {
    real gn;
    int startt = max(1, t - glength);
    for (x in startt:(t - 1)) {
      gn += predictedinfections[x] * g[(t - x)];
    }
    return gn;
  }

  real calculateexpectedI(
    array[] real predictedinfections, 
    array[] real g, 
    real rho, 
    real S, 
    int N,
    int t, 
    int glength
  ) {
    real z = calculatecumulativegenerationnumbers(
      predictedinfections, g, t, glength
    );
    return S * rho * z / N;
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
  real gammavariance;
  array[nlocations] real gamma;
  real thetavariance;
  real thetazero;
  array[(ntimes - 1)] real rawtheta;
  real tau;
  real meanimmuneduration;
  real psi;
  array[ntimes, (nseedtimes + nlocations)] real infectionrngs;
}

transformed parameters {
  array[ntimes] real theta = generatethetavector(thetazero, rawtheta, ntimes);
  array[ntimes, nlocations] real Rmat = generateRmat(
    interventions,
    alpha,
    gamma,
    gammavariance,
    theta,
    thetavariance,
    tau,
    ntimes, 
    nlocations
  );
  array[nlocations] real iota_zeros = generateiotazeros(
    incidence, nseedtimes, nlocations
  );

  /*array[ntimes, nlocations] real waned;
  for (t in 1:ntimes) {
    for (g in 1:nlocations) {
      waned[t, g] = exp(exponential_lpdf(t, 1 / meanimmuneduration));
    }
  }*/

  
  array[ntimes, (nseedtimes + nlocations)] real predictedinfections;
  for (j in 1:nlocations) {
    for (t in 1:nseedtimes) {
      real pred = calculateexpectedI_seedtimes(iota_zeros[j], t, nseedtimes);
      predictedinfections[t, j] = pred + infectionrngs[t, j] * pred;
    }
    for (t in (nseedtimes + 1):(nseedtimes + nlocations)) {
      real pred = calculateexpectedI(
        predictedinfections[1:(t - 1), j], 
        g, 
        Rmat[t, j], 
        1.0, 
        Ns[j],
        t, 
        glength
      );
      predictedinfections[t, j] = pred + infectionrngs[t, j] * pred;
    }
  }

}

model {
  alpha ~ normal(log(2), 1);
  gammavariance ~ exponential(1 / 0.4);
  gamma ~ normal(0, 1);
  thetavariance ~ exponential(1);
  thetazero ~ normal(0, 0.05);
  rawtheta ~ normal(0, 1);
  tau ~ normal(0, 1);
  meanimmuneduration ~ exponential(1.0 / 100.0);
  psi ~ beta(6, 4);
  for (t in 1:(nseedtimes + nlocations)) { 
    for (j in 1:nlocations) {
      infectionrngs[t, j] ~ normal(0.0, 1.0);
    }
  }
  for (t in 1:(nseedtimes + nlocations)) { 
    for (j in 1:nlocations) {
      real np = predictedinfections[t, j] * psi;
      real np_oneminus = np * (1 - psi) + 1e-10;
      incidence[t, j] ~ normal(np, np_oneminus);
    }
  }
}