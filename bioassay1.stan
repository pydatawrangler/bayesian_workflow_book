data {
    int J;
    vector[J] x;
    array[J] int n, y;
}
parameters {
    real a;
    real<lower=0> b;
}
model {
    {a, b} ~ normal(0, 5);
    y ~ binomial_logit(n, a + b*x);
}
