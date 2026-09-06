data {
    int J;
    vector[J] x;
    array[J] int n, y;
}
parameters {
    real a, b;
}
model {
    y ~ binomial(n, inv_logit(a + b*x));
}
