/* Implementation Dirichlet-multinomial distribution
*/
real dirichlet_multinomial_lpmf(array[] int x, vector theta) {
    real theta0 = sum(theta);
    int n = sum(x);
    real log_den = 0.0; // the log denominator
    
    if ( n == 0 ) {
        return 0.0;
    }

    for ( i in 1:size(x) ) {
        // only the non-zero counts contribute.
        if ( x[i] > 0 ) {
            log_den += log(x[i]) + lbeta(theta[i], x[i]);
        }
    }
    return log(n) + lbeta(theta0, n) - log_den;
}


/* RNG for Dirichlet-multinomial distribution */
array[] int dirichlet_multinomial_rng(vector theta, int n) {
    // sample probability vector from dirichlet distribution
    int k = num_elements(theta);
    vector[k] p = dirichlet_rng(theta);
    // sample counts from multinomial distribution
    return multinomial_rng(p, n);
}
