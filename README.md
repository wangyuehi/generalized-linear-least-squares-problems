# generalized-linear-least-squares-problems


### Comparison of three methods for computing the GLS estimator.
* Method 1 was presented in COMP 642 lecture notes, called Paige's approach.
* Method 2 was presented in "C.C. Paige. Fast numerically stable computations for generalized linear least squares problems. SIAM J. Num. Anal., 16:165-171, 1979"
* Method 3 was a modification based on method 1, but it applied Householder transformations to Σ to get Σ ̄ = QT ΣQ. Then compute the “partial” reverse Cholesky factorization 

### References:
[1] Paige, C. C. ”Computer solution and perturbation analysis of generalized linear least squares problems.” Mathematics of Computation 33.145 (1979): 171-183.

[2] Paige, Christopher C. ”Fast numerically stable computations for general- ized linear least squares problems.” SIAM Journal on Numerical Analysis 16.1 (1979): 165-171.
