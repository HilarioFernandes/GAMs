# GAMs

This repository contains implementations (in R) of some of the basic algorithms presented in the seminal paper on generalized additive models [1]. See also [2] and [3]. 

*Auxiliary_functions.R* contains useful functions, including a smoother (by moving averages or running lines), a function that returns the cross-validation error of a smoother, a function that calculates the optimal parameter of a smoother (i.e. one that minimizes the cross-validation error), and a function that adjusts a smoother for the logistic case, as described in the paper.

*Simulations.R* contains simulated examples of the basic techniques described in the [1]: two examples are Gaussian and the other two are logistic.


### References

[1] Trevor Hastie and Robert Tibshirani. Generalized additive models. Statistical Science, 1(3):297-318, 1986.

[2] Trevor Hastie and Robert Tibshirani. Generalized Additive Models. Chapman and Hall, 1990.

[3] Trevor J. Hastie and Robert Tibshirani. GENERALIZED ADDITIVE MODELS. 12-1984. Technical report.
