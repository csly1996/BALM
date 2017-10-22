# BALM

An R function for implementing a Bayesian "fill in" method for correcting for publication bias. 

Publication bias occurs when the statistical significance or direction of the results between published and unpublished studies differ after controlling for study quality, which threatens the validity of the systematic review and summary of the results on a research topic. 

Conclusions based on a meta-analysis of published studies without correcting for publication bias are often optimistic and biased toward significance or positivity. We propose a Bayesian fiLl-in Meta-analysis method (BALM) for adjusting publication bias and estimating population effect size that accommodates different assumptions for publication bias.

Simulation studies were conducted to examine the performance of BALM and compare it with several commonly used/discussed and recently proposed publication bias correction methods. The simulation results suggested BALM yielded small biases, small RMSE values, and close-to-nominal-level coverage rates in inferring the population effect size and the between-study variance, and outperformed the other examined publication bias correction methods across a wide range of simulation scenarios when the publication bias mechanism is correctly specified. 

The performance of BALM was relatively sensitive to the assumed publication bias mechanism. Even with a misspecified publication bias mechanism, BALM still outperformed the naive methods without correcting for publication in inferring the overall population effect size. BALM was applied to two meta-analysis case studies to illustrate the use of BALM in real life situations. 

An R function is provided to facilitate the implementation of BALM. Guidelines on how to specify the publication bias mechanisms in BALM and how to report overall effect size estimates are provided.
