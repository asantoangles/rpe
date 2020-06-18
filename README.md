R code for fitting computational models of reinforcement learning using behavioural and neuroimaging data
as described in Santo-Angles et al. (under review).

**Scripts**

The following scripts contain a) functions to perform model fitting and linear regression to bold signal, and b) main loop to perform cross-validation.

- ql1.R: standard Q-learning
- ql2.R: extended Q-learning assuming alpha.chosen = alpha.nochosen
- ql3.R: extended Q-learning assuming alpha.chosen ≠ alpha.nochosen
- ac1.R: standard actor-critic
- ac2.R: extended actor-critic assuming alpha.actor.chosen = alpha.actor.nochosen
- ac3.R: extended actor-critic assuming alpha.actor.chosen ≠ alpha.actor.nochosen
