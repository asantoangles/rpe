R code for fitting computational models of reinforcement learning using behavioural and neuroimaging data
(combined fitting) as described in Santo-Angles et al. (under review).

**Scripts**

The following scripts contain functions to perform model fitting and linear regression to bold signal, and the main loop to perform cross-validation.

- ql1.R: standard Q-learning
- ql2.R: extended Q-learning assuming alpha.chosen = alpha.nochosen
- ql3.R: extended Q-learning assuming alpha.chosen ≠ alpha.nochosen
- ac1.R: standard actor-critic
- ac2.R: extended actor-critic assuming alpha.actor.chosen = alpha.actor.nochosen
- ac3.R: extended actor-critic assuming alpha.actor.chosen ≠ alpha.actor.nochosen

**Data**

As an example, behavioral and neuroimaging data from one subject, e.g. timecourse of one region-of-interest, is available. To fit standard Q-learning model and perform linear regression with BOLD signal, change paths in faked_feat.fsf (lines 33, 270 and 307) and run lines 1-222 of ql1.R
