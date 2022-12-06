R code for fitting computational models of reinforcement learning using behavioural and neuroimaging data
as described in Santo-Angles A, Fuentes-Claramonte P, Argila-Plaza I, Guardiola-Ripoll M, Almodóvar-Payá C, Munuera J, McKenna PJ, Pomarol-Clotet E, Radua J. Reward and fictive prediction error signals in ventral striatum: asymmetry between factual and counterfactual processing. Brain Struct Funct. 2021 Jun;226(5):1553-1569. doi: 10.1007/s00429-021-02270-3.

The following scripts contain a) functions to perform model fitting and linear regression to bold signal, and b) the main loop to perform cross-validation.

- ql1.R: standard Q-learning
- ql2.R: extended Q-learning assuming alpha.chosen = alpha.nochosen
- ql3.R: extended Q-learning assuming alpha.chosen ≠ alpha.nochosen
- ac1.R: standard actor-critic
- ac2.R: extended actor-critic assuming alpha.actor.chosen = alpha.actor.nochosen
- ac3.R: extended actor-critic assuming alpha.actor.chosen ≠ alpha.actor.nochosen
