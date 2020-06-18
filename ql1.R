library(metafor)
library(Matrix)

# call_feat() function convolves model variables to hemodynamic response function 
# using FSL FEAT to get model regresors
call_feat = function (time, duration, value) {
  write.table(data.frame(time, duration, value), "faked_feat/ql1/rpe.txt",
              col.names = FALSE, row.names = FALSE)
  system("rm -rf faked_feat/ql1/result*.feat")
  system("feat faked_feat/ql1/faked_feat.fsf", ignore.stdout = TRUE)
  X = read.table("faked_feat/ql1/result.feat/design.mat", skip = 5)
  colnames(X) = c("x", "t")
  X
}

# feat() function performs linear regression of bold signal (as dependent variable) 
# and model regresor (as independent variable)
feat = function (subject, behav, bold) {
  data = behav$data
  rois = colnames(bold)
  V = call_feat(data$time_resp, data$time_feed - data$time_resp, data$V)
  rpe = call_feat(data$time_feed, 1, data$rpe)
  abs_rpe = call_feat(data$time_feed, 1, abs(data$rpe))
  pos_rpe = call_feat(data$time_feed[which(data$rpe > 0)], 1, data$rpe[which(data$rpe > 0)])
  neg_rpe = call_feat(data$time_feed[which(data$rpe < 0)], 1, -data$rpe[which(data$rpe < 0)])
  if (!is.null(data$rpe_nc)) {
    rpe_nc = call_feat(data$time_feed, 1, data$rpe_nc)
    abs_rpe_nc = call_feat(data$time_feed, 1, abs(data$rpe_nc))
    pos_rpe_nc = call_feat(data$time_feed[which(data$rpe_nc > 0)], 1, data$rpe_nc[which(data$rpe_nc > 0)])
    neg_rpe_nc = call_feat(data$time_feed[which(data$rpe_nc < 0)], 1, -data$rpe_nc[which(data$rpe_nc < 0)])
  }
  y = data.frame()
  for (roi in rois) {
    m = lm(bold[[roi]] ~ V$x + V$t)
    y = rbind(y, data.frame(
      subject, roi, what = "V",
      y = coef(m)[2], vy = vcov(m)[2, 2]
    ))
    m = lm(bold[[roi]] ~ rpe$x + rpe$t)
    y = rbind(y, data.frame(
      subject, roi, what = "rpe",
      y = coef(m)[2], vy = vcov(m)[2, 2]
    ))
    m = lm(bold[[roi]] ~ abs_rpe$x + abs_rpe$t)
    y = rbind(y, data.frame(
      subject, roi, what = "abs_rpe",
      y = coef(m)[2], vy = vcov(m)[2, 2]
    ))
    m = lm(bold[[roi]] ~ pos_rpe$x + pos_rpe$t)
    y = rbind(y, data.frame(
      subject, roi, what = "pos_rpe",
      y = coef(m)[2], vy = vcov(m)[2, 2]
    ))
    m = lm(bold[[roi]] ~ neg_rpe$x + neg_rpe$t)
    y = rbind(y, data.frame(
      subject, roi, what = "neg_rpe",
      y = coef(m)[2], vy = vcov(m)[2, 2]
    ))
    if (!is.null(data$rpe_nc)) {
      m = lm(bold[[roi]] ~ rpe_nc$x + rpe_nc$t)
      y = rbind(y, data.frame(
        subject, roi, what = "rpe_nc",
        y = coef(m)[2], vy = vcov(m)[2, 2]
      ))
      m = lm(bold[[roi]] ~ abs_rpe_nc$x + abs_rpe_nc$t)
      y = rbind(y, data.frame(
        subject, roi, what = "abs_rpe_nc",
        y = coef(m)[2], vy = vcov(m)[2, 2]
      ))
      m = lm(bold[[roi]] ~ pos_rpe_nc$x + pos_rpe_nc$t)
      y = rbind(y, data.frame(
        subject, roi, what = "pos_rpe_nc",
        y = coef(m)[2], vy = vcov(m)[2, 2]
      ))
      m = lm(bold[[roi]] ~ neg_rpe_nc$x + neg_rpe_nc$t)
      y = rbind(y, data.frame(
        subject, roi, what = "neg_rpe_nc",
        y = coef(m)[2], vy = vcov(m)[2, 2]
      ))
    }
  }
  y
}

# meta_z() function fits a meta-analytic random-effects model
meta_z = function (y, roi, what) {
  y_i = y[which(y$roi == roi & y$what == what), ]
  z = try(rma(yi = y_i$y, vi = y_i$vy, control=list(maxiter = 1000))$zval)
  if (class(z) == "try-error") {
    z = rma(yi = y_i$y, vi = y_i$vy, method = "DL")$zval
  }
  z
}

# penalize.alpha() function computes penalization term for alpha
penalize.alpha = function (alpha, median.alpha, penalty.alpha) {
  odd.alpha = alpha / (1 - alpha)
  odd.median.alpha = median.alpha / (1 - median.alpha)
  penalty.alpha * log(odd.alpha / odd.median.alpha)^2
}

# penalize.beta() function computes penalization term for beta
penalize.beta = function (beta, median.beta, penalty.beta) {
  penalty.beta * log(beta / median.beta)^2
}

# read_behav() function transforms raw behavioral data into useful format
read_behav = function (filename) {
  data_raw = read.table(filename, as.is = TRUE)
  data = NULL
  for (stimulus in 0:9) {
    for (trial in 0:15) {
      trial_i = data_raw[(stimulus * 16 + trial) * 3 + 1:3, ]
      if (any(trial_i$V1 != c("stim", "resp", "feed"))) {
        stop("Column 1 should be a repetition of the words [stim, resp, feed]")
      }
      data = rbind(data, data.frame(
        time_stim = as.numeric(as.character(trial_i[1, 3])) / 1000, 
        time_resp = as.numeric(as.character(trial_i[2, 3])) / 1000, 
        time_feed = as.numeric(as.character(trial_i[3, 3])) / 1000, 
        stimulus,
        response = switch(as.character(trial_i[2, 2]), A = 0, B = 1, NA),
        correct = suppressWarnings(as.numeric(as.character(trial_i[3, 2])))
      ))
    }
  }
  data$time_resp = data$time_stim + data$time_resp
  data
}

# read_bold() function creates an object with region-of-interest timeseries of BOLD signal
read_bold = function (filenames, names) {
  bold = NULL
  for (i in 1:length(filenames)) {
    y = read.table(filenames[i])
    colnames(y) = names[i]
    if (is.null(bold)) {
      bold = y
    } else {
      bold = cbind(bold, y)
    }
  }
  bold
}

# ql1_fit() function creates model variables of standard Q-learning model (ql1)
ql1_fit = function (data, alpha, beta) {
  data$ev0 = NA 
  data$ev1 = NA 
  data$p0 = NA
  data$p1 = NA
  data$rpe = NA
  for (stimulus in unique(data$stimulus)) {
    rows = which(data$stimulus == stimulus)
    ev0 = 0
    ev1 = 0
    for (row in rows) {
      data$ev0[row] = ev0
      data$ev1[row] = ev1
      data$p0[row] = exp(ev0 / beta) / (exp(ev0 / beta) + exp(ev1 / beta)) 
      data$p1[row] = exp(ev1 / beta) / (exp(ev0 / beta) + exp(ev1 / beta))
      if (data$resp[row] == 0) {
        delta = data$correct[row] - ev0 
        ev0 = ev0 + alpha * delta
      } else {
        delta = data$correct[row] - ev1
        ev1 = ev1 + alpha * delta
      }
      data$rpe[row] = delta
    }
  }
  data
}

# ql1() function performs model fitting using maximum likelihood estimation 
# in addition to a penalization term
ql1 = function (data, median.alpha, penalty.alpha, median.beta, penalty.beta) {
  data = data[which(!is.na(data$resp)), ]
  init.beta = 1;
  par = optim(
    c(median.alpha, median.beta),
    function (par, data, median.alpha, median.beta) {
      alpha = par[1]
      beta = par[2]
      if (
        alpha <= 0 || alpha >= 1 ||
        beta <= 0
      ) {
        return(Inf)
      }
      cat(".")
      data.fitted = ql1_fit(data, alpha, beta)
      -(
        sum(log(data.fitted$p0[which(data.fitted$resp == 0)]))
        + sum(log(data.fitted$p1[which(data.fitted$resp == 1)]))
        - penalize.alpha(alpha, median.alpha, penalty.alpha)
        - penalize.beta(beta, median.beta, penalty.beta)
      )
    }, gr = NULL, data, median.alpha, median.beta
  )$par;
  list(
    parameters = data.frame(
      alpha = par[1],
      beta = par[2]
    ),
    data = ql1_fit(data, par[1], par[2])
  )
}

# ql1_and_feat() function combines previous functions to perform model fitting and linear regression 
ql1_and_feat = function (subject, median.alpha, penalty.alpha, median.beta, penalty.beta) {
  behav_raw = read_behav(paste("behavioral/", subject, "_RLT.txt", sep = ""))
  behav = ql1(behav_raw, median.alpha, penalty.alpha, median.beta, penalty.beta)
  cat("\n")
  behav$data$V = apply(behav$data, 1, function(x) {
    x[which(names(x) == ifelse(x[which(names(x) == "response")] == 0, "ev0", "ev1"))]
  })
  bold = read_bold(paste("time_course/", subject, "/", rois, ".txt", sep = ""), rois)
  feat(subject, behav, bold)
}

# define regions of interest
rois <- c("lVS", "rVS", "lDS", "rDS", "lOFC", "rOFC", "vmPFC", "dACC", "lAI", "rAI", "pSMA")

# define subjects IDs
subjects <- 1:30

# main loop
for (fold in 1:10) {
  id_test = (fold - 1) * 3 + 1:3
  train = subjects[-id_test]
  test = subjects[id_test]
  X = expand.grid(
    fold = fold, 
    roi = rois, 
    what = c("V", "rpe", "abs_rpe", "pos_rpe", "neg_rpe"), 
    max_z = 0,
    median.alpha = NA,
    penalty.alpha = NA,
    median.beta = NA,
    penalty.beta = NA,
    test_y1 = NA, test_vy1 = NA,
    test_y2 = NA, test_vy2 = NA,
    test_y3 = NA, test_vy3 = NA
  )
  alpha = c()
  beta = c()
  for (subject in train) {
    behav_raw = read_behav(paste("behavioral/", subject, "_RLT.txt", sep = ""))
    behav = ql1(behav_raw, 0.5, 0, 1, 0)
    alpha = c(alpha, behav$parameters$alpha)
    beta = c(beta, behav$parameters$beta)
  }
  median.alpha = median(alpha)
  median.beta = median(beta)
  
  for (penalty.alpha in 2*(1:5) / 10) { 
    for (penalty.beta in 2*(1:5) / 10) {
      cat("=================================================\n")
      cat("Fold:", fold, "\n")
      cat("- median.alpha:", median.alpha, "- penalty.alpha:", penalty.alpha, "\n")
      cat("- median.beta:", median.beta, "- penalty.beta:", penalty.beta, "\n")
      cat("=================================================\n")
      y_train = NULL
      y_test = NULL
      
      for (subject in train) {
        y_train = rbind(y_train, ql1_and_feat(subject, median.alpha, penalty.alpha, 
                                              median.beta, penalty.beta))
        }

      for (subject in test) {
      y_test = rbind(y_test, ql1_and_feat(subject, median.alpha, penalty.alpha, 
                                          median.beta, penalty.beta))
      }
      
      for (roi in rois) {
        for (what in c("V", "rpe", "abs_rpe", "pos_rpe", "neg_rpe")) {
          z = meta_z(y_train, roi, what)
          pos_y_test = which(y_test$roi == roi & y_test$what == what)
          posX = which(X$roi == roi & X$what == what)
          if (abs(z) > X$max_z[posX]) {
            X$max_z[posX] = abs(z)
            X$median.alpha[posX] = median.alpha
            X$penalty.alpha[posX] = penalty.alpha
            X$median.beta[posX] = median.beta
            X$penalty.beta[posX] = penalty.beta
            X$test_y1[posX] = y_test$y[pos_y_test][1]
            X$test_vy1[posX] = y_test$vy[pos_y_test][1]
            X$test_y2[posX] = y_test$y[pos_y_test][2]
            X$test_vy2[posX] = y_test$vy[pos_y_test][2]
            X$test_y3[posX] = y_test$y[pos_y_test][3]
            X$test_vy3[posX] = y_test$vy[pos_y_test][3]            
          }
        }
      }
      print(X)
    }
  }
  write.table(X, paste("fold", fold, "_ql1.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
}