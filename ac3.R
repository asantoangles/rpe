library(metafor)
library(Matrix)

# call_feat() function convolves model variables to hemodynamic response function using FSL FEAT to get model regresors
call_feat = function (time, duration, value) {
  write.table(data.frame(time, duration, value), "faked_feat/ac3_new/rpe.txt",
              col.names = FALSE, row.names = FALSE)
  system("rm -rf faked_feat/ac3_new/result*.feat")
  system("feat faked_feat/ac3_new/faked_feat.fsf", ignore.stdout = TRUE)
  X = read.table("faked_feat/ac3_new/result.feat/design.mat", skip = 5)
  colnames(X) = c("x", "t")
  X
}

# feat() function performs linear regression of bold signal (as dependent variable) and model regresor (as independent variable)
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
    z = rma(yi = y_i$y, vi = y_i$vy, control=list(maxiter = 1000), method = "DL")$zval
  }
  z
}

# penalize.alpha.critic() function computes penalization term for alpha.critic
penalize_alpha.critic = function (alpha.critic, median.alpha.critic, penalty.alpha.critic) {
  odd.alpha.critic = alpha.critic / (1 - alpha.critic)
  odd.median.alpha.critic = median.alpha.critic / (1 - median.alpha.critic)
  penalty.alpha.critic * log(odd.alpha.critic / odd.median.alpha.critic)^2
}

# penalize.alpha.actor.ch() function computes penalization term for alpha.actor.ch
penalize_alpha.actor.ch = function (alpha.actor.ch, median.alpha.actor.ch, penalty.alpha.actor.ch) {
  odd.alpha.actor.ch = alpha.actor.ch / (1 - alpha.actor.ch)
  odd.median.alpha.actor.ch = median.alpha.actor.ch / (1 - median.alpha.actor.ch)
  penalty.alpha.actor.ch * log(odd.alpha.actor.ch / odd.median.alpha.actor.ch)^2
}

# penalize.alpha.actor.nc() function computes penalization term for alpha.actor.nc
penalize_alpha.actor.nc = function (alpha.actor.nc, median.alpha.actor.nc, penalty.alpha.actor.nc) {
  odd.alpha.actor.nc = alpha.actor.nc / (1 - alpha.actor.nc)
  odd.median.alpha.actor.nc = median.alpha.actor.nc / (1 - median.alpha.actor.nc)
  penalty.alpha.actor.nc * log(odd.alpha.actor.nc / odd.median.alpha.actor.nc)^2
}

# penalize.beta() function computes penalization term for beta
penalize_beta = function (beta, median.beta, penalty.beta) {
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


# ac3_fit() function creates model variables of extended actor-critic model assuming alpha.actor.ch != alpha.actor.nc (ac2)
ac3_fit = function (data, alpha.critic, alpha.actor.ch, alpha.actor.nc, beta) {
  data$V = NA 
  data$w0 = NA
  data$w1 = NA 
  data$p0 = NA
  data$p1 = NA
  data$rpe = NA
  data$rpe_nc = NA
  for (stimulus in unique(data$stimulus)) {
    rows = which(data$stimulus == stimulus)
    V = 0 
    w0 = 0.01 
    w1 = 0.01 
    for (row in rows) {
      data$V[row] = V
      data$w0[row] = w0
      data$w1[row] = w1
      data$p0[row] = exp(w0 / beta) / (exp(w0 / beta) + exp(w1 / beta)) 
      data$p1[row] = exp(w1 / beta) / (exp(w0 / beta) + exp(w1 / beta))
      delta.ch = data$correct[row] - V
      delta.nc = (1 - data$correct[row]) - V 
      V = V + alpha.critic * delta.ch
      if (data$resp[row] == 0) {
        w0 = w0 + alpha.actor.ch * delta.ch
        w1 = w1 + alpha.actor.nc * delta.nc
      } else {
        w1 = w1 + alpha.actor.ch * delta.ch
        w0 = w0 + alpha.actor.nc * delta.nc
      }
      sum_abs_w = abs(w0)+abs(w1)
      w0 = w0 / sum_abs_w
      w1 = w1 / sum_abs_w
      data$rpe[row] = delta.ch
      data$rpe_nc[row] = delta.nc
    }
  }
  data
}

# ac3() function performs model fitting using maximum likelihood estimation in addition to a penalization term
ac3 = function (data, median.alpha.critic, penalty.alpha.critic, median.alpha.actor.ch, penalty.alpha.actor.ch, median.alpha.actor.nc, penalty.alpha.actor.nc, median.beta, penalty.beta) {
  data = data[which(!is.na(data$resp)), ]
  init.beta = 1;
  par = optim(
    c(median.alpha.critic, median.alpha.actor.ch, median.alpha.actor.nc, median.beta),
    function (par, data, median.alpha.critic, median.alpha.actor.ch, median.alpha.actor.nc, median.beta) {
      alpha.critic = par[1]
      alpha.actor.ch = par[2]
      alpha.actor.nc = par[3]
      beta = par[4]
      if (
        alpha.critic <= 0 || alpha.critic >= 1 ||
        alpha.actor.ch <= 0 || alpha.actor.ch >= 1 ||
        alpha.actor.nc <= 0 || alpha.actor.nc >= 1 ||
        beta <= 0
      ) {
        return(Inf)
      }
      cat(".")
      data.fitted = ac3_fit(data, alpha.critic, alpha.actor.ch, alpha.actor.nc, beta)
      -(
        sum(log(data.fitted$p0[which(data.fitted$resp == 0)]))
        + sum(log(data.fitted$p1[which(data.fitted$resp == 1)]))
        - penalize_alpha.critic(alpha.critic, median.alpha.critic, penalty.alpha.critic)
        - penalize_alpha.actor.ch(alpha.actor.ch, median.alpha.actor.ch, penalty.alpha.actor.ch)
        - penalize_alpha.actor.nc(alpha.actor.nc, median.alpha.actor.nc, penalty.alpha.actor.nc)
        - penalize_beta(beta, median.beta, penalty.beta)
      ) 
    }, gr = NULL, data, median.alpha.critic, median.alpha.actor.ch, median.alpha.actor.nc, median.beta
  )$par;
  list(
    parameters = data.frame(
      alpha.critic = par[1],
      alpha.actor.ch = par[2],
      alpha.actor.nc = par[3],
      beta = par[4]
    ),
    data = ac3_fit(data, par[1], par[2], par[3], par[4])
  )
}

# ac3_and_feat() function combines previous functions to perform model fitting and linear regression 
ac3_and_feat = function (subject, median.alpha.critic, penalty.alpha.critic, median.alpha.actor.ch, penalty.alpha.actor.ch, median.alpha.actor.nc, penalty.alpha.actor.nc, median.beta, penalty.beta) {
  behav_raw = read_behav(paste("behavioral/", subject, "_RLT.txt", sep = ""))
  behav = ac3(behav_raw, median.alpha.critic, penalty.alpha.critic, median.alpha.actor.ch, penalty.alpha.actor.ch, median.alpha.actor.nc, penalty.alpha.actor.nc, median.beta, penalty.beta)
  cat("\n")
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
    fold = fold, roi = rois, what = c("V", "rpe", "abs_rpe", "pos_rpe", "neg_rpe", "rpe_nc", "abs_rpe_nc", "pos_rpe_nc", "neg_rpe_nc"), max_z = 0,
    median.alpha.critic = NA,
    penalty.alpha.critic = NA,
    median.alpha.actor.ch = NA,
    penalty.alpha.actor.ch = NA,
    median.alpha.actor.nc = NA,
    penalty.alpha.actor.nc = NA,
    median.beta = NA,
    penalty.beta = NA,
    test_y1 = NA, test_vy1 = NA,
    test_y2 = NA, test_vy2 = NA,
    test_y3 = NA, test_vy3 = NA
  )
  alpha.critic = c()
  alpha.actor.ch = c()
  alpha.actor.nc = c()
  beta = c()
  
  for (subject in train) {
    behav_raw = read_behav(paste("behavioral/", subject, "_RLT.txt", sep = ""))
    behav = ac3(behav_raw, 0.5, 0, 0.5, 0, 0.5, 0, 1, 0)
    alpha.critic = c(alpha.critic, behav$parameters$alpha.critic)
    alpha.actor.ch = c(alpha.actor.ch, behav$parameters$alpha.actor.ch)
    alpha.actor.nc = c(alpha.actor.nc, behav$parameters$alpha.actor.nc)
    beta = c(beta, behav$parameters$beta)
    
  }
  median.alpha.critic = median(alpha.critic)
  median.alpha.actor.ch = median(alpha.actor.ch)
  median.alpha.actor.nc = median(alpha.actor.nc)
  median.beta = median(beta)
  
  for (penalty.alpha.critic in 2*(1:5) / 10) { 
    for (penalty.alpha.actor.ch in 2*(1:5) / 10) { 
      for (penalty.alpha.actor.nc in 2*(1:5) / 10) { 
        for (penalty.beta in 2*(1:5) / 10) { 
          cat("=================================================\n")
          cat("Fold:", fold, "\n")
          cat("- median.alpha.critic:", median.alpha.critic, "- penalty.alpha.critic:", penalty.alpha.critic, "\n")
          cat("- median.alpha.actor.ch:", median.alpha.actor.ch, "- penalty.alpha.actor.ch:", penalty.alpha.actor.ch, "\n")
          cat("- median.alpha.actor.nc:", median.alpha.actor.nc, "- penalty.alpha.actor.nc:", penalty.alpha.actor.nc, "\n")
          cat("- median.beta:", median.beta, "- penalty.beta:", penalty.beta, "\n")
          cat("=================================================\n")
          y_train = NULL
          y_test = NULL
          
          for (subject in train) {
            y_train = rbind(y_train, ac3_and_feat(subject, median.alpha.critic, penalty.alpha.critic, median.alpha.actor.ch, penalty.alpha.actor.ch, median.alpha.actor.nc, penalty.alpha.actor.nc, median.beta, penalty.beta))
          print(y_train)
            }
          for (subject in test) {
            y_test = rbind(y_test, ac3_and_feat(subject, median.alpha.critic, penalty.alpha.critic, median.alpha.actor.ch, penalty.alpha.actor.ch, median.alpha.actor.nc, penalty.alpha.actor.nc, median.beta, penalty.beta))
          print(y_test)
            }
          for (roi in rois) {
            for (what in c("V", "rpe", "abs_rpe", "pos_rpe", "neg_rpe", "rpe_nc", "abs_rpe_nc", "pos_rpe_nc", "neg_rpe_nc")) {
              z = meta_z(y_train, roi, what)
              pos_y_test = which(y_test$roi == roi & y_test$what == what)
              posX = which(X$roi == roi & X$what == what)
              if (abs(z) > X$max_z[posX]) {
                X$max_z[posX] = abs(z)
                X$median.alpha.critic[posX] = median.alpha.critic
                X$penalty.alpha.critic[posX] = penalty.alpha.critic
                X$median.alpha.actor.ch[posX] = median.alpha.actor.ch
                X$penalty.alpha.actor.ch[posX] = penalty.alpha.actor.ch
                X$median.alpha.actor.nc[posX] = median.alpha.actor.nc
                X$penalty.alpha.actor.nc[posX] = penalty.alpha.actor.nc
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
    }
  }
  write.table(X, paste("fold", fold, "_ac3.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
}