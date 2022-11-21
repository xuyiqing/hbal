### R Code for hbal Tutorial
# v.1.2.6

library(ggplot2)
library(estimatr)

#####################
## Basic Usage
#####################

library(hbal)
set.seed(1984)
N <- 1500
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rbinom(N, size = 1, prob = .5)
D_star <- 0.5 * X1 + 0.3 * X2 + 0.2 * X1 * X2 - 0.5 * X1 * X3 - 1
D <- ifelse(D_star > rnorm(N), 1, 0) # Treatment indicator
y <- 0.5 * D + X1 + X2 + X2 * X3 + rnorm(N) # Outcome
dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, Y = y)
head(dat)

# ebal vs hbal
library(ebal)
ebal.out <- ebalance(Treat = dat$D, X = dat[,c('X1', 'X2', 'X3')]) # ebal
hbal.out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y',  data = dat) # hbal
# plot weights
W <- data.frame(x = ebal.out$w, y = hbal.out$weights.co)
ggplot(aes(x = x, y = y), data = W) + geom_point() + theme_bw() + 
  labs(x = "ebal weights", y="hbal weights", title = "ebal weights vs. hbal weights")


names(hbal.out)
summary(hbal.out)

# Adding higher-order terms
out <- hbal(Y = 'Y', Treat = 'D', X = c('X1', 'X2', 'X3'),  
            data = dat, expand.degree = 3)
summary(out)
summary(out, print.level = 1) # more info


# Obtaining the ATT
att(out, dr = FALSE)
att(out)
att(out, method = "lm_lin", se_type = "stata")

out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', 
            data = dat, expand.degree = 3, cv = TRUE)
summary(out)
att(out)

# Visualizing Results
plot(out)
plot(out, type='weight')

round(out$group.penalty, 2)
round(out$term.penalty, 2)



#########################################
## Additional Options
#########################################

# Controlling Exact/Approximate Balancing
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data = dat, 
            expand.degree = 3, cv = TRUE, group.exact = c(1, 1, 0, 0, 0))
summary(out)


# User-supplied Penalties
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, group.alpha = c(0, 0, 100, 100, 100, 100))
summary(out)


# Controling Serial Expansion
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, X.expand = c('X1', 'X2'))
summary(out)


# Selecting Covariates
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, ds = TRUE) 
summary(out)
att(out)


# Keeping/Excluding Covariates
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, exclude = list(c('X1', 'X2')))
summary(out) # X1.X2 and X1.X1.X2 removed from balancing scheme
att(out) 

#########################################
## Example 1: Lalonde Data
#########################################

data(hbal)
head(lalonde)
xvars <- c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75") # covariates

# hbal w/ level terms only
hbal.out <- hbal(Treat = 'nsw', X = xvars,  Y = 're78', data = lalonde) 
summary(hbal.out)
att(hbal.out)

# hbal w/ higher-order terms
hbal.full.out <- hbal(Treat = 'nsw', X = xvars, Y = 're78', data = lalonde, 
                      expand.degree = 2, cv = TRUE, exclude=list(c("educ", "nodegr")))
summary(hbal.full.out)
att(hbal.full.out)

hbal.full.out$group.penalty
plot(hbal.full.out)
plot(hbal.full.out, type = "weight")


#########################################
## Example 2: Black and Owens (2016)
#########################################

data(hbal)
str(contenderJudges)

xvars <- c("judgeJCS", "presDist", "panelDistJCS", "circmed", "sctmed", "coarevtc", "casepub")
out <- hbal(Treat = 'treatFinal0', X = xvars, Y = 'presIdeoVote', data = contenderJudges,
            expand.degree = 2, cv = TRUE)
summary(out)
att(out)
plot(out)
plot(out, type = "weight")
