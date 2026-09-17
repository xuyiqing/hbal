### R Code for hbal Tutorial
### Auto-generated via knitr::purl() from tutorial/*.Rmd -- do not edit by hand

## ----setup-01, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(fig.width = 10, fig.height = 7)
library(ggplot2)
require(ebal)


## ----message=FALSE------------------------------------------------------------
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


## ----fig.height = 6, fig.width = 6, fig.align = "left", dpi=100---------------
library(ebal)
ebal.out <- ebalance(Treat = dat$D, X = dat[,c('X1', 'X2', 'X3')]) # ebal
hbal.out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y',  data = dat) # hbal

# plot weights
W <- data.frame(x = ebal.out$w, y = hbal.out$weights.co)
ggplot(aes(x = x, y = y), data = W) + geom_point() + theme_bw() + 
  labs(x = "ebal weights", y="hbal weights", title = "ebal weights vs. hbal weights")


## -----------------------------------------------------------------------------
names(hbal.out)


## -----------------------------------------------------------------------------
summary(hbal.out)


## -----------------------------------------------------------------------------
out <- hbal(Y = 'Y', Treat = 'D', X = c('X1', 'X2', 'X3'),  
            data = dat, expand.degree = 3)
summary(out)


## ----bd-extract---------------------------------------------------------------
bd <- balanceData(out)
head(bd)


## ----bd-plot, fig.height = 6, fig.align = "left"------------------------------
bd$term <- factor(bd$term, levels = rev(unique(bd$term)))

ggplot(bd, aes(x = std.diff, y = term, fill = adjustment)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c(before = "white", after = "black")) +
  facet_wrap(~ covar.group, scales = "free_y") +
  labs(x = "Std. Diff.", y = NULL) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")

## ----setup-02, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(fig.width = 10, fig.height = 7)
library(hbal)
library(estimatr)
set.seed(1984)
N <- 1500
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rbinom(N, size = 1, prob = .5)
D_star <- 0.5 * X1 + 0.3 * X2 + 0.2 * X1 * X2 - 0.5 * X1 * X3 - 1
D <- ifelse(D_star > rnorm(N), 1, 0) # Treatment indicator
y <- 0.5 * D + X1 + X2 + X2 * X3 + rnorm(N) # Outcome
dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, Y = y)
out <- hbal(Y = 'Y', Treat = 'D', X = c('X1', 'X2', 'X3'),
            data = dat, expand.degree = 3)


## -----------------------------------------------------------------------------
att(out)
att(out, method = "lm_robust")
att(out, method = "lm_lin", se_type = "stata")
att(out, dr = FALSE)


## ----att-dr-false-legacy------------------------------------------------------
att(out, method = "lm_robust", dr = FALSE)
att(out, method = "lm_lin", dr = FALSE)
att(out, method = "elnet", dr = FALSE)


## -----------------------------------------------------------------------------
res <- att(out, displayAll = TRUE)
str(res, max.level = 1)
round(res$nuisance$coef_full, 3)


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', 
            data = dat, expand.degree = 3, cv = TRUE, seed = 94035)
summary(out)
att(out)


## ----fig.align = "left"-------------------------------------------------------
plot(out)


## -----------------------------------------------------------------------------
round(out$group.penalty, 2)
round(out$term.penalty, 2)


## ----fig.width = 6, fig.height = 6, fig.align = "left"------------------------
plot(out, type='weight')

## ----setup-03, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(fig.width = 10, fig.height = 7)
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


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data = dat, 
            expand.degree = 3, cv = TRUE, group.exact = c(1, 1, 0, 0, 0))
summary(out)


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, group.alpha = c(0, 0, 100, 100, 100, 100))
summary(out)


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, X.expand = c('X1', 'X2'))
summary(out)


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, ds = TRUE) 
summary(out)
att(out)


## -----------------------------------------------------------------------------
out <- hbal(Treat = 'D', X = c('X1', 'X2', 'X3'),  Y = 'Y', data = dat, 
            expand.degree = 3, exclude = list(c('X1', 'X2')))
summary(out) # X1.X2 and X1.X1.X2 removed from balancing scheme
att(out) 

## ----setup-04, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(fig.width = 10, fig.height = 7)
library(hbal)


## -----------------------------------------------------------------------------
data(hbal)
head(lalonde)


## -----------------------------------------------------------------------------
xvars <- c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75") # covariates
# hbal w/ level terms only
hbal.out <- hbal(Treat = 'nsw', X = xvars,  Y = 're78', data = lalonde) 
summary(hbal.out)
att(hbal.out)


## -----------------------------------------------------------------------------
hbal.full.out <- hbal(Treat = 'nsw', X = xvars, Y = 're78', data = lalonde, 
                      expand.degree = 2, cv = TRUE, exclude=list(c("educ", "nodegr")))
summary(hbal.full.out)
att(hbal.full.out)


## ----fig.height = 8, fig.align = "left"---------------------------------------
hbal.full.out$group.penalty
plot(hbal.full.out)

## ----setup-05, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(fig.width = 10, fig.height = 7)
library(hbal)


## -----------------------------------------------------------------------------
data(hbal)


## -----------------------------------------------------------------------------
str(contenderJudges)


## -----------------------------------------------------------------------------
xvars <- c("judgeJCS", "presDist", "panelDistJCS", "circmed", "sctmed", "coarevtc", "casepub")
out <- hbal(Treat = 'treatFinal0', X = xvars, Y = 'presIdeoVote', data = contenderJudges,
            expand.degree = 2, cv = TRUE)
summary(out)
att(out)


## ----fig.height = 8, fig.align = "left"---------------------------------------
plot(out)

