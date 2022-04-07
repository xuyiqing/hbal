## hbal examples
# by Yiqing Xu and Eddie Yang

#######################
## Simulated Example 
#######################

library(hbal)
set.seed(1984)
N <- 500
X1 <- rnorm(N)
X2 <- rbinom(N,size=1,prob=.5)
D_star <- 0.7 * X1 + 0.3 * X2
D <- ifelse(D_star > rnorm(N), 1, 0) # Treatment indicator
y <- 0.5 * D + X1 + X2 + rnorm(N) # Outcome
dat <- data.frame(D=D, X1 = X1, X2= X2, Y=y)
head(dat)

# regression
library(estimatr)
summary(lm_robust(Y ~ D + X1 + X2, data = dat, se_type = "stata"))

# hbal
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data=dat)
att(out, se_type = "stata")
att(out, dr = FALSE)


str(out)
head(out$weights)
length(out$weights)
summary(out)
out$penalty

class(att(out))
summary(att(out))

plot(out)
plot(out, type='weight')

#######################
## Lalonde 
#######################

data(hbal)
head(lalonde)
xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")

# ebal
library(ebal)
ebal.out <- ebalance(Treat = lalonde$nsw, X = lalonde[,xvars], print.level=-1)

# hbal w/ level terms only
hbal.out <- hbal(Treat = 'nsw', X = xvars,  Y = 're78', data=lalonde, 
                 expand.degree = 0, cv = FALSE) 

# save weights from ebal and hbal
W <- data.frame(x = hbal.out$weights*sum(D), 
                y = ebal.out$w) # store weights as x-y coordinates

# plot weights
ggplot(aes(x = x, y = y), data = W) + geom_point() + theme_bw() + 
  labs(x="hbal weights", y="ebal weights", 
   title="hbal weights vs hbal weights")

summary(hbal.out)

hbal.full.out <- hbal(Treat = 'nsw', X = xvars, Y = 're78', data=lalonde,
                      expand.degree = 2, exclude=list(c("educ", "nodegr")))
summary(hbal.full.out)$'ATT Estimate'

plot(hbal.full.out)
hbal.full.out$penalty

###########################
## Black and Owens (2016)
###########################

str(contenderJudges)

out <- hbal(Treat = 'treatFinal0', 
  X = c('judgeJCS','presDist','panelDistJCS','circmed','sctmed','coarevtc','casepub'),
  Y = 'presIdeoVote', data=contenderJudges)
summary(out)$`ATT Estimate`

plot(out)

#######################
## More Options
#######################

# set panelity for specific terms
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data=dat, 
            alpha=c('X1.X1.X2'=0))
out$penalty
plot(out)

# set panelity for groups
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data=dat, 
            group.alpha=c(0, 0, 0, 0, 0))
out$penalty
plot(out)

# excluding certain features
out <- hbal(Treat = 'D', X = c('X1', 'X2'),  Y = 'Y', data=dat, 
            exclude=list(c("X1", "X2")))
summary(att(out)) # X1.X2 and X1.X1.X2 removed from balancing scheme

