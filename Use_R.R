library("MASS")
library("mlogit")
library("ChoiceModelR")
library("dplyr")

#13 Choice Modeling

##13.1 Choice-Based Conjoint Analysis Surveys
cbc.df <- read.csv("http://goo.gl/5xQObB", colClasses = c(seat = "factor", price = "factor", choice="integer"))
cbc.df$eng <- factor(cbc.df$eng, levels = c("gas", "hyb", "elec"))
cbc.df$carpool <- factor(cbc.df$carpool, levels = c("yes", "no"))
summary(cbc.df)

##13.2 Simulating Choice Data
attrib <- list(seat = c("6", "7", "8"),
               cargo = c("2ft", "3ft"),
               eng = c("gas", "hyb", "elec"),
               price = c("30", "35", "40"))

coef.names <- NULL
for (a in seq_along(attrib)) {
  coef.names <- c(coef.names, paste(names(attrib)[a], attrib[[a]][-1], sep=""))
}
coef.names

mu <- c(-1, -1, 0.5, -1, -2, -1, -2)
names(mu) <- coef.names
mu

Sigma <- diag(c(0.3, 1, 0.1, 0.3, 1, 0.2, 0.3))
dimnames(Sigma) <- list(coef.names, coef.names)
Sigma["enghyb", "engelec"] <- Sigma["engelec", "enghyb"] <- 0.3

set.seed(33040)
resp.id <- 1:200
carpool <- sample(c("yes", "no"), size = length(resp.id), replace = T, prob = c(0.3, 0.7))
coefs <- mvrnorm(length(resp.id), mu=mu, Sigma=Sigma)
coefs[carpool=="yes", "seat8"] <- coefs[carpool=="yes", "seat8"] + 2
coefs[carpool=="yes", "seat7"] <- coefs[carpool=="yes", "seat7"] + 1.5

nques <- 15
nalt <- 3

profiles <- expand.grid(attrib)
nrow(profiles)
head(profiles)

profiles.coded <- model.matrix(~seat+cargo+eng+price, data=profiles)[, -1]
head(profiles.coded)

cbc.df <- data.frame(NULL)
for (i in seq_along(resp.id)) {
  profiles.i <- sample(1:nrow(profiles), size=nques*nalt)
  utility <- profiles.coded[profiles.i, ] %*% coefs[i, ]
  wide.util <- matrix(data = utility, ncol=nalt, byrow = T)
  probs <- exp(wide.util) /rowSums(exp(wide.util))
  choice <- apply(probs, 1, function(x) sample(1:nalt, size=1, prob=x))
  choice <- rep(choice, each=nalt)==rep(1:nalt,nques)
  conjoint.i <- data.frame(resp.id=rep(i, nques), 
                           ques = rep(1:nques, each=nalt),
                           carpool = rep(carpool[i], nques),
                           profiles[profiles.i, ],
                           choice = as.numeric(choice))
  cbc.df <- rbind(cbc.df, conjoint.i)
}
rm(a, i, resp.id, carpool, mu, Sigma, coefs, coef.names, conjoint.i, profiles, profiles.i, profiles.coded, utility, wide.util,probs, choice, nalt, nques)
head(cbc.df)
#13.3.1 Inspecting Choice Data
summary(cbc.df)
xtabs(choice ~ price, data = cbc.df)
xtabs(choice ~ cargo, data = cbc.df)

#13.3.2 Fitting Choice Models with mlogit()
cbc.mlogit <- mlogit.data(data=cbc.df, choice = "choice", shape = "long", varying = 3:6, alt.levels=paste("pos", 1:3), id.var="resp.id")
m1 <- mlogit(choice ~ 0 + seat + cargo + eng + price, data = cbc.mlogit)
summary(m1)

m2 <- mlogit(choice ~ seat + cargo + eng + price, data = cbc.mlogit)
summary(m2)

lrtest(m1, m2)

m3 <- mlogit(choice ~ 0 + seat + cargo + eng + as.numeric(as.character(price)), data = cbc.mlogit)
summary(m3)

lrtest(m1, m3)

#13.3.3 Reporting Choice Model Findings
coef(m3)["cargo3ft"]/(-coef(m3)["as.numeric(as.character(price))"]/1000)

predict.mnl <- function(model, data) {
  data.model <- model.matrix(update(model$formula,0~ .), data = data)[,-1] 
  utility <- data.model %*% model$coef
  share <- exp(utility)/sum(exp(utility))
  cbind(share, data)
}

(new.data <- expand.grid(attrib)[c(8, 1, 3, 41, 49, 26), ])

predict.mnl(m3, new.data)

predict.mnl(m1, new.data)

sensitivity.mnl <- function(model, attrib, base.data, competitor.data) {
  data <- rbind(base.data, competitor.data)
  base.share <- predict.mnl(model, data)[1,1]
  share <- NULL
  for (a in seq_along(attrib)) {
    for (i in attrib[[a]]) {
      data[1,] <- base.data
      data[1,a] <- i
      share <- c(share, predict.mnl(model, data)[1,1])
    }
  }
  data.frame(level=unlist(attrib), share = share, increase = share-base.share)
}

base.data <- expand.grid(attrib)[c(8), ]
competitor.data <- expand.grid(attrib)[c(1, 3, 41, 49, 26), ]
(tradeoff <- sensitivity.mnl(m1, attrib, base.data, competitor.data))

#13.3.4 Share Predictions for Identical Alternatives
new.data.2 <- expand.grid(attrib)[c(8, 8, 1, 3, 41, 49, 26), ]
predict.mnl(m1, new.data.2)

#13.3.5 Planning the Sample Size for a Conjoint Study
small.cbc <- mlogit.data(data = cbc.df[1:(25*15*3), ],
                         choice = "choice", shape = "long", 
                         varying = 3:6, alt.levels = paste("pos", 1:3),
                         id.var = "resp.id")
m4 <- mlogit(choice ~ 0 + seat + cargo + eng + price, data = small.cbc)

summary(m4)

cbind(predict.mnl(m4, new.data), predict.mnl(m1, new.data))

#13.4 Adding Consumer Heterogenity to Choice Models
#13.4.1 Estimating Mixed Logit Models with mlogit()
m1.rpar <- rep("n", length = length(m1$coef))
names(m1.rpar) <- names(m1$coef)
m1.rpar

m1.hier <- mlogit(choice ~ 0 + seat + eng + cargo + price, data = cbc.mlogit, panel = T, rpar = m1.rpar, correlation = F)

summary(m1.hier)

stdev(m1.hier)

m2.hier <- update(m1.hier, correlation = T)
summary(m2.hier)

cov2cor(cov.mlogit(m2.hier))

#13.4.2 Share Prediction for Heterogeneous Choice Models
predict.hier.mnl <- function(model, data, nresp = 1000) {
  data.model <- model.matrix(update(model$formula, 0 ~ .), data = data) [, -1]
  coef.Sigma <- cov.mlogit(model)
  coef.mu <- model$coef[1:dim(coef.Sigma)[1]]
  draws <- mvrnorm(n = nresp, coef.mu, coef.Sigma)
  shares <- matrix(NA, nrow = nresp, ncol = nrow(data))
  for (i in 1:nresp) {
    utility <- data.model %*% draws[i, ]
    share = exp(utility) / sum(exp(utility))
    shares[i, ] <- share
  }
  cbind(colMeans(shares), data)
}

predict.hier.mnl(m2.hier, data = new.data)

#13.5 Hierarchical Bayes Choice Models
#13.5.1 Estimating Hierarchical Bayes Choice Models with ChoiceModelR
cbc.df <- read.csv("http://goo.gl/5xQObB", colClasses = c(seat = "factor", price = "factor", choice="integer"))
choice <- rep(0, nrow(cbc.df))
choice[cbc.df[,"alt"]==1] <- cbc.df[cbc.df[,"choice"]==1, "alt"] 
head(choice)

cbc.coded <- model.matrix(~seat + eng + cargo + price, data = cbc.df)
cbc.coded <- cbc.coded[, -1]

choicemodelr.data <- cbind(cbc.df[,1:3], cbc.coded, choice)
head(choicemodelr.data)

carpool <- cbc.df$carpool[cbc.df$ques==1 & cbc.df$alt==1]=="yes"
carpool <- as.numeric(carpool)
choicemodelr.demos <- as.matrix(carpool, nrow=length(carpool))
str(choicemodelr.demos)

hb.post <- choicemodelr(data=choicemodelr.data, xcoding=rep(1,7), demos = choicemodelr.demos, mcmc = list(R = 20000, use = 10000), options = list(save = T))

names(hb.post)

hb.post$compdraw[[567]]$mu

hb.post$deltadraw[567,]

hb.post$compdraw[[567]]$rooti

crossprod(hb.post$compdraw[[567]]$rooti)

head(hb.post$betadraw[,,567])

str(hb.post$betadraw)

beta.post.mean <- apply(hb.post$betadraw, 1:2, mean)
head(beta.post.mean)

beta.post.q05 <- apply(hb.post$betadraw, 1:2, quantile, probs=c(0.05))
beta.post.q95 <- apply(hb.post$betadraw, 1:2, quantile, probs=c(0.95))
rbind(q05=beta.post.q05[1,], mean=beta.post.mean[1,], q95=beta.post.q95[1,])

#13.5.2 Share Prediction for Hierarchical Bayes Choice Models
predict.hb.mnl <- function(betadraws, data) {
  data.model <- model.matrix(~seat + eng + cargo + price, data = data)
  data.model <- data.model[,-1]
  nresp <- dim(betadraws)[1]
  ndraws <- dim(hb.post$betadraw)[3]
  shares <- array(dim=c(nresp, nrow(data), ndraws))
  for (d in 1:ndraws) {
    for (i in 1:nresp) {
      utility <- data.model %*% betadraws[i,,d]
      shares[i,,d] = exp(utility) / sum(exp(utility))
    }
  }
  shares.agg <- apply(shares, 2:3, mean)
  cbind(share = apply(shares.agg, 1, mean),
        pct = t(apply(shares.agg, 1, quantile, probs = c(0.05, 0.95))), data)
}

predict.hb.mnl(hb.post$betadraw, new.data)

#13.6 Design of Choice-Based Conjoint Surveys

#13.7 Key Points

#13.8 Data Sources

#13.9 Learning More

#13.10 Exercises
##Aufgabe 1
sportscar <- read.csv("https://goo.gl/8g7vtT")
summary(sportscar)
glimpse(sportscar)
xtabs(~trans + choice, data = sportscar)

##Aufgabe 2
sportscar$seat <- as.factor(sportscar$seat)
sportscar$price <- as.factor(sportscar$price)

sportscar.mlogit <- mlogit.data(data = sportscar, choice = "choice", shape = "long", varying = 3:6, alt.levels=paste("pos", 1:3), id.var="resp_id")
msc <- mlogit(choice ~ 0 + seat + trans + convert + price, data = sportscar.mlogit)
summary(msc)

newcars <- data.frame(seat=factor(c("2","4", "5")),
                      trans=factor(c("manual", "automatic", "automatic")),
                      convert=factor(c("no", "yes", "no")),
                      price=c("40", "30", "35"))

predict.mnl(msc, newcars)

