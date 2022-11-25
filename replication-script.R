########## REPLICATION CODE FOR PAPER ##############

sessionInfo()

#### CALLING MATRICES + ATTRIBUTE DATA #####

# Remove everything from the environment
remove(list = ls())

setwd("") # add wd here
getwd()

# call packages
#library(statnet)
library(network)
library(ergm)
library(btergm)
library(coda)

library(dplyr)
library(texreg)
library(kableExtra)
library(patchwork)
library(gridExtra)

library(ggplot2)
library(RColorBrewer)

# set random seed for reproducibility
set.seed(12345)
# turn off scientific notation
options(scipen = 999)


##### CONFLICT MATRIX #####

list.files(path = "./mydata")

# load conflict matrix
conflict.matrix <- as.matrix(
  read.csv(
    "./mydata/conflict_matrix_Feb22.csv",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
)
# check order is correct
head(conflict.matrix)
sum(rownames(conflict.matrix) != colnames(conflict.matrix))

# create conflict network object
conflictnet <-
  network(
    conflict.matrix,
    directed = F,
    ignore.eval = F,
    names.eval = 'battles'
  )

# check network summary
summary.network(conflictnet,
                print.adj = FALSE)


##### COOPERATION MATRIX #####

# load alliance matrix
alliance.matrix <- as.matrix(
  read.csv(
    "./mydata/allies_matrix_Feb22.csv",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
)
# check order is correct for alliances
sum(rownames(alliance.matrix) != colnames(alliance.matrix))

# reorder to match conflict network object
cvertex <- as.vector(network.vertex.names(conflictnet))
head(cvertex) # order of vertex in network object

# match alliances to conflict matrix
alliance.matrix <- alliance.matrix[cvertex, cvertex]
head(alliance.matrix)

# check if any difference
sum(rownames(conflict.matrix) != rownames(alliance.matrix))

# create alliance network object
alliancenet <-
  network(
    alliance.matrix,
    directed = F,
    ignore.eval = F,
    names.eval = 'cooperation'
  )
summary.network(alliancenet, print.adj = F)

# check for differences
avertex <- network.vertex.names(alliancenet)
setdiff(rownames(conflict.matrix), rownames(alliance.matrix))


#### ATTRIBUTE DATA ####

# load attribute data
attrs <- read.csv(
  "mydata/attribute_data_Feb22.csv",
  header = T,
  stringsAsFactors = F,
  row.names = 1
)
glimpse(attrs)

# rearrange attributes by conflict node order
attrs <- attrs[match(cvertex, attrs$actors), ]
head(attrs)[1:4,]

# check if difference
setdiff(cvertex, attrs$actors)


### IDENTITY SIMILARITY MATRIX ###

# similarity identity matrix (5 relevant cleavages)
library(proxy)
names(attrs)

sim <- function(x, y)
  sum(x == y) - sum(x == 0 & y == 0)

simidmat <- as.matrix(dist(attrs[c(14:18)], method = sim))

diag(simidmat) = 0
simidmat[1:5, 1:5]
rownames(simidmat) <- attrs[[2]]
colnames(simidmat) <- attrs[[2]]
# view matrices
head.matrix(simidmat)[1:5]

#write.csv(simidmat, "similarity_matrix_Feb22.csv")


##### ATTACH ATTRIBUTES TO NETWORK #####
library(statnet)

# Add the vertex attributes to conflict network
set.vertex.attribute(conflictnet, "abbrv", as.character(attrs$abbrv))

# add controls
set.vertex.attribute(conflictnet, "sponsor", as.numeric(attrs$sponsor))
set.vertex.attribute(conflictnet, "coalition", as.numeric(attrs$coalition))
set.vertex.attribute(conflictnet, "territory", as.numeric(attrs$territory))
set.vertex.attribute(conflictnet, "geopres", as.numeric(attrs$geospr))
set.vertex.attribute(conflictnet, "years", as.numeric(attrs$years))
set.vertex.attribute(conflictnet, "powproj", as.numeric(attrs$powproj))

# check vertex attributes attached properly
list.vertex.attributes(conflictnet) # 9 attributes
list.edge.attributes(conflictnet) # 2 attributes

# check networks/matrices
summary.network(conflictnet, print.adj = F)


#### STATISTICAL MODELS ####

##### MODEL 1: alliances #####
m1 <- ergm(
  conflictnet ~ edges +
    edgecov(alliancenet) +
    # controls
    nodefactor('sponsor') +
    nodefactor('coalition')  +
    nodefactor('territory') +
    nodecov('years') +
    nodecov('geopres') +
    nodecov('powproj') +
    # endogenous dependencies
    isolates() +
    gwdegree(0.2, fixed = T) +
    gwesp(0, fixed = T),
  control = control.ergm(
    seed = 123,
    MCMC.burnin = 10000,
    MCMC.samplesize = 20000
  )
)
summary(m1) # AIC 1186
m1$mle.lik[1] # -582.2243

# check MCMC diagnostics
#mcmc.diagnostics(m1)

##### MODEL 2: shared identity #####

m2 <- ergm(
  conflictnet ~ edges +
    # controls
    nodefactor('sponsor') +
    nodefactor('coalition')  +
    nodefactor('territory') +
    nodecov('years') +
    nodecov('geopres') +
    nodecov('powproj') +
    # dyadic cov = co-identity
    edgecov(simidmat) +
    # endogenous terms
    isolates() +
    gwdegree(0.2, fixed = T) +
    gwesp(0, fixed = T),
  control = control.ergm(
    seed = 123,
    MCMC.burnin = 10000,
    MCMC.samplesize = 20000
  )
)
summary(m2) # AIC 1182
m2$mle.lik[1] # -580.0996

# check MCMC diagnostics
#mcmc.diagnostics(m2)

##### MODEL 3: all + shared identity #####
m3 <- ergm(
  conflictnet ~ edges +
    nodefactor('sponsor') +
    nodefactor('coalition')  +
    nodefactor('territory') +
    nodecov('years') +
    nodecov('geopres') +
    nodecov('powproj') +
    # dyadic covariates
    edgecov(alliancenet) +
    edgecov(simidmat) +
    # endogenous
    isolates() +
    gwdegree(0.2, fixed = T) +
    gwesp(0, fixed = T),
  # increase MCMC burnin and sample size
  control = control.ergm(
    seed = 123,
    MCMC.burnin = 10000,
    MCMC.samplesize = 20000
  )
)
summary(m3) # AIC 1180
m3$mle.lik[1] # -577.7896

# check MCMC diagnostics
#mcmc.diagnostics(m3)


#### ROBUSTNESS MODELS FOR APPENDIX ONLY ####

##### MODEL 4: power projection ####
m4 <- ergm(
  conflictnet ~ edges +
    nodefactor('sponsor') +
    nodefactor('coalition')  +
    nodefactor('territory') +
    nodecov('years') +
    nodecov('geopres') +
    # covariate of interest
    nodecov('powproj') +
    absdiff('powproj') +
    # endogenous dependencies
    isolates() +
    gwdegree(0.2, fixed = T) +
    gwesp(0, fixed = T),
  control = control.ergm(
    seed = 123,
    MCMC.burnin = 10000,
    MCMC.samplesize = 20000
  )
)
summary(m4) # 1191
m4$mle.lik[1] # -584.7242

# check MCMC diagnostics
#mcmc.diagnostics(m4)

##### MODEL 5: all + power projection #####
m5 <- ergm(
  conflictnet ~ edges +
    nodefactor('sponsor') +
    nodefactor('coalition')  +
    nodefactor('territory') +
    nodecov('years') +
    nodecov('geopres') +
    # covariates of interest
    nodecov('powproj') +
    absdiff('powproj') +
    # dyadic covariates
    edgecov(alliancenet) +
    edgecov(simidmat) +
    # endogenous dependencies
    isolates() +
    gwdegree(0.2, fixed = T) +
    gwesp(0, fixed = T),
  control = control.ergm(
    seed = 123,
    MCMC.burnin = 10000,
    MCMC.samplesize = 20000
  )
)
summary(m5) # 1181  
m5$mle.lik[1] # -577.5788

# check MCMC diagnostics
#mcmc.diagnostics(m4)


#### MODEL TABLES + COMPARISON ####

##### Main table of results #####
library(texreg)
mytable = screenreg(
  list(m1, m2, m3),
  single.row = F,
  ci.force = F,
  digits = 2
)
mytable

##### Model 3 odds + 95% CI #####

# manual table for model 3
or <- exp(coef(m3)) # odds ratio
ste <- sqrt(diag(m3$covar)) # st error
lci <- exp(coef(m3) - 1.96 * ste) # low ci
uci <- exp(coef(m3) + 1.96 * ste) # high ci
oddsratios <- rbind(round(or, digits = 4), # bind
                    round(lci, digits = 4),
                    round(uci, digits = 4))
oddsratios <- t(oddsratios) # table
colnames(oddsratios) <- c("Odds Ratio",
                          "Lower",
                          "Upper") # cols
oddsr <- as.data.frame(oddsratios)
oddsr

# extract CI values
ci.low <- exp(coef(m3) - 1.96 * (sqrt(diag(m3$covar))))
ci.high <- exp(coef(m3) + 1.96 * (sqrt(diag(m3$covar))))

# odds for all models
odds1 <- exp(coef(m1)) #informal cooperation
odds2 <- exp(coef(m2)) #shared identity
odds3 <- exp(coef(m3)) #cooperation+identity


##### MCMC DIAGNOSTICS #####

# using coda
#install.packages("coda")
library(coda)
mcmc.m3 <- mcmc(m3$sample)
str(mcmc.m3)

# run centered diagnostics
codamenu() # geweke, crosscorr, autocorr

dev.off()

##### GOODNESS-OF-FIT COMPARISON #####

# gof plots for model 3 - use btergm
m5_gof <- btergm::gof(m3,
                      nsim = 5000,
                      statistics = c(deg,
                                     esp,
                                     dsp,
                                     geodesic,
                                     rocpr))

# set graphic parameters: mai = c(bottom, left, top, right)
par(mfrow = c(2, 2), mai = c(.8, .7, .1, .1))
plot(m5_gof[[1]], main = "")
plot(m5_gof[[2]], main = "")
plot(m5_gof[[3]], main = "")
plot(m5_gof[[4]], main = "")

dev.off()


##### MODEL 3 NETWORK SIMULATIONS #####
set.seed(1234)
m3.sims <- simulate(m3, nsim = 10)

par(mfrow = c(1, 2), mar = c(0.8, 0.2, 0.8, 0.2))
sapply(m3.sims[1], plot, vertex.cex = 1, vertex.col = "tomato")
plot(conflictnet, vertex.cex = 1, vertex.col = "blue")

dev.off()

##### PREDICTED PROBABILITY PLOT: IDENTITY #####

# frequency table for shared identity
ep <- edgeprob(m3)
names(ep)

# frequency table for identity
tab <-
  data.frame(table(ep$`edgecov.simidmat[[i]]`)) # interpretation: shared dimensions across each dyad
tab

# plot frequency and range barplot
fr <- ggplot(tab,
             aes(Var1, Freq)) +
  geom_bar(stat = "identity",
           fill = "gray55", alpha = 0.8) +
  xlab("Shared identity") +
  ylab("Frequency") +
  labs(title = "Frequency of shared identity between armed actors") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 12)) +
  coord_cartesian(ylim = c(1, 10500),
                  expand = T)
print(fr)

# predicted probabilities
pp <- ggplot(ep,
             aes(x = `edgecov.simidmat[[i]]`,
                 y = probability)) +
  stat_summary(
    geom = "ribbon",
    fun.data = mean_cl_normal,
    fun.args = (conf.int = 0.95),
    fill = "gray65",
    alpha = 0.4
  ) +
  stat_summary(geom = "line", fun = mean) +
  stat_summary(
    geom = "point",
    fun = mean,
    size = 0.8,
    color = "black"
  ) +
  ylab("Predicted probability") +
  xlab("Shared identity") +
  labs(title = "Predicted probabilities for shared identity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 12)) +
  coord_cartesian(ylim = c(0.0, 0.03),
                  xlim = c(0, 4),
                  expand = T)
print(pp)

# save figures to pdf
grid.arrange(fr, pp, ncol = 2)

dev.off()


##### MARGINAL EFFECTS PLOT #####

alldyads <- edgeprob(m3)
allnames <- matrix(data = names(alldyads), ncol = 1)
colnames(allnames) <- "Variable name"
kable_styling(kable(allnames),
              bootstrap_options = c("condense"),
              full_width = FALSE)

colnames(alldyads)
head(alldyads[, c("tie", "i", "j", "probability")])

# plotting
marplot <- ggplot(alldyads,
                  aes(
                    x = `edgecov.simidmat[[i]]`,
                    y = probability,
                    colour = factor(`edgecov.alliancenet[[i]]`)
                  )) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm",
              alpha = 0.1,
              fullrange = F) +
  labs(
    title = "Marginal Effect of Identity by Tactical Cooperation",
    x = "Shared Identity",
    y = "Probability of Conflict Tie",
    colour = "Cooperation"
  ) +
  theme(legend.position = "bottom")
marplot


# marginal effect of identity by coalition-presence
marplot2 <- ggplot(alldyads,
                   aes(
                     x = `edgecov.simidmat[[i]]`,
                     y = probability,
                     color = factor(nodefactor.coalition.1)
                   )) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm",
              alpha = 0.1,
              fullrange = F) +
  coord_cartesian(ylim = c(0.00, 1.0),
                  expand = T) +
  labs(
    title = "Marginal Effect of Identity by Coalition Membership",
    x = "Shared Identity",
    y = NULL,
    colour = "Membership"
  ) +
  theme(legend.position = "bottom")
marplot2

# save figures to pdf
fullmar <-
  marplot + scale_colour_brewer(palette = "Set2", labels = c("no", "yes")) +
  marplot2 + scale_colour_brewer(palette = "Set2", labels = c("none", "one", "both")) &
  theme_bw() + theme(legend.position = "bottom")
fullmar

dev.off()

### END OF CODE ###