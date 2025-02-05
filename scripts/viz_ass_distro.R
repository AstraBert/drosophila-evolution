#### Initial cleaning and general options #################
# Clean all objects in the global environment
rm(list = ls())
# Free memory used by deleted objects
gc()
# Clean all open graphics (if using RStudio)
if (!is.null(dev.list())) dev.off()
# Reset options
options(default = TRUE)
#options(prompt = FALSE)
options(ask = FALSE)
# Clean history (optional)
cat("", file = ".Rhistory")
#####################################################

output.file = "distributions_priors.txt"
# Differents quantiles
q <- c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)

sink(file = output.file, split = TRUE)

# POUR times

cat("\n")
cat("TIMES", "\n")

# distribution Unif
cat("d <- runif(1000000,0,20000)", "\n") # Added
d <- runif(1000000, 0, 20000)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(1),log(20000))", "\n") # Added
d <- runif(1000000, log(100000), log(10000000))
d <- exp(d)
# plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(2000),log(20000))", "\n") # Added
d <- runif(1000000, log(2000), log(20000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# POUR DB
cat("\n")
cat("DBN", "\n")

# distribution Unif
cat("d <- runif(1000000,2,200)", "\n") # Added
d <- runif(1000000, 0, 200)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(2),log(200))", "\n") # Added
d <- runif(1000000, log(2), log(200))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# POUR NBN
cat("\n")
cat("NBN", "\n")

# distribution Unif
cat("d <- runif(1000000,2,1000)", "\n") # Added
d <- runif(1000000, 2, 1000)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(2),log(1000))", "\n") # Added
d <- runif(1000000, log(2), log(1000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(20),log(2000))", "\n") # Added
d <- runif(1000000, log(20), log(2000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(10),log(1000))", "\n") # Added
d <- runif(1000000, log(10), log(1000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# POUR NA
cat("\n")
cat("NA", "\n")

# distribution Unif
cat("d <- runif(1000000,10000,1000000)", "\n") # Added
d <- runif(1000000, 10000, 1000000)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(10000),log(1000000))", "\n") # Added
d <- runif(1000000, log(10000), log(1000000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogUnif
cat("d <- runif(1000000,log(10000),log(10000000))", "\n") # Added
d <- runif(1000000, log(10000), log(10000000))
d <- exp(d)
plot(density(d))
summary(d)
quantile(d, probs = q)







# POUR NBN/DBD (bottleneck intensity WITH LU)
cat("\n")
cat("NBN/DBD (bottleneck intensity)", "\n")

# distribution LogUnif
cat("NBN <- runif(1000000,log(10),log(1000))", "\n") 

NBN <- runif(1000000, log(10), log(1000))
NBN <- exp(NBN)
plot(density(NBN))
summary(NBN)
quantile(NBN, probs = q)

cat("\n")
# Distribution LogUnif
cat("DBN <- runif(1000000,log(2),log(200))", "\n") 
DBN <- runif(1000000, log(2), log(200))
DBN <- exp(DBN)
plot(density(DBN))
summary(DBN)
quantile(DBN, probs = q)

cat("\n")
cat("NBN/DBD (bottleneck intensity)", "\n")
NBN_DBN <- NBN/DBN
plot(density(NBN_DBN))
summary(NBN_DBN)
quantile(NBN_DBN, probs = q)

# POUR NBN/DBD (bottleneck intensity WITH UNI)
cat("\n")
cat("NBN/DBD (bottleneck intensity)", "\n")

# distribution UNIF
cat("NBN <- runif(1000000,2),1000)", "\n") 

NBN <- runif(1000000, 2, 1000)
plot(density(NBN))
summary(NBN)
quantile(NBN, probs = q)

cat("\n")
# Distribution UNIF
cat("DBN <- runif(1000000,0,200)", "\n") 
DBN <- runif(1000000, 0, 200)
plot(density(DBN))
summary(DBN)
quantile(DBN, probs = q)

cat("\n")
cat("NBN/DBD (bottleneck intensity)", "\n")
NBN_DBN <- NBN/DBN
plot(density(NBN_DBN))
summary(NBN_DBN)
quantile(NBN_DBN, probs = q)

sink()














############ OLD #################
############################################

# distribution LogU
d <- runif(1000000,log(73),log(500))
d<-exp(d)
x11()
plot(density(d))
summary(d)
quantile(d, probs = q)

# distribution LogU
d <- runif(1000000,log(0.05),log(0.95))
d<-exp(d)
x11()
plot(density(d))
q <-c(0.01,0.05,0.10,0.25,0.5,0.75,0.90,0.95,0.99)
summary(d)
quantile(d, probs = q)


# TRUC BIZARE
# distribution Unif
d1 <- runif(100,0,3000)
# distribution LogU
d2 <- runif(100,log(73),log(500))
d2<-exp(d2)

d3=d2 if {d2<d1}


