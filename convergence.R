# convergence diagonostics

library(coda)
### load the files
chain1 <- read.coda(output.file =  "coda1.txt", index.file="codaIndex.txt", quiet=FALSE)
dim(chain1)  # 1000 * 25105  iteration * varible

chain2 <- read.coda(output.file =  "coda2.txt", index.file="codaIndex.txt", quiet=FALSE)

chainboth <- mcmc.list(chain1,chain2)

################   pick a subset of parameters
c1 <- chain1[,1:20]
c2 <- chain2[,1:20]
cc12 <- mcmc.list(c1,c2)
###
# Diagonostics
###
pdf(file = "trace.pdf")
plot(cc12, density=TRUE, smooth=FALSE)  # trace plots & density plots
dev.off()


densplot(cc12)  # density plots
autocorr.plot(cc12)  # autocorrelation plot
effectiveSize(cc12)  # effective sample size
gelman.diag(cc12)    # PSRF
gelman.plot(cc12)

summary(cc12)
