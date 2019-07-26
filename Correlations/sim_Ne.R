data1 <- read.csv("Lynch_data.csv", header = TRUE)

sp <- unlist(data1$sp)
num <- unlist(data1$num)
mu <- unlist(data1$mu)
pi <- unlist(data1$nmu)
hord <- unlist(as.integer(data1$hord))
Lynch_ne <- unlist(data1$ne)
GS <- unlist(data1$GS)
CDS <- unlist(data1$CDS)
muperCDS <- unlist(data1$muperCDS)
phyl <-unlist(data1$muperCDS)




myNe <- pi/(hord*mu - hord*pi*mu)

data1$newNe <- myNe

dev <- abs((myNe - Lynch_ne)/((myNe+Lynch_ne)/2))


#test whether Ne mu are independent

real.cor <- cor.test(myNe, mu, method="spearman")
real.rho <- real.cor$estimate

real.cor.cds <- cor.test(myNe, phyl, method="spearman")
real.rho.cds <- real.cor.cds$estimate

r.rholist  <- c()
r.rholist.cds <- c()

r.plist  <- c()
r.plist.cds <- c()


nsims <- 100000

for (i in c(1:nsims)) {

r.pi <- sample(pi)

r.Ne <- r.pi/(hord*mu - hord*r.pi*mu)

r.cor <-cor.test(r.Ne, mu, method = "spearman")
r.cor.cds <- cor.test(r.Ne, phyl, method = "spearman")

r.rho <- r.cor$estimate
r.rho.cds <-r.cor.cds$estimate

r.p <- r.cor$p.value
r.p.cds <- r.cor.cds$p.value

r.rholist <- c(r.rholist, r.rho)
r.rholist.cds <- c(r.rholist.cds, r.rho.cds)


r.plist <- c(r.plist, r.p)
r.plist.cds <- c(r.plist.cds, r.p.cds)
}

moreneg <- r.rholist[r.rholist <= real.rho]

numless <- length(moreneg)
P <- numless/nsims

moreneg.cds <- r.rholist.cds[r.rholist.cds <= real.rho.cds]

numless.cds <- length(moreneg.cds)
P.cds <- numless.cds/nsims
par(mfrow=c(1,2))
hist(r.rholist)
hist(r.rholist.cds)



numsig <- length(r.plist[r.plist <0.05])
numsig.cds <- length(r.plist.cds[r.plist.cds <0.05])