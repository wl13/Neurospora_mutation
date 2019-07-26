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
col <- unlist(as.character(data1$colour))
size = unlist(data1$size)




myNe <- pi/(hord*mu - hord*pi*mu)

data1$newNe <- myNe

dev <- abs((myNe - Lynch_ne)/((myNe+Lynch_ne)/2))

pdf("Fig1.pdf", height = 5, width = 15)

par(mfrow=c(1,3))
#pdf("1a_GenomeSizevMu.pdf")
par(pty="s")
plot(log10(GS), log10(mu), pch=0,  xlab = "", ylab="", cex = size*3)
text(x=log10(GS), y=log10(mu), labels = num, cex =0.9, col=col)
mtext("Log(Genome Size in Mb)", side=1, cex=1, line =3)
mtext("Log(Mutation rate per bp per generation)", side=2, line =3)
mtext("a", side=3, line = 1, adj =0)
#dev.off()


#pdf("1b_NevMu.pdf")

plot(log10(myNe), log10(mu), pch=0, xlab="", ylab = "", cex = size*3)
text(x=log10(myNe), y=log10(mu), labels = num, cex=0.9, col=col)
mtext(expression("Log"(italic("N")["e"])), side=1, cex=1, line = 3)
mtext("Log(Mutation rate per bp per generation)", side=2, line =3)
mtext("b", side=3, line = 1, adj =0)
#dev.off()

#pdf("1c_NevMupercds.pdf")

plot(log10(myNe), log10(muperCDS), pch=0,  xlab="", ylab = "", cex = size*3)
text(x=log10(myNe), y=log10(muperCDS),labels = num, cex=0.9, col=col)
mtext(expression("Log"(italic("N")["e"])), side=1, cex=1, line = 3)
mtext("Log(mutations in CDS per generation)", side=2, line =3)
mtext("c", side=3, line = 1, adj =0)



dev.off()

write.csv(data1, file="New_Data.csv")

