
library(ggplot2)
library(ggthemes)
library(ggrepel)


data1 <- read.csv("Lynch_data_with_newNe.csv", header = TRUE)

ggplot(data1, aes(x=log10(GS), y=log10(mu), colour=colour)) + geom_point() +
  scale_color_manual(values=c("black","orangered3","green4","dodgerblue3")) +
  geom_text_repel(aes(log10(GS), log10(mu), label=num))

ggsave("Fig2a.GS_vs_mu.pdf", width=8.8, height=7.2)

ggplot(data1, aes(x=log10(newNe), y=log10(mu), colour=colour)) + geom_point() +
  scale_color_manual(values=c("black","orangered3","green4","dodgerblue3")) +
  geom_text_repel(aes(log10(newNe), log10(mu), label=num))

ggsave("Fig2b.newNe_vs_mu.pdf", width=8.8, height=7.2)

ggplot(data1, aes(x=log10(newNe), y=log10(muperCDS), colour=colour)) + geom_point() +
  scale_color_manual(values=c("black","orangered3","green4","dodgerblue3")) +
  geom_text_repel(aes(log10(newNe), log10(muperCDS), label=num))

ggsave("Fig2c.newNe_vs_muperCDS.pdf", width=8.8, height=7.2)


ggplot(data1, aes(x=log10(GS), y=log10(muperCDS), colour=colour)) + geom_point() +
  scale_color_manual(values=c("black","orangered3","green4","dodgerblue3")) +
  geom_text_repel(aes(log10(GS), log10(muperCDS), label=num))

ggsave("Fig2d.GS_vs_muperCDS.pdf", width=8.8, height=7.2)



data1 <- read.csv("Lynch_data_with_newNe_and_muCDS.csv", header = TRUE)

ggplot(data1, aes(x=log10(GS), y=log10(CDSmu), colour=colour)) + geom_point() +
  scale_color_manual(values=c("black","orangered3","green4","dodgerblue3")) +
  geom_text_repel(aes(log10(GS), log10(CDSmu), label=num))

ggsave("Fig2e.GS_vs_CDSmu.pdf", width=8.8, height=7.2)


dev.off()

