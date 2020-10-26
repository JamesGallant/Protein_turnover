d <- read.delim("halflifeperrep.txt")
require(dplyr)
require(ggplot2)
require(limma)

d.CDC1551 <- data.frame(row.names = d$Protein.IDs,
                        CDC1551.count = d$CDC1551.count.z,
                        CDC1551.RSQ = d$CDC1551.rsq.z,
                        CDC1551.halflife = d$CDC1551.halflife.z)

d.dPPE38 <- data.frame(row.names = d$Protein.IDs,
                       dPPE38.count = d$dPPE38.count.z,
                       dPPE38.RSQ = d$dPPE38.rsq.z,
                       dPPE38.halflife = d$dPPE38.halflife.z)

d.Comp <- data.frame(row.names = d$Protein.IDs,
                     Comp.count = d$Comp.count.z,
                     Comp.RSQ = d$Comp.rsq.z,
                     Comp.halflife = d$Comp.halflife.z)

d.Uninfected <- data.frame(row.names = d$Protein.IDs,
                           Uninfected.count = d$Uninfected.count.z,
                           Uninfected.RSQ = d$Uninfected.rsq.z,
                           Uninfected.halflife = d$Uninfected.halflife.z)

d.CDC1551 <- d.CDC1551[which(d.CDC1551$CDC1551.count > 2 & d.CDC1551$CDC1551.RSQ > 0.85),] 
d.dPPE38 <- d.dPPE38[which(d.dPPE38$dPPE38.count > 2 & d.dPPE38$dPPE38.RSQ > 0.85),] 
d.Comp <- d.Comp[which(d.Comp$Comp.count > 2 & d.Comp$Comp.RSQ > 0.85),]
d.Uninfected <- d.Uninfected[which(d.Uninfected$Uninfected.count > 2 & d.Uninfected$Uninfected.RSQ > 0.85),] 

d.t1 <- merge(d.CDC1551,d.dPPE38, 
              by = 0, all = T)
rownames(d.t1) <- d.t1$Row.names; d.t1$Row.names <- NULL

d.t2 <- merge(d.t1, d.Comp, by = 'row.names', all = T)
rownames(d.t2) <- d.t2$Row.names; d.t2$Row.names <- NULL

d.t3 <- merge(d.t2, d.Uninfected, by = 'row.names', all = T)
rownames(d.t3) <- d.t3$Row.names

d <- read.delim("halflifeperrep.txt")

#####
#reps
#Zero included 

d.halflife <- data.frame(row.names = d$Protein.IDs,
                         CDC1551.1.halflife = d$CDC1551.1.halflife.z,
                         CDC1551.2.halflife = d$CDC1551.2.halflife.z,
                         dPPE38.1.halflife = d$dPPE38.1.halflife.z,
                         dPPE38.2.halflife = d$dPPE38.2.halflife.z,
                         Comp.1.halflife = d$Comp.1.halflife.z,
                         Comp.2.halflife = d$Comp.2.halflife.z,
                         Uninfected.1.halflife = d$Uninfected.1.halflife.z,
                         Uninfected.2.halflife = d$Uninfected.2.halflife.z)

d.zero <- data.frame(ID=d$Protein.IDs,
                     CDC1551.1.halflife=d$CDC1551.1.halflife.z,
                     CDC1551.2.halflife=d$CDC1551.2.halflife.z,
                     CDC1551.av.halflife=rowMeans(d.halflife[,1:2]),
                     CDC1551.1.count=d$CDC1551.1.count.z,
                     CDC1551.2.count=d$CDC1551.2.count.z,
                     CDC1551.1.rsq=d$CDC1551.1.rsq.z,
                     CDC1551.2.rsq=d$CDC1551.2.rsq.z,
                     dPPE38.1.halflife=d$dPPE38.1.halflife.z,
                     dPPE38.2.halflife=d$dPPE38.2.halflife.z,
                     dPPE38.av.halflife=rowMeans(d.halflife[,3:4]),
                     dPPE38.1.count=d$dPPE38.1.count.z,
                     dPPE38.2.count=d$dPPE38.2.count.z,
                     dPPE38.1.rsq=d$dPPE38.1.rsq.z,
                     dPPE38.2.rsq=d$dPPE38.2.rsq.z,
                     Comp.1.halflife=d$Comp.1.halflife.z,
                     Comp.2.halflife=d$Comp.2.halflife.z,
                     Comp.av.halflife=rowMeans(d.halflife[,5:6]),
                     Comp.1.count=d$Comp.1.count.z,
                     Comp.2.count=d$Comp.2.count.z,
                     Comp.1.rsq=d$Comp.1.rsq.z,
                     Comp.2.rsq=d$Comp.2.rsq.z,
                     Uninfected.1.halflife=d$Uninfected.1.halflife.z,
                     Uninfected.2.halflife=d$Uninfected.2.halflife.z,
                     Uninfected.av.halflife=rowMeans(d.halflife[,7:8]),
                     Uninfected.1.count=d$Uninfected.1.count.z,
                     Uninfected.2.count=d$Uninfected.2.count.z,
                     Uninfected.1.rsq=d$Uninfected.1.rsq.z,
                     Uninfected.2.rsq=d$Uninfected.2.rsq.z)

d.CDC1551.1 <- data.frame(row.names = d.zero$ID,
                          CDC1551.1.count = d.zero$CDC1551.1.count,
                          CDC1551.1.rsq = d.zero$CDC1551.1.rsq,
                          CDC1551.1.halflife = d.zero$CDC1551.1.halflife)

d.CDC1551.2 <- data.frame(row.names = d.zero$ID,
                          CDC1551.2.count = d.zero$CDC1551.2.count,
                          CDC1551.2.rsq = d.zero$CDC1551.2.rsq,
                          CDC1551.2.halflife = d.zero$CDC1551.2.halflife)

d.dPPE38.1 <- data.frame(row.names = d.zero$ID,
                         dPPE38.1.count = d.zero$dPPE38.1.count,
                         dPPE38.1.rsq = d.zero$dPPE38.1.rsq,
                         dPPE38.1.halflife = d.zero$dPPE38.1.halflife)

d.dPPE38.2 <- data.frame(row.names = d.zero$ID,
                         dPPE38.2.count = d.zero$dPPE38.2.count,
                         dPPE38.2.rsq = d.zero$dPPE38.2.rsq,
                         dPPE38.2.halflife = d.zero$dPPE38.2.halflife)

d.Comp.1 <- data.frame(row.names = d.zero$ID,
                       Comp.1.count = d.zero$Comp.1.count,
                       Comp.1.rsq = d.zero$Comp.1.rsq,
                       Comp.1.halflife = d.zero$Comp.1.halflife)

d.Comp.2 <- data.frame(row.names = d.zero$ID,
                       Comp.2.count = d.zero$Comp.2.count,
                       Comp.2.rsq = d.zero$Comp.2.rsq,
                       Comp.2.halflife = d.zero$Comp.2.halflife)

d.Uninfected.1 <- data.frame(row.names = d.zero$ID,
                             Uninfected.1.count = d.zero$Uninfected.1.count,
                             Uninfected.1.rsq = d.zero$Uninfected.1.rsq,
                             Uninfected.1.halflife = d.zero$Uninfected.1.halflife)

d.Uninfected.2 <- data.frame(row.names = d.zero$ID,
                             Uninfected.2.count = d.zero$Uninfected.2.count,
                             Uninfected.2.rsq = d.zero$Uninfected.2.rsq,
                             Uninfected.2.halflife = d.zero$Uninfected.2.halflife)





d.CDC1551.1 <- d.CDC1551.1[which(d.CDC1551.1$CDC1551.1.count > 2 & d.CDC1551.1$CDC1551.1.rsq > 0.85),]
d.CDC1551.2 <- d.CDC1551.2[which(d.CDC1551.2$CDC1551.2.count > 2 & d.CDC1551.2$CDC1551.2.rsq > 0.85),]

d.dPPE38.1 <- d.dPPE38.1[which(d.dPPE38.1$dPPE38.1.count > 2 & d.dPPE38.1$dPPE38.1.rsq > 0.85),]
d.dPPE38.2 <- d.dPPE38.2[which(d.dPPE38.2$dPPE38.2.count > 2 & d.dPPE38.2$dPPE38.2.rsq > 0.85),]

d.Comp.1 <- d.Comp.1[which(d.Comp.1$Comp.1.count > 2 & d.Comp.1$Comp.1.rsq > 0.85),]
d.Comp.2 <- d.Comp.2[which(d.Comp.2$Comp.2.count > 2 & d.Comp.2$Comp.2.rsq > 0.85),]

d.Uninfected.1 <- d.Uninfected.1[which(d.Uninfected.1$Uninfected.1.count > 2 & d.Uninfected.1$Uninfected.1.rsq > 0.85),]
d.Uninfected.2 <- d.Uninfected.2[which(d.Uninfected.2$Uninfected.2.count > 2 & d.Uninfected.2$Uninfected.2.rsq > 0.85),]

d.t1 <- merge(d.CDC1551.1, d.CDC1551.2, by = 0, all = T)
rownames(d.t1) <- d.t1$Row.names; d.t1$Row.names <- NULL

d.t2 <- merge(d.t1, d.dPPE38.1, by = 'row.names', all = T )
rownames(d.t2) <- d.t2$Row.names; d.t2$Row.names <- NULL

d.t3 <- merge(d.t2, d.dPPE38.2, by = 'row.names', all = T )
rownames(d.t3) <- d.t3$Row.names; d.t3$Row.names <- NULL

d.t4 <- merge(d.t3, d.Comp.1, by = 'row.names', all = T )
rownames(d.t4) <- d.t4$Row.names; d.t4$Row.names <- NULL

d.t5 <- merge(d.t4, d.Comp.2, by = 'row.names', all = T )
rownames(d.t5) <- d.t5$Row.names; d.t5$Row.names <- NULL

d.t6 <- merge(d.t5, d.Uninfected.1, by = 'row.names', all = T )
rownames(d.t6) <- d.t6$Row.names; d.t6$Row.names <- NULL

d.t7 <- merge(d.t6, d.Uninfected.2, by = 'row.names', all = T )
rownames(d.t7) <- d.t7$Row.names; d.t7$Row.names <- NULL

d.zero.filtered

d.t1 <- merge(d.CDC1551,d.dPPE38, 
              by = 0, all = T)
rownames(d.t1) <- d.t1$Row.names; d.t1$Row.names <- NULL

d.t2 <- merge(d.t1, d.Comp, by = 'row.names', all = T)
rownames(d.t2) <- d.t2$Row.names; d.t2$Row.names <- NULL

d.t3 <- merge(d.t2, d.Uninfected, by = 'row.names', all = T)
rownames(d.t3) <- d.t3$Row.names

d.zero.halflife <- data.frame(row.names = rownames(d.t7),
                              CDC1551.1.halflife = d.t7$CDC1551.1.halflife,
                              CDC1551.2.halflife = d.t7$CDC1551.2.halflife,
                              dPPE38.1.halflife = d.t7$dPPE38.1.halflife,
                              dPPE38.2.halflife = d.t7$dPPE38.2.halflife,
                              Comp.1.halflife = d.t7$Comp.1.halflife,
                              Comp.2.halflife = d.t7$Comp.2.halflife,
                              Uninfected.1.halflife = d.t7$Uninfected.1.halflife,
                              Uninfected.2.halflife = d.t7$Uninfected.2.halflife)


#average halflife

d.zero.halflife.av <- data.frame(row.names = rownames(d.halflife),
                                 CDC1551.halflife = rowMeans(d.halflife[,1:2]),
                                 dPPE38.halflife = rowMeans(d.halflife[,3:4]),
                                 Comp.halflife = rowMeans(d.halflife[,5:6]),
                                 Uninfected.halflife = rowMeans(d.halflife[,7:8]))


lm.R2.CDC1551.dPPE38 <- lm(CDC1551.halflife ~ dPPE38.halflife, data = d.zero.halflife.av)
R2.CDC1551.dPPE38 <- summary(lm.R2.CDC1551.dPPE38)$r.squared

ggplot(d.zero.halflife.av, aes(log(CDC1551.halflife), log(dPPE38.halflife))) +
  geom_point(pch=21, size=5, fill="#f0a33e", alpha = 0.8, color = 'black') +
  geom_text(aes(x = 2, y = 5.5, label = round(R2.CDC1551.dPPE38, digits = 2))) + 
  theme_classic(base_size = 14) +
  xlab("Log(CDC1551 half-life)") + ylab("Log(dPPE38 half-life)")

lm.R2.CDC1551.Comp <- lm(CDC1551.halflife ~ Comp.halflife, data = d.zero.halflife.av)
R2.CDC1551.Comp <- summary(lm.R2.CDC1551.Comp)$r.squared

ggplot(d.zero.halflife.av, aes(log(CDC1551.halflife), log(Comp.halflife))) +
  geom_point(pch=21, size=5, fill="#f0a33e", alpha = 0.8, color = 'black') +
  geom_text(aes(x = 2, y = 5.5, label = round(R2.CDC1551.Comp, digits = 2))) + 
  theme_classic(base_size = 14) +
  xlab("Log(CDC1551 half-life)") + ylab("Log(Comp half-life)")

lm.R2.CDC1551.Uninfected <- lm(CDC1551.halflife ~ Uninfected.halflife, data = d.zero.halflife.av)
R2.CDC1551.Uninfected <- summary(lm.R2.CDC1551.Uninfected)$r.squared

ggplot(d.zero.halflife.av, aes(log(CDC1551.halflife), log(Uninfected.halflife))) +
  geom_point(pch=21, size=5, fill="#f0a33e", alpha = 0.8, color = 'black') +
  geom_text(aes(x = 2, y = 5.5, label = round(R2.CDC1551.Uninfected, digits = 3))) + 
  theme_classic(base_size = 14) +
  xlab("Log(CDC1551 half-life)") + ylab("Log(Uninfected half-life)")


lm.R2.dPPE38.Uninfected <- lm(dPPE38.halflife ~ Uninfected.halflife, data = d.zero.halflife.av)
R2.dPPE38.Uninfected <- summary(lm.R2.dPPE38.Uninfected)$r.squared

ggplot(d.zero.halflife.av, aes(log(dPPE38.halflife), log(Uninfected.halflife))) +
  geom_point(pch=21, size=5, fill="#f0a33e", alpha = 0.8, color = 'black') +
  geom_text(aes(x = 2, y = 5.5, label = round(R2.dPPE38.Uninfected, digits = 2))) + 
  theme_classic(base_size = 14) +
  xlab("Log(dPPE38 half-life)") + ylab("Log(Uninfected half-life)")

lm.R2.dPPE38.Comp <- lm(dPPE38.halflife ~ Comp.halflife, data = d.zero.halflife.av)
R2.dPPE38.Comp <- summary(lm.R2.dPPE38.Comp)$r.squared

ggplot(d.zero.halflife.av, aes(log(dPPE38.halflife), log(Comp.halflife))) +
  geom_point(pch=21, size=5, fill="#f0a33e", alpha = 0.8, color = 'black') +
  geom_text(aes(x = 2, y = 5.5, label = round(R2.dPPE38.Comp, digits = 2))) + 
  theme_classic(base_size = 14) +
  xlab("Log(dPPE38 half-life)") + ylab("Log(Comp half-life)")


####3
#heatmap
require(gplots)

my_palette <- colorRampPalette(brewer.pal(n=7, name= "RdYlBu"))



dat2 <- na.omit(d.zero.halflife.av)
dat2 <- log2(dat2)
breaks <- seq(min(dat2, na.rm = T), max(dat2, na.rm = T), length.out = 21)

dist_no_na <- function(dat2) {
  edist <- dist(dat2)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

heatmap.2(as.matrix(dat2), trace="none", na.color = "#595757", scale="none", 
          col = my_palette, breaks=breaks, distfun = dist_no_na, labRow = F, cexCol = 1, )

##
#limma

dat.CDC1551vsdPPE38 <- read.delim("avHalflfie_CDCvsdPPE.txt", row.names = 1)
dat <- log2(dat.CDC1551vsdPPE38)

meta.dat <- data.frame(Group = c(rep("CDC1551", 2), rep("dPPE38", 2)))

rownames(meta.dat) <- colnames(dat)


#make df for limma
f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(dat, design)

cont.matrix <- makeContrasts(CDC1551vsdPPE38="CDC1551-dPPE38", 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)


volc.plot.1 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.1 = mutate(volc.plot.1, 
                     sig = ifelse(volc.plot.1$EffectSize > 0.01 & volc.plot.1$p.adj < 0.05, "upregulated",
                                  ifelse(volc.plot.1$EffectSize < -0.01 & volc.plot.1$p.adj < 0.05, "downregulated", "non-significant")))


ggplot(volc.plot.1, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_point(aes(fill=sig), alpha = 0.8, size = 5, colour="black", pch=21) +
  scale_fill_manual(values=c("non-significant" = "#595757",
                             "downregulated" = "#5757f9", 
                             "upregulated" = "#f0a33e")) +
  theme(legend.position = "top", legend.title = element_blank())

#CDCvsComp
dat.CDC1551vsComp <- read.delim("avHalflfie_CDCvsComp.txt", row.names = 1)
dat <- log2(dat.CDC1551vsComp)

meta.dat <- data.frame(Group = c(rep("CDC1551", 2), rep("Comp", 2)))

rownames(meta.dat) <- colnames(dat)



f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(dat, design)

cont.matrix <- makeContrasts(CDC1551vsComp="CDC1551-Comp", 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)


volc.plot.2 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.2 = mutate(volc.plot.2, 
                     sig = ifelse(volc.plot.2$EffectSize > 0.01 & volc.plot.2$p.adj < 0.05, "upregulated",
                                  ifelse(volc.plot.2$EffectSize < -0.01 & volc.plot.2$p.adj < 0.05, "downregulated", "non-significant")))


ggplot(volc.plot.2, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_point(aes(fill=sig), alpha = 0.8, size = 5, colour="black", pch=21) +
  scale_fill_manual(values=c("non-significant" = "#595757",
                             "downregulated" = "#5757f9", 
                             "upregulated" = "#f0a33e")) +
  theme(legend.position = "top", legend.title = element_blank())

#CDCvsuninf
dat.CDC1551vsuninf <- read.delim("avHalflfie_CDCvsuninf.txt", row.names = 1)
dat <- log2(dat.CDC1551vsuninf)

meta.dat <- data.frame(Group = c(rep("CDC1551", 2), rep("uninf", 2)))

rownames(meta.dat) <- colnames(dat)



f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(dat, design)

cont.matrix <- makeContrasts(CDC1551vsuninf="CDC1551-uninf", 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)


volc.plot.3 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.3 = mutate(volc.plot.3, 
                     sig = ifelse(volc.plot.3$EffectSize > 0.01 & volc.plot.3$p.adj < 0.05, "upregulated",
                                  ifelse(volc.plot.3$EffectSize < -0.01 & volc.plot.3$p.adj < 0.05, "downregulated", "non-significant")))


ggplot(volc.plot.3, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_point(aes(fill=sig), alpha = 0.8, size = 5, colour="black", pch=21) +
  scale_fill_manual(values=c("non-significant" = "#595757",
                             "downregulated" = "#5757f9", 
                             "upregulated" = "#f0a33e")) +
  theme(legend.position = "top", legend.title = element_blank())

#dPPE38vsuninf
dat.dPPE38vsuninf <- read.delim("avHalflife_dPPE38vsuninf.txt", row.names = 1)
dat <- log2(dat.dPPE38vsuninf)

meta.dat <- data.frame(Group = c(rep("dPPE38", 2), rep("uninf", 2)))

rownames(meta.dat) <- colnames(dat)



f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(dat, design)

cont.matrix <- makeContrasts(dPPE38vsuninf="dPPE38-uninf", 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)


volc.plot.4 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.4 = mutate(volc.plot.4, 
                     sig = ifelse(volc.plot.4$EffectSize > 0.01 & volc.plot.4$p.adj < 0.05, "upregulated",
                                  ifelse(volc.plot.4$EffectSize < -0.01 & volc.plot.4$p.adj < 0.05, "downregulated", "non-significant")))


ggplot(volc.plot.4, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_point(aes(fill=sig), alpha = 0.8, size = 5, colour="black", pch=21) +
  scale_fill_manual(values=c("non-significant" = "#595757",
                             "downregulated" = "#5757f9", 
                             "upregulated" = "#f0a33e")) +
  theme(legend.position = "top", legend.title = element_blank())


#Comp vs uninf
dat.Compvsuninf <- read.delim("avHalflife_Compvsuninf.txt", row.names = 1)
dat <- log2(dat.Compvsuninf)

meta.dat <- data.frame(Group = c(rep("Comp", 2), rep("uninf", 2)))

rownames(meta.dat) <- colnames(dat)



f.df <- factor(meta.dat$Group)
design <- model.matrix(~0+f.df)
colnames(design) <- levels(f.df)
#fit model
fit <- lmFit(dat, design)

cont.matrix <- makeContrasts(Compvsuninf="Comp-uninf", 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

#this is the final results array
fit2 <- eBayes(fit2)


volc.plot.5 <- data.frame(ID = names(fit2$coefficients[,1]),
                          p = fit2$p.value[,1],
                          p.adj = p.adjust(fit2$p.value[,1], "BH"),
                          EffectSize = fit2$coefficients[,1],
                          threshold = as.factor(-log10(p.adjust(fit2$p.value[,1], "BH")) <= -log10(0.01)),
                          comparison = 1)

volc.plot.5 = mutate(volc.plot.5, 
                     sig = ifelse(volc.plot.5$EffectSize > 0.01 & volc.plot.5$p.adj < 0.05, "upregulated",
                                  ifelse(volc.plot.5$EffectSize < -0.01 & volc.plot.5$p.adj < 0.05, "downregulated", "non-significant")))


ggplot(volc.plot.5, aes(x=EffectSize, y=-log10(p))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_point(aes(fill=sig), alpha = 0.8, size = 5, colour="black", pch=21) +
  scale_fill_manual(values=c("non-significant" = "#595757",
                             "downregulated" = "#5757f9", 
                             "upregulated" = "#f0a33e")) +
  theme(legend.position = "top", legend.title = element_blank())



write.table(volc.plot.1, file = "stats_CDCvsdPPE38.txt", sep = "\t", row.names = F)
write.table(volc.plot.2, file = "stats_CDCvsComp.txt", sep = "\t", row.names = F)
write.table(volc.plot.3, file = "stats_CDCvsuninfected.txt", sep = "\t", row.names = F)
write.table(volc.plot.4, file = "stats_dPPE38vsuninfected.txt", sep = "\t", row.names = F)
write.table(volc.plot.5, file = "stats_Compvsuninfected.txt", sep = "\t", row.names = F)
#without zero
d.noZ <- data.frame(ID=d$Protein.IDs,
                    CDC1551.1.halflife=d$CDC1551.1.halflife,
                    CDC1551.2.halflife=d$CDC1551.2.halflife,
                    CDC1551.1.count=d$CDC1551.1.count,
                    CDC1551.2.count=d$CDC1551.2.count,
                    CDC1551.1.rsq=d$CDC1551.1.rsq,
                    CDC1551.2.rsq=d$CDC1551.2.rsq,
                    dPPE38.1.halflife=d$dPPE38.1.halflife,
                    dPPE38.2.halflife=d$dPPE38.2.halflife,
                    dPPE38.1.count=d$dPPE38.1.count,
                    dPPE38.2.count=d$dPPE38.2.count,
                    dPPE38.1.rsq=d$dPPE38.1.rsq,
                    dPPE38.2.rsq=d$dPPE38.2.rsq,
                    Comp.1.halflife=d$Comp.1.halflife,
                    Comp.2.halflife=d$Comp.2.halflife,
                    Comp.1.count=d$Comp.1.count,
                    Comp.2.count=d$Comp.2.count,
                    Comp.1.rsq=d$Comp.1.rsq,
                    Comp.2.rsq=d$Comp.2.rsq,
                    Uninfected.1.halflife=d$Uninfected.1.halflife,
                    Uninfected.2.halflife=d$Uninfected.2.halflife,
                    Uninfected.1.count=d$Uninfected.1.count,
                    Uninfected.2.count=d$Uninfected.2.count,
                    Uninfected.1.rsq=d$Uninfected.1.rsq,
                    Uninfected.2.rsq=d$Uninfected.2.rsq)



####3
#GO

d <- read.delim("GO clusters.txt")

ggplot(d, aes(reorder(d$Description, Ratio), d$Ratio, fill = d$cluster)) +
  geom_bar(stat = "identity", color = "black", lwd=1.3, width = 0.6, position = "dodge" ) +
  scale_fill_manual(values = c("#5757f9", "#f0a33e")) +
  coord_flip() +
  theme_classic(base_size = 14) +
  xlab(NULL) + ylab("Enrichment ratio") + labs(fill = NULL) +
  theme(legend.position = "top")


