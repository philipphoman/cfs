#-----------------------------------------------------------------------
# Reversal learning CFS
#
# This is the actual script, producing tables and figures
#-----------------------------------------------------------------------
# PH, 3/23/18
#-----------------------------------------------------------------------
source("../src/cfs_load.R")
# Create file marker
system("touch ../output/R/cfs_do.Rout")
#-----------------------------------------------------------------------
# Plot group responses and regressions 
#-----------------------------------------------------------------------
dfm <- df %>% group_by(manipulation, stage, spider, ctrial) %>%
  dplyr::summarize(mean=mean(logdcm2, na.rm=TRUE),
                   sd=sd(logdcm2, na.rm=TRUE),
                   n=sum(!is.na(logdcm2)),
                   se=sd/sqrt(n),
                   ci=se * 1.96)

# response traces
dfm$manipulation <- as.character(dfm$manipulation)
dfm$manipulation[dfm$manipulation=="noCFS"] <- "no-CFS"
p1 <- ggplot(dfm, aes(x=ctrial, y=mean, group=spider)) +
  geom_line() +
  geom_point(size=3, aes(shape=spider, fill=spider)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.0) +
  geom_vline(xintercept=8.5, linetype="dashed") + 
  facet_wrap(~manipulation, ncol=1, scale="free") +
  #facet_wrap(~manipulation, ncol=1, scale="free") +
  theme_classic(base_size=23) +
  theme(
    axis.title=element_text(face="bold"),
    #axis.ticks.x=element_blank(),
    legend.title=element_blank(),
    #legend.position="right",
    legend.position=c(0.88, 1.0),
    strip.text.x=element_text(face="bold", hjust=0),
    strip.background = element_blank()
    ) +
  scale_x_continuous(breaks=c(4, 12), labels=c("Acquisition",
                                               "Reversal")) +
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c("white", "black"))+
  xlab("Trial") +
  ylab("SNA Amplitude") +
  ylim(0.1, 0.95)

# repeat response traces with meanperc == 0
dfm2 <- df %>% filter(meanperc==0) %>%
  group_by(manipulation, stage, spider, ctrial) %>%
  dplyr::summarize(mean=mean(logdcm2, na.rm=TRUE),
                   sd=sd(logdcm2, na.rm=TRUE),
                   n=sum(!is.na(logdcm2)),
                   se=sd/sqrt(n),
                   ci=se * 1.96)
p1_1 <- ggplot(dfm2, aes(x=ctrial, y=mean, group=spider)) +
  geom_line() +
  geom_point(size=3, aes(shape=spider, fill=spider)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.0) +
  geom_vline(xintercept=8.5, linetype="dashed") + 
  facet_wrap(~manipulation, ncol=1, scale="free") +
  #facet_wrap(~manipulation, ncol=1, scale="free") +
  theme_classic(base_size=23) +
  theme(
    axis.title=element_text(face="bold"),
    #axis.ticks.x=element_blank(),
    legend.title=element_blank(),
    legend.position="top",
    strip.text.x=element_text(face="bold", hjust=0),
    strip.background = element_blank()
    ) +
  scale_x_continuous(breaks=c(4, 12), labels=c("Acquisition",
                                               "Reversal")) +
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c("white", "black"))+
  xlab("Trial") +
  ylab("SNA Amplitude") +
  ylim(0.1, 0.95)
ggsave(plot=p1_1, filename="../output/figures/cfs_timecourses_mperc0.pdf",
       width=9, height=4)

# Plot reversal index
dfri <- df %>% filter(trial==1) %>%
  group_by(manipulation) %>%
  dplyr::summarize(mean=mean(revindex, na.rm=TRUE),
                   sd=sd(revindex, na.rm=TRUE),
                   n=sum(!is.na(revindex), na.rm=TRUE),
                   se=sd/sqrt(n),
                   ci=se * 1.96) %>% mutate(x="Reversal index")

p2 <- ggplot(dfri, aes(x=manipulation, y=mean)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.0) +
  facet_wrap(~x, ncol=1, scale="free") +
  theme_classic(base_size=23) +
  theme(
    axis.line.y=element_blank(),
    axis.line.x=element_blank(),
    axis.ticks=element_blank(),
    #axis.text=element_blank(),
    plot.title=element_blank(),
    axis.title=element_text(face="bold"),
    legend.title=element_blank(),
    #strip.text.x=element_blank(),
    strip.text.x=element_text(face="bold"),
    strip.background = element_blank()) +
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c("white", "black"))+
  xlab("") +
  ylab("") + 
  geom_hline(yintercept=0, linetype="dashed") 

dfl <- df %>% filter(trial==1&manipulation=="CFS") %>%
  dplyr::select(id, stait, stais, fsq, revindex) %>%
  rename(STAIT=stait, STAIS=stais, FSQ=fsq) %>%
  gather(key=scale, value=score, STAIT, STAIS, FSQ) 

# Create data frame with labels of correlation strength for plot
rsl <- by(dfl, list(dfl$scale), lm, formula=formula("revindex ~ score"))
rl <- data.frame(r=sapply(rsl, function(x) sqrt(summary(x)$r.squared)),
                 sign=sapply(rsl, function(x) sign(coef(summary(x))[2])),
                 p=sapply(rsl, function(x) coef(summary(x))[8]),
                 scale=c("FSQ", "STAIS", "STAIT"),
                 x=c(85, 45, 45))
rl$r <- rl$r * rl$sign
rl$rstr <- paste("r=", round(rl$r, 2), starsfromp(rl$p), sep="")

#dfl$r <- rep(rs[c(3, 2, 1)], each=nrow(dfl)/3)
#dfl$p <- rep(ps[c(3, 2, 1)], each=nrow(dfl)/3)
                                      
p3 <- ggplot(dfl, aes(x=score, y=revindex)) +
  geom_point(size=4) +
  geom_smooth(method="lm") +
  facet_wrap(~scale, scale="free") +
  theme_gray(base_size=23) +
  theme(
    axis.title=element_text(face="bold"),
    strip.text=element_text(face="bold")
  ) +
  geom_text(data=rl, size=7, aes(x=x, y=5, label=rstr)) +
  xlab("Score") +
  ylab("Reversal index")

pg1 <- plot_grid(p1, plot_grid(p2, nrow=1, rel_heights=c(1, 1)),
                 ncol=2, rel_widths=c(2.0, 1), labels="AUTO",
                 scale=0.9, label_size=35)

pg2 <- plot_grid(p3, ncol=1, nrow=2, rel_heights=c(2.1, 1),
                 labels="C",
                 scale=0.9, label_size=35)

p <- plot_grid(pg1, pg2, ncol=1, axis="lgbt", align="hv")
ggsave(plot=p, filename="../output/figures/cfs_timecourses.pdf",
       width=12, height=13)


# alternative solutation as requested by DC
pg1 <- plot_grid(p1, plot_grid(p2, nrow=1, rel_heights=c(1, 1)),
                 ncol=2, rel_widths=c(2.0, 1), labels="AUTO",
                 scale=0.9, label_size=35)
pg2 <- plot_grid(p3, ncol=1, nrow=2, rel_heights=c(2.1, 1),
                 labels="C",
                 scale=0.94, label_size=35)
p <- plot_grid(pg1, pg2, ncol=1, axis="lgbt", align="hv",
               scale=1.0,
               rel_widths=c(1, 1, 2))
ggsave(plot=p, filename="../output/figures/cfs_timecourses.pdf",
       width=12, height=13)

#-----------------------------------------------------------------------
# n-tile plots 
#-----------------------------------------------------------------------
fm <- c(
  "revindex ~ meanpercsub",
  "revindex ~ stait + meanpercsub",
  "revindex ~ stait + hasproblem + meanpercsub",
  "revindex ~ stait + hasproblem + randomguess + meanpercsub"
  )

# Plot awareness even vs awareness odd
df$mp <- df$manipulation
snt <- df %>% filter(trial==1, manipulation=="CFS") %>%
  group_by(as.factor(ntile_e)) %>%
  dplyr::summarize(snte=sum(!is.na(ntile_e))) %>%
  mutate(q=c("UA", "T1", "T2", "T3")) %>%
  mutate(r1=cor.test(df$meanperceven[df$trial==1&df$mp=="CFS"],
                    df$meanpercodd[df$trial==1&df$mp=="CFS"])$estimate,
         p1=cor.test(df$meanperceven[df$trial==1&df$mp=="CFS"],
                    df$meanpercodd[df$trial==1&df$mp=="CFS"])$p.val,
         label1=paste("r=", round(r1, 2), starsfromp(p1), sep=""),
         r2=cor.test(df$meanperc[df$trial==1&df$mp=="CFS"],
                    df$revindex[df$trial==1&df$mp=="CFS"])$estimate,
         p2=cor.test(df$meanperc[df$trial==1&df$mp=="CFS"],
                    df$revindex[df$trial==1&df$mp=="CFS"])$p.val,
         label2=paste("r=", round(r2, 2), starsfromp(p2), sep=""))
         
write.csv(snt, "../output/tables/cfs_ntile_even.csv")

#scale_labels=paste(snt$q, "(N=", snt$snte, ")", sep="")
scale_labels=snt$q

xy_pairs <- list(
  c("meanperceven", "meanpercodd", "Awareness index (even)",
    "Awareness index (odd)", snt$label1[1], "ntile_e",
    scale_labels, "top", TRUE),
  c("meanperceven", "revindex_even", "Awareness index (even)",
    "Reversal index (even)", "", "ntile_e",
    scale_labels, "top", NULL),
  c("meanpercodd", "revindex_odd", "Awareness index (odd)",
    "Reversal index (odd)", "", "ntile_e",
    scale_labels, "top", NULL),
  c("meanperc", "revindex", "Awareness index",
    "Reversal index", snt$label2[1], "group",
    scale_labels, "", TRUE, TRUE))

plist <- lapply(xy_pairs,
                function(x) {
                  px <- ggplot(df %>%
                               filter(trial==1&manipulation=="CFS") %>%
                               mutate(ntile_e=as.factor(ntile_e)),
                               aes_string(x=x[1], y=x[2])) +
                    theme_classic(base_size=25) +
                    theme(
                      legend.title=element_blank(),
                      legend.position=x[11],
                      #legend.text=element_text(size=12),
                      axis.title=element_text(face="bold")
                    ) +
                    scale_color_hue(labels=x[7:10]) +
                    annotate("text", 0.8 *max(df[, x[1]]), max(df[, x[2]]),
                             label=x[5], size=7) +
                    xlab(x[3]) +
                    ylab(x[4])
                  if (!is.na(x[12])) {
                    px <- px + geom_smooth(method="lm")
                  } else {
                    px <- px + ggtitle("")
                  }
                  if (!is.na(x[13])) {
                    px <- px + geom_point(size=3, col="black")
                  } else {
                    px <- px + geom_point(size=3, aes_string(col=(x[6])))
                  }
                  })



fm <- c(
  "revindex ~ meanperc",
  "revindex ~ stait + meanperc",
  "revindex ~ stait + hasproblem + meanperc",
  "revindex ~ stait + hasproblem + randomguess + meanperc"
  )

lmfits <- lapply(fm, function(x)
  lm(x, data=df %>% filter(trial==1, manipulation=="CFS")))
# Create tables of regression results
coeftabslist <- lapply(lmfits, function(x) as.data.frame(coef(summary(x))))
coeftabs <- do.call("rbind", coeftabslist)
write.csv(coeftabs, "../output/tables/cfs_interceptmodels.csv")
    
coefs <- lapply(lmfits, function(x) coef(summary(x))[1, ])
coefsdf <- as.data.frame(do.call("rbind", coefs)) %>%
  mutate(model=c("Awareness index",
                 "Awareness index + STAIT",
                 "Awareness index + STAIT + tracking",
                 "Awareness index + STAIT + tracking + random guess"),
         modelshort=c("Model 1", "Model 2", "Model 3", "Model 4"),
         modelnum=c(1, 2, 3, 4),
         ci=`Std. Error` * 1.96,
         x=0.2)
pd <- position_dodge(width=0.9)
px <- ggplot(coefsdf[-4, ], aes(x=modelnum, y=Estimate)) +
  geom_point(size=3, position=pd) +
  geom_errorbar(aes(ymin=Estimate-ci, ymax=Estimate+ci, width=0),
                position=pd) +
  geom_hline(yintercept=0, linetype="dashed") + 
  theme_classic(base_size=25) +
  theme(
    axis.title=element_text(face="bold") 
    #legend.position=c(1, 0.5),
    #legend.position="right",
    #legend.title=element_blank(),
    #legend.text=element_text(size=20),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank()
  ) +
  #scale_color_manual(values=c("darkblue", "darkred", "darkgreen")) +
  xlab("Model") +
  ylab("Intercept with 95% CI")

# Extract legend
#leg <- g_legend(px)
#px <- px + theme(legend.position="")


p <- plot_grid(plist[[1]], plist[[2]], plist[[3]],
               plist[[4]], px, labels="AUTO",
               label_size=35, align="hv", axis="lgbt",
               ncol=3, scale=0.95, rel_widths=c(2, 2, 2, 3, 3))

ggsave(plot=p, filename="../output/figures/cfs_regressions.pdf",
       width=15.2, height=10)

# alternative solution as requested by DC
r1 <- plot_grid(plist[[1]], plist[[2]],
                plist[[3]], ncol=3, labels="AUTO",
                label_size=35, scale=0.9, align="hv", axis="lgbt")
r2 <- plot_grid(NULL, plist[[4]], px, NULL, ncol=4,
                rel_widths=c(1, 4, 4, 1),
                labels=c("", "D", "E", ""), label_size=35,
                scale=0.7, align="hv", axis="lgbt")
p <- plot_grid(r1, r2, nrow=2, rel_heights=c(1, 1))
ggsave(plot=p, filename="../output/figures/cfs_regressions-alt.pdf",
       width=15.2, height=10)
#-----------------------------------------------------------------------
# Single subject choice trace
#-----------------------------------------------------------------------
#dfl <- df %>% filter(id %in% c(1866, 7003)) %>%
dfl1 <- df %>% filter(id %in% c(1866)) %>%
  dplyr::select(id, trial, visresp, spidernum)  %>%
  mutate(spidernum=spidernum+1) %>%
  gather(key=tracetype, value=trace, spidernum, visresp)

dfl2 <- df %>% filter(id %in% c(1866)) %>%
  dplyr::select(id, trial, visresp, shockapplied)  %>%
  mutate(shockapplied=shockapplied+1) %>%
  gather(key=tracetype, value=trace, shockapplied, visresp)

dfl <- rbind(dfl1, dfl2) %>% mutate(comparison=rep(c("resp vs stim",
                                                     "resp vs shock"),
                                                   each=64))

p1 <- ggplot(dfl %>% filter(comparison=="resp vs stim"),
             aes(x=trial, y=visresp)) +
  geom_line(aes(x=trial, y=trace, col=tracetype), size=3) +
  facet_wrap(~comparison, scale="free", nrow=2) +
  scale_y_continuous(breaks=c(1, 2), labels=c("Flower", "Spider")) +
  scale_colour_manual(name="", values=c("red", "blue"),
                      labels=c("Stimulus trace",
                               "Choice trace")) +
  theme_classic(base_size=25) +
  theme(axis.title=element_text(face="bold"),
        strip.text=element_blank()) +
  xlab("Trial") +
  ylab("Choice")

png1 <- readPNG("../output/figures/cfs_binomsim-crop.png")
p0 <- rasterGrob(png1, interpolate=FALSE)

p <- plot_grid(plot_grid(p1, scale=0.8, labels=NULL),
               p0, labels="AUTO", label_size=35,
               rel_widths=c(3, 1), scale=0.9)
print(p)
ggsave(plot=p, filename="../output/figures/cfs_sub_binomsim.pdf",
       width=18, height=7)

#-----------------------------------------------------------------------
# Model comparison
#-----------------------------------------------------------------------
pdf("../output/figures/cfs_modelfits.pdf", paper="a4r", height=0, width=0)
par(mar=c(2, 1, 2, 1), oma=c(5, 5, 5, 5))
layout(matrix(c(1:16), nrow=4, byrow=TRUE))
summarize_models()
dev.off()




