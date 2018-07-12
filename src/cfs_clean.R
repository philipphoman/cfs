#-----------------------------------------------------------------------
# Reversal learning CFS
#
# Data cleaning
#-----------------------------------------------------------------------
# PH, 3/25/18
#-----------------------------------------------------------------------
# read data frame
df <- read.csv("../data/cfs.csv")

# Logarithmize the outcome measure
df$logdcm2 <- log(df$dcm2 + 1)

# Resrict to right amount of trials and full data set
df  <- df %>% filter(trial<=32, !is.na(age), !is.na(gender),
                     !is.na(meanperc), ntrials!=0) %>%
  mutate(trialeven=ifelse(trial %% 2 == 0, TRUE, FALSE))

# Summarize mean shocks per subject
dfu <- df %>% group_by(id) %>%
  dplyr::summarize(meanur=mean(dcm3, na.rm=TRUE)) %>%
  left_join(df)
df <- dfu

# Summarize confidence per subject and correct/incorrect resp
#dfi <- df %>% group_by(visresp, id) %>% filter(!is.na(visresp)) %>%
#  filter(visresp==1) %>%
#  dplyr::summarize(meanconf_incorr=mean(confresp, na.rm=TRUE)) %>% 
#  dplyr::select(-visresp) %>%
#  left_join(df)
#df <- dfi
  #dplyr::select(id, visresp, confresp, meanconf_corrincorr) 

dfi <- df %>% group_by(id) %>% filter(!is.na(visresp)) %>%
  filter(visresp==2) %>%
  dplyr::summarize(meanconf_incorr=mean(confresp, na.rm=TRUE)) %>%
  right_join(df)
df <- dfi

dfc <- df %>% group_by(id) %>% filter(!is.na(visresp)) %>%
  filter(visresp==2) %>%
  dplyr::summarize(meanconf_corr=mean(confresp, na.rm=TRUE))  %>%
  right_join(df)
df <- dfc

dfm <- df %>% group_by(id) %>% filter(!is.na(visresp)) %>%
  dplyr::summarize(meanconf_new=mean(confresp, na.rm=TRUE))  %>%
  right_join(df)
df <- dfm

# Calculate reversal index
dfm <- df %>% group_by(stage, spider, id) %>%
    dplyr::summarize(meanstim=mean(dcm2, na.rm=TRUE)) %>%
    group_by(stage, id) %>%
    dplyr::summarize(diffstim=diff(meanstim, na.rm=TRUE)) %>%
    group_by(id) %>%
    dplyr::summarize(revindex=diff(diffstim, na.rm=TRUE))  %>%
    left_join(df)
df <- dfm

# Calculate reversal index for even/odd trials
dfm <- df %>% filter(trialeven==TRUE) %>%
  group_by(stage, spider, id) %>%
  dplyr::summarize(meanstim=mean(dcm2, na.rm=TRUE)) %>%
  group_by(stage, id) %>%
  dplyr::summarize(diffstim=diff(meanstim, na.rm=TRUE)) %>%
  group_by(id) %>%
  dplyr::summarize(revindex_even=diff(diffstim, na.rm=TRUE))  %>%
  left_join(df)
df <- dfm

dfm <- df %>% filter(trialeven==FALSE) %>%
  group_by(stage, spider, id) %>%
  dplyr::summarize(meanstim=mean(dcm2, na.rm=TRUE)) %>%
  group_by(stage, id) %>%
  dplyr::summarize(diffstim=diff(meanstim, na.rm=TRUE)) %>%
  group_by(id) %>%
  dplyr::summarize(revindex_odd=diff(diffstim, na.rm=TRUE))  %>%
  left_join(df)
df <- dfm

df <- preparedataset(df)

# Add mean perception index for even/odd trials
dfm <- df %>% group_by(trialeven, id) %>%
  dplyr::summarize(mp=mean(bleedlevel, na.rm=TRUE))  %>%
  spread(key=trialeven, value=mp) %>%
  rename(meanpercodd=`FALSE`,
         meanperceven=`TRUE`) %>%
  left_join(df) 
df <- dfm

# Add ntiles 
dfm <- df %>% filter(trial==1) %>%
  dplyr::select(id, manipulation, meanperc, meanperceven, meanpercodd) %>%
  gather(key=trialtype, value=mp, meanperc, meanperceven, meanpercodd) %>%
  arrange(manipulation, trialtype, id)


nt <- by(dfm, list(dfm$trialtype, dfm$manipulation),
         function(x) ntile(x$mp[x$mp>0], n=3)+1)
dfm$ntile <- 1
dfm$ntile[dfm$mp>0] <- unlist(nt)

dfm <- dfm %>%
  dplyr::select(-mp) %>%
  spread(key=trialtype, value=ntile) %>%
  rename(ntile=meanperc,
         ntile_e=meanperceven,
         ntile_o=meanpercodd) %>%
  dplyr::select(-manipulation) %>%
  left_join(df)
df <- dfm

#dfm$ntile[which(!is.na(dfm$ntile))] <- dfm$ntile[which(!is.na(dfm$ntile))] + 1
#dfm$ntile[which(is.na(dfm$ntile))] <- 1 

# Add "problem flag", i.e., likely breakthroughs
dfh <- createprobflag(df=df, ollim=8)
df  <- merge(df, subset(dfh, !is.na(id)), all=TRUE)
write.csv("../data/cfs.csv")

# Restrict to complete baseline cases
#df <- df %>% filter(!is.na(stait), !is.na(stais), !is.na(fsq))
#-----------------------------------------------------------------------

write.csv(df, "../data/cfs.csv")

