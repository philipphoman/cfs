# CFS preprocessing
#
# 11/9/16, PH

rm(list=ls())
library("data.table")
library("plyr")
setwd(".")

# load files
cs  <- read.csv("cfs_subjectdata.csv")
ce  <- read.csv("cfs_eprime.csv")
co  <- read.csv("cfs_orders.csv")
cd  <- read.csv("cfs_dcm3.csv")

# merge eprime with subject-data
ces  <- merge(ce,cs,all=T)

# also merge dcm with subject-data
cds  <- merge(cd,cs,all=T)

# merge both with orders
ceso  <- merge(ces,co,by=c("order", "trial"),all=T)
cdso  <- merge(cds,co,by=c("order", "trial"),all=T)

# merge ceso with dcm-data
c  <- merge(ceso,cdso,all=T)

# add missing age, gender, group if possible
c$age[c$id==8021] <- 19
c$gender[c$id==8021] <- "M" 
c$age[c$id==7382] <- 28 
c$gender[c$id==7382] <- "F" 
c$group[c$id==7382] <- "Unaware"
c$age[c$id==6534] <- 23 
c$gender[c$id==6534] <- "F" 
c$age[c$id==524] <- 28 
c$age[c$id==8688] <- 22 
c$gender[c$id==8688] <- "F" 
c$age[c$id==3763] <- 19 
c$gender[c$id==3763] <- "F" 
c$age[c$id==7791] <- 25 
c$gender[c$id==7791] <- "F" 
c$age[c$id==9247] <- 28 
c$gender[c$id==9247] <- "F" 
c$age[c$id==6793] <- 19 
c$gender[c$id==6793] <- "F" 
c$age[c$id==9076] <- 34 
c$gender[c$id==9076] <- "M" 
c$age[c$id==3019] <- 33 
c$gender[c$id==3019] <- "M" 
c$age[c$id==7271] <- 29 
c$gender[c$id==7271] <- "M" 
c$age[c$id==9247] <- 28 
c$gender[c$id==9247] <- "M" 
c$age[c$id==5765] <- 19 
c$gender[c$id==5765] <- "F" 



# set bleedlevel
c$bleedlevel <- NA
c$bleedlevel[c$visresp==2&c$confresp==3&!is.na(c$spider)] <- 1
c$bleedlevel[c$visresp==2&c$confresp==2&!is.na(c$spider)] <- 0.5 
c$bleedlevel[c$visresp==2&c$confresp==1&!is.na(c$spider)] <- 0
c$bleedlevel[c$visresp==1&c$confresp==3&!is.na(c$spider)] <- 0
c$bleedlevel[c$visresp==1&c$confresp==2&!is.na(c$spider)] <- 0
c$bleedlevel[c$visresp==1&c$confresp==1&!is.na(c$spider)] <- 0
#c$bleedlevel[c$group=="Aware"&!is.na(c$spider)] <- 1
#c$bleedlevel <- c$bleedlevel+1

# set shock expectancy
#c$shockprob <- factor(c$shockprob)
#c$shockprob[c$shockprob=="random"&c$spider=="Spider A"&c$stage=="Acquisition"] <- "random-high"
#c$shockprob[c$shockprob=="random"&c$spider=="Spider B"&c$stage=="Acquisition"] <- "random-low"
#c$shockprob[c$shockprob=="random"&c$spider=="Spider A"&c$stage=="Reversal"] <- "random-low"
#c$shockprob[c$shockprob=="random"&c$spider=="Spider B"&c$stage=="Reversal"] <- "random-high"

# create mean awareness
ca <- aggregate(c$bleedlevel,by=list(c$id),FUN=mean,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$meanperc <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)



# create median awareness
ca <- aggregate(c$bleedlevel,by=list(c$id),FUN=median,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$medianperc <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)

# create median and mean response 
ca <- aggregate(c$visresp,by=list(c$id),FUN=median,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$medianresp <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)
ca <- aggregate(c$visresp,by=list(c$id),FUN=mean,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$meanresp <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)


# create median and mean confidence 
ca <- aggregate(c$confresp,by=list(c$id),FUN=median,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$medianconf <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)
ca <- aggregate(c$confresp,by=list(c$id),FUN=mean,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$meanconf <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)

# create mean US-response 
ca <- aggregate(c$dcm1,by=list(c$id),FUN=mean,na.rm=T)
ca$id <- ca$Group.1
ca$id <- factor(ca$id)
ca$meanus <- ca$x
ca <- ca[c(3,4)]
c$id <- factor(c$id)
c <- merge(c,ca,all=T)

# create no correct responses
c$iscorrect <- ifelse(c$visresp==2, 1, 0)
ca <- ddply(c, "id", summarize,
            correctresp=sum(iscorrect, na.rm=TRUE)/sum(!is.na(visresp)),
            correctrespn=sum(iscorrect, na.rm=TRUE))
c <- merge(c, ca, all=TRUE)



# correct manipulations
c$manipulation <- factor(c$manipulation, levels=c("CFS", "No", "noCFS"))
c$manipulation[c$group=="Aware"] <- "noCFS"
c$manipulation <- factor(c$manipulation)

# correct stages
c$stage[c$stage=="AcquisitionOpp"] <- "Acquisition"
c$stage[c$stage=="ReversalOpp"] <- "Reversal"
c$stage[c$procedure=="ThankYou"] <- NA
c$stage[c$trial<0] <- NA
c$trial[!c$stage %in% c("Acquisition", "Reversal")] <- NA

# correct remaining stages
#c$stage[is.na(c$stage)&c$trial<17&c$trial>0&!is.na(c$spider)]  <- "Acquisition"
#c$stage[is.na(c$stage)&c$trial>16&!is.na(c$spider)]  <- "Reversal"

# correct stim
#c$stim[c$procedure=="CSminus"&c$stage=="Acquisition"] <- "CS-"
#c$stim[c$procedure=="CSplusUS"&c$stage=="Acquisition"] <- "CS+"
#c$stim[c$procedure=="CSminus"&c$stage=="Reversal"] <- "CS+"
#c$stim[c$procedure=="CSplusUS"&c$stage=="Reversal"] <- "CS-"
#c$stim[c$procedure=="CSminusRev"&c$stage=="Reversal"] <- "CS+"
#c$stim[c$procedure=="CSplusUSRev"&c$stage=="Reversal"] <- "CS-"
#c$stim[c$procedure=="CSminusRev"&c$stage=="Acquisition"] <- "CS-"
#c$stim[c$procedure=="CSplusUSRev"&c$stage=="Acquisition"] <- "CS+"

c$stim[c$procedure=="CSminus"&c$stage=="Acquisition"] <- "CSminus"
c$stim[c$procedure=="CSplusUS"&c$stage=="Acquisition"] <- "CSplus"
c$stim[c$procedure=="CSminus"&c$stage=="Reversal"] <- "CSplus"
c$stim[c$procedure=="CSplusUS"&c$stage=="Reversal"] <- "CSminus"
c$stim[c$procedure=="CSminusRev"&c$stage=="Reversal"] <- "CSplus"
c$stim[c$procedure=="CSplusUSRev"&c$stage=="Reversal"] <- "CSminus"
c$stim[c$procedure=="CSminusRev"&c$stage=="Acquisition"] <- "CSminus"
c$stim[c$procedure=="CSplusUSRev"&c$stage=="Acquisition"] <- "CSplus"

# The spider-terminology is really confusing
# so get this straight in the dataset
# "Spider A" is always CS+, "Spider B" is CS-
# irrespective of stage (acquisition or reversal)
# So, set spider-type that was used as initial CS+
# vertical spider type: X 
# horizontal spider type: >< 
c$firstunsafespidertype[c$order=="A"&c$preforder=="A"] <- "X"
c$firstunsafespidertype[c$order=="A"&c$preforder=="B"] <- "X"
c$firstunsafespidertype[c$order=="B"&c$preforder=="A"] <- "X"
c$firstunsafespidertype[c$order=="B"&c$preforder=="B"] <- "X"
c$firstunsafespidertype[c$order=="C"&c$preforder=="A"] <- "><"
c$firstunsafespidertype[c$order=="C"&c$preforder=="B"] <- "><"
c$firstunsafespidertype[c$order=="D"&c$preforder=="A"] <- "><"
c$firstunsafespidertype[c$order=="D"&c$preforder=="B"] <- "><"

# set ctrial
# sort c first
c <- c[order(c$id,c$trial),]
dt <- data.table(c)
#dt[,ctrial:=1:.N,by=c("id","stim")]
dt[,ctrial:=1:.N,by=c("id","spider")]
c <- as.data.frame(dt)
c <- c[order(c$id,c$trial),]

# this is kind of ugly but necessary
#c$ctrial[!(c$spider=="Spider A"|c$spider=="Spider B")]  <- NA
#c$ctrial[!c$spider=="Spider A"|spider]  <- NA
c$ctrial[is.na(c$spider)]  <- NA
c$ctrial[c$trial<1]  <- NA

c$id <- factor(c$id)
c$spider <- factor(c$spider)
c$group <- factor(c$group)
c$stage <- factor(c$stage)
c$visresp <- factor(c$visresp)
c$confresp <- factor(c$confresp)
#c$bleedlevel <- factor(c$bleedlevel)

# remove practice trials
c <- subset(c, trial>0)

# calculate # of trials per id
dfm <- ddply(c, "id", summarize, ntrials=sum(!is.na(dcm1)))

# merge
c <- merge(c, dfm, all=TRUE)

# sort file
c <- c[order(c$id,c$trial),] 

# write final data file
write.csv(c,"../data/cfs.csv",na="",row.names=F)
