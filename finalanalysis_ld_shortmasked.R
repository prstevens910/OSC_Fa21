# Load packages
library(tidyverse)

# Read in the data 
setwd("~/Documents/MorphReadingDev/202012_OSCemergence/lexDecExpt_fixedfillerissues/shortmasked/ld_data_shortmasked")
temp = list.files(pattern="*_complete_data.csv")
myfiles = lapply(temp, read.csv)


# Merge the csv files for all subjects
raw <- bind_rows(myfiles)

sonaIDs <- raw$responses[which(sapply(raw$responses,grepl,pattern="SonaID"))]
if (sum(duplicated(sonaIDs)) > 0) {
  print(paste("duplicate: ",sonaIDs[duplicated(sonaIDs)]))
  }

# read in stimlist
stimlist <- read.table("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen/finalstim_andchars_101_maxentropy.txt")

# read in the other OSC values (tcent) so I can compare 
osc.chars <- read.table("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen/primingandELPtargs_OSCcalc_compare_updated_CELEXandorth.txt",header=TRUE)
# stimlist$osc <- NULL
colnames(stimlist)[which(colnames(stimlist)=="osc")] <- "old.OSCccent" # this is because the new ccent accounts for cumulative freqs from CELEX. 
stimlist <- left_join(stimlist,select(osc.chars,c(targ,OSC.ccent,OSC.tcent)),by=c("target"="targ"))

### Read in Marelli's OSC values from 2018 paper
oscM18 <- read.table("~/Documents/MorphReadingDev/202003_LexiconProject/Marelli18_OSCcalculation/English_Marelli18_OSC.txt",header=TRUE)
oscM18$log_frequency <- NULL
colnames(oscM18) <- c("target","osc")
stimlist <- left_join(stimlist,oscM18,by="target")
# OSC.ccent and old.OSCccent are pretty similar so I don't need to worry too much about this difference 


# setwd so plots get saved to analysis folder
setwd("~/Documents/MorphReadingDev/202012_OSCemergence/analysis")

# go through and insert SonaIDs in their proper places 
start.inds <- which(sapply(raw$stimulus,grepl,pattern="Before we begin, please take a moment"))
if (length(sonaIDs) != length(start.inds)){print("alert! mismatch of sonaIDs and start inds")}
start.inds <- c(start.inds,nrow(raw)+1) # add final element for indexing
raw$participant <- ""
raw$dup <- ""
raw$list <- "list1"
for (i in 1:length(sonaIDs)) {
  raw$participant[start.inds[i]:(start.inds[i+1]-1)] <- sonaIDs[i]
  raw$dup[start.inds[i]:(start.inds[i+1]-1)] <- ifelse(i %in% which(duplicated(sonaIDs)),"_dup","")
  if (raw$target[which(raw$prime == "ALLURE" & raw$participant == sonaIDs[i])] == "hunt") {raw$list[start.inds[i]:(start.inds[i+1]-1)] <- "list2"}
}
getID <- function(x){strsplit(x,split='\"')[[1]][4]}
raw$participant <- sapply(raw$participant,getID) 
raw$participant <- paste(raw$participant,raw$dup,sep="")

# Quantify how many participants there are, and how many viewed
# list 2 versus list 1: 
print(paste("Participants:",length(unique(raw$participant))))
print(paste("List 1:",
            sum(raw$target[which(raw$prime == "ALLURE")] == "all")))
print(paste("List 2:", 
            sum(raw$target[which(raw$prime=="ALLURE")] == "hunt")))

# check the vision, colorblind, native questionnaire questions 
questionnaire <- raw %>% 
  filter(!is.na(response))
for (i in 1:nrow(questionnaire)) {
  if (grepl("vision",questionnaire$stimulus[i]) & questionnaire$response[i] == 1) {
    print(paste(questionnaire$participant[i],"said NO to vision Q"))
  }
  # if (grepl("colorblind",questionnaire$stimulus[i]) & questionnaire$button_pressed[i] == 0) {
  #   print(paste(questionnaire$participant[i],"said YES to colorblind Q"))
  # }
  if (grepl("native",questionnaire$stimulus[i]) & questionnaire$response[i] == 1) {
    print(paste(questionnaire$participant[i],"said NO to native English Q"))
  }
} # for 63 participants (final sample) just 149827 said no to vision Q 

# Get rid of all rows that don't have a target saved
trial.rows <- which(raw$target != "")
data <- raw[trial.rows,]

# Take out duplicate participants and pattest data
practice.data <- filter(data,sequence_type == "practice" & dup =="")

data <- filter(data,sequence_type!="practice" & dup == "") # there are no duplicates in this data (as of 7/31/22) so this just filters practice dat.

data$prime <- tolower(data$prime) # lowercase 

## some general checks

data %>% group_by(participant) %>% summarise(num = n()) %>% 
  filter(num != sum(data$participant==unique(data$participant)[1])) %>% nrow() 
# ^^ this should be 0 - all participants should have same number of data rows 


data$rt <- as.numeric(data$rt) # NAs are from missing responses (too slow!)

# What is the overall distribution of RTs? 
ggplot(data,aes(x=rt)) + geom_histogram() + 
  facet_wrap(vars(participant),nrow=4) + ggtitle("RT distributions by participant")
# weird dists: 148834, 148960**, 149821, 150001, 150064, 148834, 150016, 150331, 149290

## practice accuracy, trial accuracy, average RT, list assignment and key assignment for each participant. 
prac_acc <- practice.data %>%
  # filter(is_word=="yes") %>%
  group_by(participant) %>%
  summarise(mn = mean(correct=="true"), n_trials = n())
# identify LAST ITERATION for every participant - what was FINAL accuracy? 
final.accs <- c()
for (pu in unique(practice.data$participant)) {
  jo <- filter(practice.data,participant == pu)
  final.accs <- c(final.accs, mean(jo[c((nrow(jo)-16+1):nrow(jo)),]$correct == "true"))
}
prac_acc$finalacc <- final.accs
bad.prac <- filter(prac_acc,final.accs < 0.85)
sum( !(unique(bad.prac$participant) %in% p.rm))
min(bad.prac$finalacc[which(!(bad.prac$participant %in% p.rm))])
print("Practice Accuracy per Participant")
prac_acc
print("Participants with practice accuracy below 75% : ")
print(prac_acc[which(prac_acc$mn < 0.75 & prac_acc$n_trials < 18),]) # low acc and didn't repeat practice 
# 149080,149326,149560,149668,150046,150103,150160 # cool, only added one more person to this with final data collection. 

## TIMING ANALYSIS - How long did they take on average? 
bk1.times <- c()
bk2.times <- c()
bk3.times <- c()
total.times <- c()
for (pa in unique(data$participant)) { 
  total.times <- c(total.times,max(data$time_elapsed[data$participant == pa]))
  bk1.times <- c(bk1.times,max(data$time_elapsed[data$participant==pa & data$sequence_type == "bk1"]) - 
                   min(data$time_elapsed[data$participant==pa & data$sequence_type == "bk1"]))
  bk2.times <- c(bk2.times,max(data$time_elapsed[data$participant==pa & data$sequence_type == "bk2"]) - 
                   min(data$time_elapsed[data$participant==pa & data$sequence_type == "bk2"]))
  bk3.times <- c(bk3.times,max(data$time_elapsed[data$participant==pa & data$sequence_type == "bk3"]) - 
                   min(data$time_elapsed[data$participant==pa & data$sequence_type == "bk3"]))
}
# note - code below for calculating average total times
# filters out 3 participants who probably opened window in 
# browser and then didn't complete it until later that day (16+hours)
print(paste("average expt duration:",round(mean(total.times[which((total.times /6000)< 1000)] / 60000),2),"minutes"))
print(paste("average block 1 duration:",round(mean(bk1.times / 60000),2),"minutes"))
print(paste("average block 2 duration:",round(mean(bk2.times / 60000),2),"minutes"))
print(paste("average block 3 duration:",round(mean(bk3.times / 60000),2),"minutes"))

#### NOW look at their overall accuracy of responses - identify participants for removal 
acc <- data %>% 
  filter(is_word == "no") %>% # false positive 
  group_by(participant) %>%
  summarise(mn.acc = mean(correct=="true"))
print("Actual Accuracy per Participant")
acc
print("accuracies under 0.75")
paste(acc$participant[acc$mn.acc < 0.75],acc$mn.acc[acc$mn.acc < 0.75]) # false pos rate of greater than 25% 
p.rm <- acc$participant[acc$mn.acc < 0.75]
# 148960, 150355, 149791, 150016, 150394, 150064 # NOTE - before I had it as accuracy instead of false pos rate, and 150160 was included in bad parts. 

targ.acc <- data %>% 
  # filter(is_word=="yes") %>%
  group_by(target,is_word) %>%
  summarise(mn.acc = mean(correct=="true"))
print("Low-Accuracy Real Word Targets")
low.acc.wrds <- targ.acc$target[which(targ.acc$mn.acc < 0.5 & targ.acc$is_word=="yes")]
length(low.acc.wrds)
low.acc.wrds
print("Low-Accuracy Real Word Targets in Experimental Conditions")
low.acc.wrds[which(low.acc.wrds %in% data$target[which(data$prime_transparency %in% c("transparent","opaque"))])]
# need to remove 8 targets which had less that 50% accurate responses. one from transp (ZEAL) and 
# 7 from opaque (BRIG, DOLE, GALL, GLUT, POSIT, REND, WELT)
print("Low-Accuracy Nonword Targets")
low.acc.nonwrds <- targ.acc$target[which(targ.acc$mn.acc < 0.5 & targ.acc$is_word=="no")]
length(low.acc.nonwrds)
low.acc.nonwrds
targ.acc <- left_join(filter(targ.acc,is_word=="yes"),
                      select(stimlist,c("target","len","targfreq")),
                      by="target")
# remove words that are low-accuracy from analysis ?
data$lowaccwrd <- data$target %in% c(low.acc.wrds,low.acc.nonwrds)


avRT <- data %>%
  group_by(participant) %>%
  summarise(mn.rt=mean(rt,na.rm=TRUE)) # some NA RTs, probably from when participant took too long to respond 
print("Mean RT per Participant")
avRT$mn.rt
paste(avRT$participant[which(avRT$mn.rt > 900)],avRT$mn.rt[which(avRT$mn.rt > 900)]) # 150001 is aberrant, def needs to go
# participant 150001 has mean RT of >900, add to p.rm
p.rm <- c(p.rm, avRT$participant[which(avRT$mn.rt > 900)])

# Relationship between accuracy and speed? May help us identify ones who rushed, or weren't trying as hard
accRT <- left_join(acc,avRT,by="participant")
ggplot(accRT,aes(x=mn.rt,y=mn.acc,label=participant)) + 
  geom_point() + geom_text()
# verifies that the participants below 0.7 accuracy are baddies - outside the cluster, weirdly fast, weirdly slow or reeeeeeally inaccurate. 

# # Distribution of lists - list1 versus list2 
# lst1 <- sum(data$target=="brand" & data$related_prime == "yes")
# lst2 <- sum(data$target=="brand" & data$related_prime == "no")
# print("WORD LISTS")
# print(paste("list1: ",lst1,", list2: ",lst2,sep=""))

##  Last thing - distribution of assigned "yes" key. 
inst.inds <- which(sapply(raw$stimulus,grepl,pattern="key if the string is NOT a word"))
getkey <- function(s) {strsplit(strsplit(s,split="Press the ")[[1]][2], split=" key ")[[1]][1]}
kys <- data.frame(participant=raw$participant[inst.inds], yeskey = sapply(raw$stimulus[inst.inds],getkey))
print(paste("Proportion of participants assigned J as yes key:",round(mean(kys$yeskey == "J"),4)))
# 60% of participants were assigned J as "yes" key

# look at how often they pushed "yes" versus "no" 
if ("null" %in% unique(data$key_press)) {
  key.ops <- unique(data$key_press)[-which(unique(data$key_press)=="null")]
} else {key.ops <- unique(data$key_press)}
kprs <- data %>% 
  group_by(participant) %>%
  filter(key_press != "null") %>%
  summarise(mn=mean(key_press==key.ops[1]))
print("Key press per Subject")
as.numeric(kprs$mn) - 0.5 # (should all be close to 0.5 - half yes and half no)
if (sum(abs(kprs$mn - 0.5) > 0.15) > 0) {print("some participants are favoring one key")}
print(paste(kprs$participant[abs(kprs$mn-0.5)>0.15],kprs$mn[abs(kprs$mn-0.5)>0.15]))
# These participants are also ones with low accuracy and / or weird RT dists. 
print(paste(sum(!(kprs$participant[abs(kprs$mn-0.5)>0.15] %in% p.rm)), "of one-sided responders are not already in p.rm"))

# filter all of these because they're low accuracy / weird RT distributions / too slow responding.
p.rm <- c("222222","222223", p.rm) # add me to the list, since I'm a tester 
paste("participants in raw data:",length(unique(data$participant))) # 61 participants to start 
paste("participants removed:",length(p.rm)) # 7 participants removed due to low acc or too slow. Leaves *** 54 participants ***

# process final Qs : 
### Look at responses to final Qs 
# Question 1 : Press S if you saw anything else, K if you didn\'t.
sum(raw$key_press == "k" & !(raw$participant %in% p.rm))
sum(raw$key_press == "s" & !(raw$participant %in% p.rm))

nosee <- raw$participant[which(raw$key_press == "k" & !(raw$participant %in% p.rm))]
didsee <- raw$participant[which(raw$key_press == "s" & !(raw$participant %in% p.rm))]
# what did nosee participants say for questions 2, 3, 4? 
get.string.resp <- function(x){strsplit(x,":")[[1]][2]}
nosee.Q2.inds <- unname(which(sapply(raw$responses,grepl,pattern="descripQ2") & raw$participant %in% nosee))
# unname(sapply(raw$responses[nosee.Q2.inds],get.string.resp)) # all nas 
nosee.Q3.inds <- unname(which(sapply(raw$responses,grepl,pattern="descriptQ3") & raw$participant %in% nosee))
# unname(sapply(raw$responses[nosee.Q3.inds],get.string.resp)) # all nas 
nosee.Q4.inds <- unname(which(sapply(raw$responses,grepl,pattern="descriptQ4") & raw$participant %in% nosee))
# unname(sapply(raw$responses[nosee.Q4.inds],get.string.resp)) 
# ^ for Q4, a few realized they had seen an extra word on a few trials (5) , or "jumbled letters" 
# thought it was part of word generation process (1) 
didsee.Q2.inds <- unname(which(sapply(raw$responses,grepl,pattern="descripQ2") & raw$participant %in% didsee))
# unname(sapply(raw$responses[didsee.Q2.inds],get.string.resp)) # all nas 
didsee.Q3.inds <- unname(which(sapply(raw$responses,grepl,pattern="descriptQ3") & raw$participant %in% didsee))
# unname(sapply(raw$responses[didsee.Q3.inds],get.string.resp)) # all nas 
didsee.Q4.inds <- unname(which(sapply(raw$responses,grepl,pattern="descriptQ4") & raw$participant %in% didsee))
# unname(sapply(raw$responses[didsee.Q4.inds],get.string.resp)) 
# two said a few trials, two said "all", rest said "no" or "never" or "part of the word for a few trials but not often" 
# these are the participants who reported seeing "all" the prime words 
saw.prime <- raw$participant[didsee.Q4.inds[25]] # 25 reported consistently seeing word before target (so did 14 but changed answer from Q3 to Q4, so I don't they actually did )


###############################
#### BEGIN DISPLAY ANALYSIS ###
###############################

## NOW do some additional filtering for bad prime presentation estimates, and participants with too-slow frame rates. 
# I also need to remove participants with weirdly slow frame rates, and trials where display seems likely incorrect: 
data <- data %>% group_by(participant,sequence_type) %>% 
  mutate(trial_dur = time_elapsed - lag(time_elapsed, default=first(time_elapsed), 
                                        order_by = trial_index),
         starts_block = time_elapsed == min(time_elapsed)) %>% ungroup() # this estimates each trial duration - groups by participant and block so estimates aren't subtracting across breaks / participants
# turn rt into numeric so we can use RT and trial duration to estimate prime duration 

data$trial_dur[which(data$starts_block)] <- NA # if it's the first trial of the block, can't estimate trial duration because unclear how long prev trial was 
data$prime_dur <- data$trial_dur - data$rt - 2000 # 2000 for iti, fix, hash

# ok now question - for what proportion of trials is the prime dur estimate within 25 ms of the assigned SOA? 

if (sum(is.na(data$prime_dur)) != sum(is.na(data$trial_dur)) + sum(is.na(data$rt)) - sum(is.na(data$trial_dur) & is.na(data$rt))) {
  print("uh oh : mysterious NAs for prime_dur calculation")
}

sum(data$trial_dur < 2000 & ! data$starts_block) # there is one trial where the trial_dur is less than 2000, other than the ones where starts_block is true 
sum(data$prime_dur < 0, na.rm=TRUE)
# Ok I need a way to label the ones we're gonna throw out: 
data$weirddisp <- FALSE
data$weirddisp[which((data$prime_dur > 80 |  # 80 for 50 (predet SOA) + 30
                        data$prime_dur < 20 ) & # 20 for 50 - 30
                       !(data$starts_block))] <- TRUE # the starts_block trials can't have prime duration estimates calculated
mean(data$weirddisp) # 3.2% of data is less than 20 ms or more than 80 ms for prime duration estimates 

# We also need to identify the participants with too-slow frame rates: 
# soa.filt = 43
data %>% group_by(participant) %>% 
  summarize(mn.frame_time = mean(avg_frame_time),sd.frame_time = sd(avg_frame_time)) %>% # filter(mn.frame_time > 24) %>% .$participant
  ggplot(., aes(x=participant,y=mn.frame_time,ymin=mn.frame_time-sd.frame_time,ymax=mn.frame_time+sd.frame_time)) + 
  geom_pointrange() + # ylim(0,40) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Participant Average Frame Times")

# # remove all participants with avg frame durations above 17 ms (removes 7 participants, 2 of whom already in p.rm)
# weirddisp.parts <- data %>% 
#   group_by(participant) %>% 
#   summarize(mn.frame_time = mean(avg_frame_time),sd.frame_time = sd(avg_frame_time)) %>% 
#   filter(mn.frame_time > 25) %>% .$participant # replaced 17 ms with 25 ms for cut-off for frame time

data %>% # filter(participant %in% weirddisp.parts) %>% 
  group_by(participant) %>%
  summarize(mn.prime.dur = mean(prime_dur,na.rm=TRUE),sd.prime.dur = sd(prime_dur,na.rm=TRUE)) %>%
  ggplot(.,aes(x=participant,y=mn.prime.dur,ymin=mn.prime.dur-sd.prime.dur,ymax=mn.prime.dur+sd.prime.dur)) + 
  geom_pointrange() + ylim(0,100) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Participant Average Prime Times")

weirddisp.parts <- data %>% # filter(participant %in% weirddisp.parts) %>% 
  filter(prime_dur > 0) %>% # only two of these
  group_by(participant) %>%
  # mutate(mn.prime.dur = mean(prime_dur,na.rm=TRUE),sd.prime.dur = sd(prime_dur,na.rm=TRUE)) %>% 
  # ungroup() %>%
  # mutate(part.out = mn.prime.dur > mn.prime.dur + 3.5*sd.prime.dur | mn.prime.dur < mn.prime.dur  - 3.5*sd.prime.dur) %>% 
  # filter(!part.out) %>% # get rid of outliers 
  # group_by(participant) %>% # recalculate means with outliers removed 
  summarize(mn.prime.dur = mean(prime_dur,na.rm=TRUE),sd.prime.dur = sd(prime_dur,na.rm=TRUE)) %>% 
  filter(mn.prime.dur > 75 | mn.prime.dur < 50) %>% .$participant

# Are the ones who reported seeing "all" the prime words also in weirddisp.parts? 
sum(saw.prime %in% weirddisp.parts) # no they're not :(


# # below is for visualizing individual participants' distributions of prime_durations
# p = "150046"
# # p = sample(unique(data$participant)[which(!(unique(data$participant) %in% weirddisp.parts))])
# ggplot(filter(data,participant==p & prime_dur > 0),aes(x=prime_dur)) + geom_histogram() + ggtitle(p)
# # * 149080 has large spread, centered at 90 
# # 148930 might not be that bad - choppy but centered at 57 
# # 149668 also choppy but centered at 62, almost completely below 70 ms 
# # * 149869 centered at 80, widespread 
# # 149905 has a few VERY  negative values, otherwise all prime_durs between 50 and 80 
# # * 150046 is centered at 80, widespread (up to 100)


# note: 149905 has way more variable prime durs than anyone else, even though their frame rate is normal 
# looking at it more closely, they've just got a large sd because of one prime dur value that is super negative,
# which will def be filtered out, so no worries 

###############################
#### END DISPLAY ANALYSIS #####
###############################

# filter out the non-experimental trials, incorrect trials, too-slow trials, targets with low response acc 
e.data.msk <- filter(data,is_word == "yes" & # only want real words
                   correct == "true" & # only want correct answers
                   !is.na(rt) & # & # only want non-missing RT values
                   !participant %in% p.rm & # & get rid of questionably-performing participants
                   !lowaccwrd & 
                   !(prime_transparency %in% c("","undefined")))
nrow(e.data.msk)
e.data.msk$logrt <- log(e.data.msk$rt + 1)
# Get rid of by-participant outliers (lose 120 trials - 15 to <200ms and 105 to > 3.5*sd)
e.data.msk %>%
  group_by(participant) %>%
  mutate(mn.rt=mean(logrt),sd.rt=sd(logrt)) %>%
  ungroup() %>%
  mutate(part.out = logrt > mn.rt + 3.5*sd.rt | logrt < mn.rt - 3.5*sd.rt) -> e.data.msk
ggplot(e.data.msk,aes(x=rt,fill=part.out)) + geom_dotplot(dotsize = 0.3)
print(paste("Are you sure you want to get rid of",sum(e.data.msk$part.out),"by-participant outliers?"))
e.data.msk <- filter(e.data.msk,part.out == FALSE & rt >= 200)
e.data.msk$part.out <- NULL

# get rid of other now-unnecessary variables
e.data.msk$lowaccwrd <- NULL
e.data.msk$correct <- NULL
e.data.msk$is_word <- NULL
e.data.msk$dup <- NULL

# plot the prime dur ests, figure out what proportion of trials still has to be removed after filters applied above 
theme_set(theme_light(base_size = 20))
e.data.msk %>% filter(!(starts_block | trial_dur < 2000)) %>% # select(c(rt,trial_index,prime_dur,correct)) %>% head()
  ggplot(.,aes(x=prime_dur)) + geom_histogram()  + geom_vline(xintercept = 50,col="red",linetype="dashed",size = 1) +
  ggtitle("Fidelity of prime presentation durations") + xlim(0,120) + 
  labs(x = "Prime presentation duration estimates (ms)")
# ggsave("mskprime_presfid.png",last_plot(),device="png",width=8,height=6)
# filter out participants wit too long avg_frame_rates and weird disp trials
nrow(e.data.msk)
e.data.msk <- filter(e.data.msk,!participant %in% weirddisp.parts) # lose 3 participants 
nrow(e.data.msk)
e.data.msk <- filter(e.data.msk, !(participant %in% saw.prime)) # removes guy who said he saw every prime 
nrow(e.data.msk)
e.data.msk <- filter(e.data.msk, !weirddisp) # lose 112 trials for weird disp ests. 
nrow(e.data.msk)

# # update on list assignment after removing participants
# lst1 <- sum(e.data$target=="brand" & e.data$related_prime == "yes")
# lst2 <- sum(e.data$target=="brand" & e.data$related_prime == "no")
# print("WORD LISTS")
# print(paste("list1: ",lst1,", list2: ",lst2,sep=""))

# relationship between trial_index and RT 
ggplot(e.data,aes(x=trial_index,y=rt)) + geom_point()

# relationship between transparency and RT 
ggplot(e.data,aes(x=rt,col=prime_transparency)) + geom_density(alpha=0.6) + theme_light()

# relationship between relatedness and RT 
ggplot(e.data,aes(x=rt,col=related_prime)) + geom_density(alpha=0.6) + theme_light()


###################################################################
# MAIN ANALYSIS ###################################################
###################################################################

# add OSC to e.data info
e.data.msk <- left_join(e.data.msk,select(stimlist,c(target,numemb,len,MLfamsize,old,Afffamsize,primefreq,targfreq,targbg,OSC.ccent,osc)),by="target")
e.data.msk$numembscale <- e.data.msk$numemb / 1000


# save the data - this is what I'm posting for the data set to kilthub
e.data.msk.towrite <- select(e.data.msk,c(rt,response,trial_index,time_elapsed,avg_frame_time,
                                          sequence_type,related_prime,prime_transparency,target,prime,
                                          participant,list,prime_dur,len,MLfamsize,old,primefreq,
                                          targfreq,targbg,OSC.ccent))
e.data.msk.towrite$rt <- round(e.data.msk.towrite$rt,4)
e.data.msk.towrite$avg_frame_time <- round(e.data.msk.towrite$avg_frame_time,4)
colnames(e.data.msk.towrite)[6] <- "block"
colnames(e.data.msk.towrite)[13] <- "prime_dur_est"
write.table(e.data.msk.towrite,"study2_adultOSCeffects_data.txt",quote=FALSE,row.names=FALSE)


e.data <- e.data.msk
# split up data, since these stimuli weren't designed to be analyzed together 
op.dat <- filter(e.data,prime_transparency == "opaque")
nrow(op.dat) # 4107 trials , 2081 related
tr.dat <- filter(e.data,prime_transparency == "transparent")
nrow(tr.dat) # 4637 trials , 2320 related

# load necessary packages 
library(lme4) # use help(package = "lme4") to see version
library(lmerTest)
library(boot)

# some helper functions for visualization and model fitting 
mn.fun <- function(data, idx) { mean( data[ idx ] ) }

get.ci <- function(data) {
  bootstrp <- boot(data, mn.fun, R = 10000)
  ci <- boot.ci(boot.out = bootstrp, type="perc")$percent[1,4:5]
  return(ci)
}

se <- function(x){sd(x)/sqrt(length(x))}

invfn <- function() {
  ## link
  linkfun <- function(y) -1000/y
  ## inverse link
  linkinv <- function(eta)  -1000/eta
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 1000/(eta^2) }
  valideta <- function(eta) TRUE
  link <- "-1000/y"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

# what would be the most proper thing to do? 
# I should maybe first account for all the other things that are relevant 

########################################################################################
# fit models - first two-way and then three-way, for both 
# first opaque
start.time <- Sys.time()
opmod.IG.id <- glmer(rt ~ OSC.ccent*related_prime*targfreq + len + MLfamsize + 
                        (related_prime|target) + (related_prime|participant) , # maximal random effects structure
                       data = op.dat, family = inverse.gaussian( link = "identity" ),
                       glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
end.time <- Sys.time() # took 2.6 minutes 
print(end.time - start.time) # converged after intxn removed from RE
# summary(opmod.IG.id) 

# this exact model shows overal neg intxn of OSC and rel, pos 3-way of OSC, rel, targfreq 

# AIC      BIC   logLik deviance df.resid 
# 49281.6  49388.7 -24623.8  49247.6     4012 
# Number of obs: 4029, groups:  target, 93; participant, 51
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)                          685.412     16.810  40.774  < 2e-16 ***
#   OSC.ccent                            155.823     27.907   5.584 2.35e-08 ***
#   related_primeyes                      36.416     14.897   2.445  0.01450 *  
#   targfreq                              23.410      7.653   3.059  0.00222 ** 
#   len                                  -10.509      4.330  -2.427  0.01522 *  
#   MLfamsize                            -20.001     10.971  -1.823  0.06831 .  
# OSC.ccent:related_primeyes           -86.567     28.920  -2.993  0.00276 ** 
#   OSC.ccent:targfreq                   -94.611     14.765  -6.408 1.48e-10 ***
#   related_primeyes:targfreq            -19.119      8.498  -2.250  0.02446 *  
#   OSC.ccent:related_primeyes:targfreq   44.469     17.894   2.485  0.01295 *

# # looking into what happens with lots of predictors - very hard to converge, and results are weird (pos priming effect)
# fit not as good (AIC,BIC diff = 130 bigger than opmod.IG.id) 
# start.time <- Sys.time()
# opmod.IG.id.mega <- glmer(rt ~ OSC.ccent*related_prime*targfreq + OSC.ccent*related_prime*primefreq + len + MLfamsize + old + # Afffamsize + numembscale + 
#                        # (1|target) + 
#                          (1|participant) , # maximal random effects structure
#                      data = op.dat, family = inverse.gaussian( link = "identity" ),
#                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
# end.time <- Sys.time() # took 2.6 minutes 
# print(end.time - start.time)
# summary(opmod.IG.id.mega)

# # Checking results for inverse link function 
# start.time <- Sys.time()
# opmod.IG.inv <- glmer(rt ~ OSC.ccent*related_prime*targfreq + len + MLfamsize + 
#                        (related_prime|target) + (related_prime|participant) , # maximal random effects structure
#                      data = op.dat, family = inverse.gaussian( link = invfn() ),
#                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
# end.time <- Sys.time() # took 2.6 minutes 
# print(end.time - start.time) # converged after intxn removed from RE
# # summary(opmod.IG.inv) 
# 
# AIC(opmod.IG.id,opmod.IG.inv)
# BIC(opmod.IG.id,opmod.IG.inv)

# then transparent 
start.time <- Sys.time()
trmod.IG.id <- glmer(rt ~ OSC.ccent*related_prime*targfreq + len + MLfamsize + # primefreq + # old + Afffamsize + 
                        (related_prime|target) + (related_prime|participant) , # maximal random effects structure
                      data = tr.dat, family = inverse.gaussian( link = "identity" ),
                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
end.time <- Sys.time()
print(end.time - start.time)
summary(trmod.IG.id) # interaction of target frequency is weaker but in same dir as in opaque mod...

# AIC      BIC   logLik deviance df.resid 
# 55878.7  55988.2 -27922.3  55844.7     4620 
# Number of obs: 4637, groups:  target, 99; participant, 52
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
#    (Intercept)                          718.979     16.827  42.728   <2e-16 ***
#   OSC.ccent                             35.589     17.222   2.066   0.0388 *  
#   related_primeyes                      -7.226     14.335  -0.504   0.6142    
#   targfreq                             -17.338      6.950  -2.495   0.0126 *  
#   len                                   -9.480      4.291  -2.209   0.0272 *  
#   MLfamsize                              4.947     12.377   0.400   0.6894    
#   OSC.ccent:related_primeyes           -17.605     17.913  -0.983   0.3257    
#   OSC.ccent:targfreq                   -15.195     10.432  -1.457   0.1452    
#   related_primeyes:targfreq             -6.712      6.690  -1.003   0.3157    
#   OSC.ccent:related_primeyes:targfreq   17.608     11.001   1.601   0.1095 

# # Checking results for inverse link function 
# start.time <- Sys.time()
# trmod.IG.inv <- glmer(rt ~ OSC.ccent*related_prime*targfreq + len + MLfamsize + # primefreq + # old + Afffamsize + 
#                        (related_prime|target) + (related_prime|participant) , # maximal random effects structure
#                      data = tr.dat, family = inverse.gaussian( link = invfn() ),
#                      glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
# end.time <- Sys.time()
# print(end.time - start.time)
# # summary(trmod.IG.inv)

AIC(trmod.IG.id,trmod.IG.inv)
BIC(trmod.IG.id,trmod.IG.inv)



# ################ FOR TWO-WAY INTERACTION MODELS ###################
# For opaque 
start.time <- Sys.time()
opmod.IG.id.tw <- glmer(rt ~ OSC.ccent*related_prime + targfreq + len + MLfamsize +
                        (related_prime|target) + 
                         (related_prime|participant) , # maximal random effects structure
                     data = op.dat, family = inverse.gaussian( link = "identity" ),
                     glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
end.time <- Sys.time() # took 2.6 minutes
print(end.time - start.time) # converged after intxn removed from RE
summary(opmod.IG.id.tw) 


# AIC      BIC   logLik deviance df.resid 
# 50253.1  50341.6 -25112.6  50225.1     4093 
# Number of obs: 4107, groups:  target, 93; participant, 52
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
#   (Intercept)                 784.741     24.880  31.541  < 2e-16 ***
#   OSC.ccent                   -36.096     32.010  -1.128 0.259465    
#   related_primeyes             -3.633     16.172  -0.225 0.822231    
#   targfreq                    -17.517      5.180  -3.381 0.000721 ***
#   len                         -13.943      5.937  -2.349 0.018842 *  
#   MLfamsize                   -17.388     12.558  -1.385 0.166158    
#   OSC.ccent:related_primeyes    5.143     28.380   0.181 0.856187  

# # Check results for inv link function 
# start.time <- Sys.time()
# opmod.IG.inv.tw <- glmer(rt ~ OSC.ccent*related_prime + targfreq + len + MLfamsize + 
#                         (related_prime|target) + (related_prime|participant) , # maximal random effects structure
#                       data = op.dat, family = inverse.gaussian( link = invfn() ),
#                       glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
# end.time <- Sys.time() # took 2.6 minutes 
# print(end.time - start.time) # converged after intxn removed from RE
# # summary(opmod.IG.inv.tw) 
# 
# AIC(opmod.IG.id.tw,opmod.IG.inv.tw)
# BIC(opmod.IG.id.tw,opmod.IG.inv.tw)
# 
# For transparent 
start.time <- Sys.time()
trmod.IG.id.tw <- glmer(rt ~ OSC.ccent*related_prime + targfreq + len + MLfamsize +
                          # (1|target) + 
                          (1|participant) , # maximal random effects structure - had to decrease from rel | target
                        data = tr.dat, family = inverse.gaussian( link = "identity" ),
                        glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
end.time <- Sys.time() # took 2.6 minutes
print(end.time - start.time)
summary(trmod.IG.id.tw) # marginal interaction of OSC.ccent and prime_related 
# AIC      BIC   logLik deviance df.resid 
# 55888.7  55966.0 -27932.3  55864.7     4625 
# Number of obs: 4637, groups:  target, 99; participant, 52
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
#   (Intercept)                 737.577     21.597  34.153  < 2e-16 ***
#   OSC.ccent                    -5.432     15.342  -0.354   0.7233    
#   related_primeyes            -22.950      9.351  -2.454   0.0141 *  
#   targfreq                    -24.265      4.936  -4.916   8.83e-07 ***
#   len                         -10.025      4.654  -2.154   0.0312 *  
#   MLfamsize                     6.249     12.482   0.501   0.6167    
#   OSC.ccent:related_primeyes   28.584     15.923   1.795   0.0726 . 


# # Checking results for inv link function 
# start.time <- Sys.time()
# trmod.IG.inv.tw <- glmer(rt ~ OSC.ccent*related_prime + targfreq + len + MLfamsize + 
#                            (1|target) + (related_prime|participant) , # maximal random effects structure
#                          data = tr.dat, family = inverse.gaussian( link = invfn() ),
#                          glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
# end.time <- Sys.time() # took 2.6 minutes 
# print(end.time - start.time) # converged after intxn removed from RE
# summary(trmod.IG.inv.tw) 
# 
# AIC(trmod.IG.id.tw,trmod.IG.inv.tw)
# BIC(trmod.IG.id.tw,trmod.IG.inv.tw)

########################################################################################
# DATA VISUALIZATION 
# plot both data sets to show OSC * related_prime OR OSC * related_prime * targetfreq relationship

get.terc <- function(x, tercs) { # tercile func for plotting by targfreq 
  if (x == min(tercs)) { 
    i = 2
  } else if (x == max(tercs)) {
    i = length(tercs)
  } else { 
    i = min((which(x <= tercs)))
  }
  return(tercs[i-1] + (tercs[i] - tercs[i-1])/2) # this returns a centered value for the quintile instead of upper bound 
}

tf.tc <-round(quantile(op.dat$targfreq,c(0,0.5,1)),3)
TF.labs <- c("Low target frequency", "High target frequency")
names(TF.labs) <- c( "1.076", "3.4805")
osc.tc <- quantile(op.dat$OSC.ccent,c(0,0.2,0.4,0.6,0.8,1))
theme_set(theme_light(base_size = 20))
op.dat %>% select(related_prime,rt,OSC.ccent,target,targfreq) %>% 
  group_by(related_prime,target,OSC.ccent,targfreq) %>% 
  summarise(mn.rt = mean(rt)) %>%
  pivot_wider(id_cols=c(target,OSC.ccent,targfreq),names_from = related_prime,values_from=mn.rt) %>%
  mutate(prime.mag = no - yes, 
         targ.tile = get.terc(targfreq,tf.tc),
         osc.chunk = get.terc(OSC.ccent,osc.tc)) %>% 
  group_by(targ.tile,osc.chunk) %>% 
  mutate(mn.prime.mag = mean(prime.mag), ci = list(get.ci(prime.mag)), cinames = list(c("ci1","ci2"))) %>%
  ungroup() %>%
  unnest(cols = c(targ.tile,osc.chunk, mn.prime.mag, ci, cinames,prime.mag,OSC.ccent)) %>% 
  pivot_wider(id_cols = c(targ.tile, osc.chunk, mn.prime.mag,prime.mag,OSC.ccent), 
              names_from = cinames, 
              values_from = ci) %>%
  ggplot(.,aes(x=OSC.ccent,y=prime.mag)) + 
  labs(x= "OSC",y="Priming magnitude (Unrelated - related, ms)",title="Effect of OSC on opaque masked priming") +
  # geom_point(size=0.3) + 
  geom_smooth(method="glm",se=FALSE,col="red") + # This is fitted to actual OSC.ccent and prime.mag vals, not quintile means 
  geom_point(aes(x=osc.chunk,y=mn.prime.mag)) + geom_line(aes(x=osc.chunk,y=mn.prime.mag)) + 
  geom_errorbar(aes(x=osc.chunk,y=mn.prime.mag,ymin=ci1,ymax=ci2),width=0.05) + 
  facet_grid(.~targ.tile , labeller = labeller(targ.tile = TF.labs)) +  xlim(0.2,0.85)
ggsave("OSCmskprimemags_opaque_bytarg.png",last_plot(),device="png",width=9,height=7)

## Post-hoc analysis: filtering by target frequency, do we see sig OSC*rel intxn for either lower-freq or higher-freq half? 
start.time <- Sys.time()
opmod.IG.id.tw.ltf <- glmer(rt ~ OSC.ccent*related_prime + len + MLfamsize +
                          # (1|target) + 
                            (1|participant) , # maximal random effects structure
                        data = filter(op.dat,targfreq <= tf.tc[2]), family = inverse.gaussian( link = "identity" ),
                        glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
end.time <- Sys.time() # took 2.6 minutes
print(end.time - start.time) # converged after intxn removed from RE
summary(opmod.IG.id.tw.ltf) # neither the lower half nor the upper half yielded sig results (both p's > 0.4)


# plot rts to show OSC * related_prime * targetfreq relationship
tf.tc <-round(quantile(tr.dat$targfreq,c(0,0.5,1)),2)
TF.labs <- c("Low target frequency", "High target frequency")
names(TF.labs) <- c( "1.37", "3.355")
osc.tc <- quantile(tr.dat$OSC.ccent,c(0,0.2,0.4,0.6,0.8,1))
theme_set(theme_light(base_size = 20))
tr.dat %>% select(related_prime,rt,OSC.ccent,target,targfreq) %>%
  group_by(related_prime,target,OSC.ccent,targfreq) %>%
  summarise(mn.rt = mean(rt)) %>%
  pivot_wider(id_cols=c(target,OSC.ccent,targfreq),names_from = related_prime,values_from=mn.rt) %>%
  mutate(prime.mag = no - yes,
         targ.tile = get.terc(targfreq,tf.tc),
         osc.chunk = get.terc(OSC.ccent,osc.tc)) %>%
  group_by(targ.tile,osc.chunk) %>%
  mutate(mn.prime.mag = mean(prime.mag), ci = list(get.ci(prime.mag)), cinames = list(c("ci1","ci2"))) %>%
  ungroup() %>%
  unnest(cols = c(targ.tile,osc.chunk, mn.prime.mag, ci, cinames,prime.mag,OSC.ccent)) %>%
  pivot_wider(id_cols = c(targ.tile, osc.chunk, mn.prime.mag,prime.mag,OSC.ccent),
              names_from = cinames,
              values_from = ci) %>%
  ggplot(.,aes(x=OSC.ccent,y=prime.mag)) +
  labs(x= "OSC",y="Priming magnitude (Unrelated - related, ms)",title="Effect of OSC on transparent masked priming") +
  # geom_point(size=0.3) +
  geom_smooth(method="glm",se=FALSE,col="red") + # This is fitted to actual OSC.ccent and prime.mag vals, not quintile means
  geom_point(aes(x=osc.chunk,y=mn.prime.mag)) + geom_line(aes(x=osc.chunk,y=mn.prime.mag)) +
  geom_errorbar(aes(x=osc.chunk,y=mn.prime.mag,ymin=ci1,ymax=ci2),width=0.05) +
  facet_grid(.~targ.tile, labeller = labeller(targ.tile = TF.labs)) + xlim(0.2,0.85)
ggsave("OSCmskprimemags_transparent_bytarg.png",last_plot(),device="png",width=9,height=7)

# plot transp prime mags to show OSC * related_prime relationship
osc.tc <- quantile(tr.dat$OSC.ccent,c(0,0.2,0.4,0.6,0.8,1))
theme_set(theme_light(base_size = 20))
tr.dat %>% select(related_prime,rt,OSC.ccent,target) %>%
  group_by(related_prime,target,OSC.ccent) %>%
  summarise(mn.rt = mean(rt)) %>%
  pivot_wider(id_cols=c(target,OSC.ccent),names_from = related_prime,values_from=mn.rt) %>%
  mutate(prime.mag = no - yes,
         osc.chunk = get.terc(OSC.ccent,osc.tc)) %>%
  group_by(osc.chunk) %>%
  mutate(mn.prime.mag = mean(prime.mag), ci = list(get.ci(prime.mag)), cinames = list(c("ci1","ci2"))) %>%
  ungroup() %>%
  unnest(cols = c(osc.chunk, mn.prime.mag, ci, cinames,prime.mag,OSC.ccent)) %>%
  pivot_wider(id_cols = c(osc.chunk, mn.prime.mag,prime.mag,OSC.ccent),
              names_from = cinames,
              values_from = ci) %>%
  ggplot(.,aes(x=OSC.ccent,y=prime.mag)) +
  labs(x= "OSC",y="Priming magnitude (Unrelated - related, ms)",title="Effect of OSC on transparent masked priming") +
  # geom_point(size=0.3) +
  geom_smooth(method="glm",se=FALSE,col="red") + # This is fitted to actual OSC.ccent and prime.mag vals, not quintile means
  geom_point(aes(x=osc.chunk,y=mn.prime.mag)) + geom_line(aes(x=osc.chunk,y=mn.prime.mag)) +
  geom_errorbar(aes(x=osc.chunk,y=mn.prime.mag,ymin=ci1,ymax=ci2),width=0.05)  + scale_y_continuous(limits=c(-20,45)) + scale_x_continuous(limits=c(0.2,0.85))# +
  # facet_grid(.~targ.tile, labeller = labeller(targ.tile = TF.labs)) + xlim(0.2,0.85)
ggsave("OSCmskprimemags_transparent.png",last_plot(),device="png",width=8,height=6)

# plot transp prime mags to show OSC * related_prime relationship
osc.tc <- quantile(op.dat$OSC.ccent,c(0,0.2,0.4,0.6,0.8,1))
theme_set(theme_light(base_size = 20))
op.dat %>% select(related_prime,rt,OSC.ccent,target) %>%
  group_by(related_prime,target,OSC.ccent) %>%
  summarise(mn.rt = mean(rt)) %>%
  pivot_wider(id_cols=c(target,OSC.ccent),names_from = related_prime,values_from=mn.rt) %>%
  mutate(prime.mag = no - yes,
         osc.chunk = get.terc(OSC.ccent,osc.tc)) %>%
  group_by(osc.chunk) %>%
  mutate(mn.prime.mag = mean(prime.mag), ci = list(get.ci(prime.mag)), cinames = list(c("ci1","ci2"))) %>%
  ungroup() %>%
  unnest(cols = c(osc.chunk, mn.prime.mag, ci, cinames,prime.mag,OSC.ccent)) %>%
  pivot_wider(id_cols = c(osc.chunk, mn.prime.mag,prime.mag,OSC.ccent),
              names_from = cinames,
              values_from = ci) %>%
  ggplot(.,aes(x=OSC.ccent,y=prime.mag)) +
  labs(x= "OSC",y="Priming magnitude (Unrelated - related, ms)",title="Effect of OSC on opaque masked priming") +
  # geom_point(size=0.3) +
  geom_smooth(method="glm",se=FALSE,col="red") + # This is fitted to actual OSC.ccent and prime.mag vals, not quintile means
  geom_point(aes(x=osc.chunk,y=mn.prime.mag)) + geom_line(aes(x=osc.chunk,y=mn.prime.mag)) +
  geom_errorbar(aes(x=osc.chunk,y=mn.prime.mag,ymin=ci1,ymax=ci2),width=0.05) + scale_y_continuous(limits=c(-20,45)) + scale_x_continuous(limits=c(0.2,0.85))
# facet_grid(.~targ.tile, labeller = labeller(targ.tile = TF.labs)) + xlim(0.2,0.85)
ggsave("OSCmskprimemags_opaque.png",last_plot(),device="png",width=8,height=6)
# 
# osc.tc <- quantile(tr.dat$OSC.ccent,c(0,0.2,0.4,0.6,0.8,1))
# theme_set(theme_light(base_size = 20))
# tf.tc <-round(quantile(tr.dat$targfreq,c(0,0.5,1)),2)
# TF.labs <- c("Low target frequency", "High target frequency")
# names(TF.labs) <- c( "1.37", "3.355")
# tr.dat %>% select(related_prime,rt,OSC.ccent,target,targfreq) %>%
#   group_by(related_prime,target,OSC.ccent,targfreq) %>%
#   summarise(mn.rt = mean(rt)) %>%
#   pivot_wider(id_cols=c(target,OSC.ccent,targfreq),names_from = related_prime,values_from=mn.rt) %>%
#   mutate(prime.mag = no - yes,
#          targ.tile = get.terc(targfreq,tf.tc),
#          osc.chunk = get.terc(OSC.ccent,osc.tc)) %>%
#   group_by(targ.tile,osc.chunk) %>%
#   mutate(mn.prime.mag = mean(prime.mag), se.prime.mag = se(prime.mag)) %>%
#   ungroup() %>%
#   ggplot(.,aes(x=OSC.ccent,y=prime.mag)) +
#   labs(x= "OSC",y="Priming magnitude (Unrelated - related, ms)",title="Effect of OSC on transparent masked priming") +
#   # geom_point(size=0.3) +
#   geom_smooth(method="glm",se=FALSE,col="red") + # This is fitted to actual OSC.ccent and prime.mag vals, not quintile means
#   geom_point(aes(x=osc.chunk,y=mn.prime.mag)) + geom_line(aes(x=osc.chunk,y=mn.prime.mag)) +
#   geom_errorbar(aes(x=osc.chunk,y=mn.prime.mag,ymin=(mn.prime.mag-se.prime.mag),ymax=(mn.prime.mag+se.prime.mag)),width=0.05) +
#   facet_grid(.~targ.tile , labeller = labeller(targ.tile = TF.labs)) + xlim(0.2,0.85)
# ggsave("OSCprimemags_transparent_bytarg_se.png",last_plot(),device="png",width=9,height=7)

