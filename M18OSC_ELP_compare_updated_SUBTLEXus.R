library(tidyverse)

### Read in Marelli's OSC values from 2018 paper
oscM18 <- read.table("~/Documents/MorphReadingDev/202003_LexiconProject/Marelli18_OSCcalculation/English_Marelli18_OSC.txt",header=TRUE)
oscM18$log_frequency <- NULL
colnames(oscM18) <- c("Word","OSC")

### Read in priming data from studies. 
# Read in Rastle 2004 data
rast04 <- read.csv("/Users/patience/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/MorphNet/testing_stim/Rastle-stim.csv")
rast.pem <- filter(rast04,prime.type!="control")
rast.pem$RTms <- rast04$RTms[seq(2,nrow(rast04),2)] - rast04$RTms[seq(1,nrow(rast04),2)] # fixed this - was rel - unrel before 
colnames(rast.pem)[3] <- "PEM"
take.prime <- function(x){strsplit(x,"-")[[1]][1]}
take.targ <- function(x){strsplit(x,"-")[[1]][2]}
rast.pem$prime <- sapply(rast.pem$name,take.prime)
# Take unrelated means as well so I can look at unrelated prime RTs in addition to ELP word-alone RTs. 
rast.pem$unrel.RT <- rast04$RTms[seq(2,nrow(rast04),2)]

# Read in Rastle 2000 data
rast00 <- read.csv("/Users/patience/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/MorphNet/testing_stim/Rastle2000-stim.csv")
rast00$pem43 <- rast00$contrl - rast00$SOA1.43ms
rast00$pem72 <- rast00$contrl.1 - rast00$SOA2.72ms
rast00$pem230 <- rast00$contrl.2 - rast00$SOA3.230ms

pem.toadd <- rast00 %>% 
  select(c("condition","prime","target","pem43","contrl")) %>% 
  mutate(name=paste(prime,target,sep="-"))
colnames(pem.toadd) <- c("prime.type","prime","targ","PEM","unrel.RT","name")

# Read in Jared 2017 data! 
jar17 <- read.csv("/Users/patience/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/MorphNet/testing_stim/Jared-Jouravlev-Joanisse-AllLDT.csv")
# fill in prime words from the stimulus set I already have transcribed - match by prime word. 
jar17$Target <- tolower(jar17$Target)
jar17$TargetType <- tolower(jar17$TargetType)
jar17 <- jar17[-duplicated(select(jar17,c(Experiment,Subject,List,Condition,Target,TargetType))),] # remove one repeated row. 
jar.mns <- jar17 %>% 
  group_by(Condition,Target,TargetType) %>%
  summarise(mn.rt = mean(RTLME, na.rm=TRUE)) %>% 
  spread(Condition,mn.rt) %>% 
  mutate(pem = unrelated - related)

# read in Jared stimuli so I can insert the prime word 
jarstim <- read.csv("/Users/patience/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/MorphNet/testing_stim/jared-stim.csv")
jarstim$Prime <- sapply(jarstim$name,take.prime)
jarstim$Target <- sapply(jarstim$name,take.targ)
jar.mns <- left_join(jar.mns,select(jarstim,c(name,Prime,Target)),by="Target")
jar.mns$Prime[which(is.na(jar.mns$Prime))] <- c("bargain","capsize","cashew","colourful","diverse",
                                                "early","favourable","lapse","lately","pastel","",
                                                "season","venture")
jar.mns <- jar.mns[-which(jar.mns$Prime == ""),]
colnames(jar.mns) <- c("targ","prime.type","related","unrel.RT","PEM","name","prime")
jar.mns$name <- paste(jar.mns$prime,jar.mns$targ,sep="-")
jar.mns$prime.type[which(jar.mns$prime.type == "control")] <- "orthographic"

# read in Andrews & Lo stimuli 
and.mns <- read.csv("/Users/patience/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/MorphNet/testing_stim/AL2013_itemaverages_summary_data.csv")
and.mns$X..Error <- NULL
and.mns.rel <- filter(and.mns,PrimeType == "Rel")
and.mns.rel$unrelated <- filter(and.mns,PrimeType == "Unrel")$Average.RT
and.mns.rel$PrimeType <- NULL
and.mns.rel$Excluded <- NULL
colnames(and.mns.rel) <- c("prime.type","targ","prime","related","unrel.RT")
and.mns.rel$PEM <- and.mns.rel$unrel.RT - and.mns.rel$related
and.mns.rel$targ <- sapply(and.mns.rel$targ,tolower)
and.mns.rel$name <- paste(and.mns.rel$prime,and.mns.rel$targ,sep="-")
and.mns <- and.mns.rel
and.mns$prime.type[which(and.mns$prime.type=="Form")] <- "orthographic"
and.mns$prime.type[which(and.mns$prime.type=="Opaque")] <- "opaque"
and.mns$prime.type[which(and.mns$prime.type=="Transp")] <- "transparent"

# merge Rastle 2004, Rastle 2000 43 ms data, and Jared 2017 into one big beautiful df. 
# combine 2000 and 2004 data, get rid of all but morphological conditions 
rast.pem$targ <- sapply(rast.pem$name,take.targ)
rast.pem$study <- "rast04"
pem.toadd$study <- "rast00"
jar.mns$study <- "jar17"
and.mns$study <- "and13"
rast.pem <- rbind(rast.pem,pem.toadd, select(and.mns,c(targ,prime.type,PEM,unrel.RT,name,prime,study)),
                  select(jar.mns,c(targ,prime.type,PEM,unrel.RT,name,prime,study))) %>% filter(prime.type %in% c("orthographic","opaque","quasi","transparent"))


### Read in SUBTLEX vocabulary
voc_SUBTLEX <- read.csv("~/Documents/MorphReadingDev/202003_LexiconProject/SUBTLEXus/SUBTLEX-US.csv") # SUBTLEX- retrieved from  
vocS <- select(voc_SUBTLEX,c(Word,SUBTLWF)) # take just the columns I need: Word, COB Frequency (raw, not per million words, for more precision)

# Also read in CELEX vocabulary 
voc_CELEX <- read.csv("~/Documents/MorphReadingDev/202003_LexiconProject/celex2/english/efw/efw.cd",
                      header=FALSE,sep="\\") # COBUILD wordform freqs from CELEX, based on 17.9 million words. 
vocC <- select(voc_CELEX,c(V2,V4)) # take just the columns I need: Word, COB Frequency (raw, not per million words, for more precision)
colnames(vocC) <- c("Word","COBfreq") 
# change voc so that there's just one of each word, and the frequencies are sum of freqs for each word
sumfreq <- function(x){return(sum(vocC$COBfreq[vocC$Word==x]))}
upcfs <- sapply(unique(vocC$Word),sumfreq) # unique prime COB freq sums 
vocC <- data.frame(Word=unique(vocC$Word),COBfreq = upcfs)

# only get words that are in stim$target and stim$prime - those are the only ones I care about. 
# get rid of words with non-letter characters
spec <- unlist(sapply(vocS$Word,grepl,pattern="[^a-z]")) # instead of "[^a-z]"
if (sum(spec)>0){vocS <- vocS[-which(spec),]} # leaves 74285 words 

spec <- unlist(sapply(vocC$Word,grepl,pattern="[^a-z]"))
if (sum(spec)>0){vocC<-vocC[-which(spec),]} # leaves 64735 words 

# the vocC and vocS vocabularies have 46986 words in common...

# ok - lets just use CELEX , and focus on the words that CELEX and SUBTLEX have in common

voc <- vocC %>% filter(Word %in% vocS$Word)# make CELEX the main frequency, filter out words not in SUBTLEX

### Identify all the words that have words from prime-target pairs embedded
# For each target, figure out which words in voc have it as a substring
# NOTE - I'm using the smaller number of words to identify ortho neighbors and the bigger CELEX voc to calculate freqs 
print("identify orthographic neighbors...")
targ.list <- unique(c(rast.pem$targ,oscM18$Word))
# targ.list <- c("whisk")
cont <- rep(0,length(targ.list)*nrow(voc))
for (t in 1:length(targ.list)) {
  if (t %% 100 == 0) {print(paste(t, "out of",length(targ.list)))}
  cont[((t-1)*nrow(voc)+1):(t*nrow(voc))] <- as.numeric(sapply(voc$Word,grepl,pattern=targ.list[t]))
}

# # now, let's turn this into a matrix. 

# setwd("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen")
# write.table(cont,"OSCcont_allELPcalc.txt") # uncomment this if there's a need for it to be written but note it takes a long time. 

cont.mat <- matrix(cont,nrow=length(targ.list),byrow=TRUE)
dimnames(cont.mat) <- list(targ.list,voc$Word)
ngdf <- voc[which(colSums(cont.mat)>0),]
ngdf <- rbind(ngdf,vocC[which(vocC$Word %in% targ.list[!(targ.list %in% ngdf$Word)]),]) # adding the target words that aren't in voc but are in vocC

# ##########################################################################
# ####******* The code below was used to generate "whisk" example of ***####
# ####******* difference in reference point for the two methods.     ***####
# # Note: you'll need glovevecs read in for this, but that's it. 

# whisk.words <- c("whisk","whisked","whisker","whiskered","whiskers","whiskies","whisking","whisks","whisky")
# whisk.centroid <- colSums(glovevecs$logCOBfreq[which(glovevecs$Word %in% whisk.words)]*
#                             glovevecs[which(glovevecs$Word %in% whisk.words),2:301])/
#   sum(glovevecs$logCOBfreq[which(glovevecs$Word %in% whisk.words)])
# ind.df <- data.frame(cmdscale(dist(
#   rbind(
#                 glovevecs[which(glovevecs$Word %in% whisk.words),2:301],
#                 whisk.centroid
#                 ))))
# ind.df$Word <- c(glovevecs$Word[which(glovevecs$Word %in% whisk.words)],"centroid")
# ind.df <- left_join(ind.df,select(glovevecs,c(Word,logCOBfreq)),by="Word")
# colnames(ind.df) <- c("Dim1","Dim2","Word","Frequency")
# # wiggle some of the words a bit so they don't overlap. 
# ind.df$Dim2[ind.df$Word=="whisks"] <- ind.df$Dim2[ind.df$Word=="whisks"] -0.2
# ind.df$Dim2[ind.df$Word=="whisky"] <- ind.df$Dim2[ind.df$Word=="whisky"] - 0.4
# ind.df$Dim2[ind.df$Word=="whisking"] <- ind.df$Dim2[ind.df$Word=="whisking"] - 0.2
# ind.df$Dim1[ind.df$Word=="whisky"] <- ind.df$Dim1[ind.df$Word=="whisky"] - 0.15
# ggplot(filter(ind.df,Word!="centroid"),aes(x=Dim1,y=Dim2,size=Frequency)) +
#   theme_light() + 
#   theme(legend.position = "none",axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
#   geom_point() +
#   xlim(-7.1,4) +
#   geom_text(size=8,nudge_y=0.4,nudge_x=-0.8,aes(label=Word)) +
#   geom_point(shape="x", color="red",size=11,show.legend=FALSE,
#                             aes(x=ind.df$Dim1[which(ind.df$Word=="centroid")],
#                                 y=ind.df$Dim2[which(ind.df$Word=="centroid")])) +
#   geom_point(shape="+", color="red",size=11,show.legend=FALSE,
#              aes(x=ind.df$Dim1[which(ind.df$Word=="whisk")],
#                  y=ind.df$Dim2[which(ind.df$Word=="whisk")]))
# 
# ggsave("whiskplot_OSCreferencepointscomparison.png",last_plot(),width=5.0,height=5.0)

#################******** END of whiskplot Code ********##################
##########################################################################


### Extract Glove vectors for all the words which contain one of these words
print("Extract glovevecs...")
glove.dims<-300
glovevecs <- data.frame(Word=as.character(ngdf$Word),matrix(rep(0.0,glove.dims*nrow(ngdf)),nrow=nrow(ngdf),ncol=glove.dims))
inputFile <- "~/Documents/MorphReadingDev/201902_ModelingDevMorphEffects_authenticvocabulary/glove.42B.300d.txt"
con  <- file(inputFile, open = "r")

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  myVector <- strsplit(oneLine, " ")[[1]]
  if (myVector[1] %in% glovevecs$Word) {
    glovevecs[which(glovevecs$Word == as.character(myVector[1])),2:(glove.dims+1)] <- as.numeric(myVector[2:(glove.dims+1)])
  }
}
close(con)

# get rid of the ones that didn't appear in glove reps
empty <- which(rowSums(glovevecs[,2:301]) == 0.0)
glovevecs_all <- glovevecs
glovevecs <- glovevecs[-empty,]
glovevecs <- left_join(glovevecs,vocC,by="Word")
glovevecs$COBfreq <- glovevecs$COBfreq + 1 # add one so you're never dividing by 0 while calculating tcent
glovevecs$logCOBfreq <- log10(glovevecs$COBfreq) + 0.0001 # add 0.0001 for same problem while calculating ccent.
setwd("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen")
# write.table(glovevecs,paste("glovevecs_forprimingandELPtargs_",nrow(glovevecs),"words_plusCOBfreq.txt",sep=""))
glovevecs <- read.table("glovevecs_forprimingandELPtargs_47383words_plusCOBfreq.txt")

# define functions for osc calculations
# # cosine similarity
cos.sim <- function(A,B) {return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )}

# # wtd cosine similarity
wtd.cos <- function(w,t){c(cos.sim(as.numeric(glovevecs[which(glovevecs$Word==t),2:301]),
                                   as.numeric(glovevecs[which(glovevecs$Word==w),2:301])),
                           glovevecs$COBfreq[which(glovevecs$Word==w)])}

# for every word in OSCM18, use cont.mat and glovevecs to calculate that word's OSC - target-centered and cluster-centered methods. 
tocalc <- which(targ.list %in% glovevecs$Word) # if the embedded word isn't in glovevecs, no point in trying to calculate osc. 
OSC.targs <- data.frame(targ = targ.list[tocalc],OSC.tcent = rep(NA,length(tocalc)),
                        OSC.ccent = rep(NA,length(tocalc)),numembed = rep(NA,length(tocalc)))
# colnames(oscM18)[1] <- "targ"
# rast.pem <- left_join(rast.pem,oscM18,by="targ")
# rast.pem$OSC.tcent <- NA
# oscM18$OSC.tcent <- NA # target-centered OSC calculation
# rast.pem$OSC.ccent <- NA
# oscM18$OSC.ccent <- NA # cluster-centered OSC calculation
# oscM18$numembed <- NA # number of words with the target embedded 
# rast.pem$numembed <- NA


# Go through and calculate all of these words' OSC using both target-centered and cluster-centered methods. (note: log freq for targ-cent?)
print("calculate OSC...")
for (s in tocalc) { # for every stimulus target I've got data for 
  si = which(tocalc == s) # figure out where in the list we are. 
  # progress report
  if (si %% 50 == 0){print(paste(si,"out of",length(tocalc)))}
  # get all the words that contain the root 
  if (targ.list[s] %in% glovevecs$Word) { # leave them both as NA if the embedded word is NOT in glovevecs. 
    wrds <- as.character(colnames(cont.mat)[which(cont.mat[targ.list[s],]==1)]) # find the CELEX words that contain the embedded string. 
    if (!(targ.list[s] %in% wrds)) { wrds <- c(wrds,targ.list[s]) } # if target isn't in there (bc not in constrained voc), add it in - it must be in glovevecs 
    targlove <- glovevecs[which(glovevecs$Word %in% wrds),2:301]
    if (nrow(targlove)==1) { # if there's only one word, we know it's the embedded word. both OSCs should be 1
      # oscM18$OSC.tcent[s] <- 1
      # rast.pem$OSC.tcent[s] <- 1
      OSC.targs$OSC.tcent[si] <- 1
      # oscM18$OSC.ccent[s] <- 1
      # rast.pem$OSC.ccent[s] <- 1
      OSC.targs$OSC.ccent[si] <- 1
    } else {# otherwise do distinct calculations for each. 
      # Cluster-centered calculation 
      wtd.mn <- colSums(glovevecs$logCOBfreq[which(glovevecs$Word %in% wrds)]*targlove)/
        sum(glovevecs$logCOBfreq[which(glovevecs$Word %in% wrds)])
      # oscM18$OSC.ccent[s] <- mean(apply(targlove,1,cos.sim,B=wtd.mn))
      # rast.pem$OSC.ccent[s] <- mean(apply(targlove,1,cos.sim,B=wtd.mn))
      OSC.targs$OSC.ccent[si] <- mean(apply(targlove,1,cos.sim,B=wtd.mn))
      
      # Target-centered calculation 
      # cosmat <- sapply(wrds[which(wrds != oscM18$Word[s] & wrds %in% glovevecs$Word)],wtd.cos,t=oscM18$Word[s]) # target not included
      cosmat <- sapply(wrds[which(wrds %in% glovevecs$Word)],wtd.cos,t=targ.list[s]) # target included
      # oscM18$OSC.tcent[s] <- sum(cosmat[2,]*cosmat[1,])/sum(cosmat[2,])
      # rast.pem$OSC.tcent[s] <- sum(cosmat[2,]*cosmat[1,])/sum(cosmat[2,])
      OSC.targs$OSC.tcent[si]<- sum(cosmat[2,]*cosmat[1,])/sum(cosmat[2,])
    }
    OSC.targs$numembed[si] <- nrow(targlove)
  }
  
  # if (sum(is.na(oscM18$OSC.ccent[s]),is.na(oscM18$OSC.tcent[s]))) {# if one or the other is NA but not both. 
  #   # print out some info
  #   print(paste(oscM18$Word[s],", glovevecs:",nrow(targlove),", wrds:",length(wrds)))
  # }
  
  if (sum(is.na(OSC.targs$OSC.ccent[si]),is.na(OSC.targs$OSC.tcent[si]))) {# if one or the other is NA but not both. 
    # print out some info
    print(paste(OSC.targs$Word[si],", glovevecs:",nrow(targlove),", wrds:",length(wrds)))
  }
}
write.table(OSC.targs,"primingandELPtargs_OSCcalc_compare_updated_CELEXandorth.txt",quote = FALSE,row.names = FALSE)
# write.table(oscM18,"allELP_OSCcalc_compare_updated.txt",quote = FALSE,row.names = FALSE)

# read in OSC targs instead of generating them .
setwd("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen")
OSC.targs <- read.table("primingandELPtargs_OSCcalc_compare_updated_CELEXandorth.txt",header=TRUE)

# Then get the word characteristics from ELP website, put them here - this has all the words in OSC.targs that exist in ELP
setwd("~/Documents/MorphReadingDev/202003_LexiconProject/English/ELP_wordchars")
ELP.chars <- read.csv("Items.csv")
setwd("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen") # back to home base
ELP.chars$SUBTLWF <- as.numeric(gsub(",","",ELP.chars$SUBTLWF))
ELP.chars$BG_Mean <- as.numeric(gsub(",","",ELP.chars$BG_Mean))
ELP.chars$I_Mean_RT <- as.numeric(gsub(",","",ELP.chars$I_Mean_RT))
ELP.chars$logRT <- log10(ELP.chars$I_Mean_RT + 1)
ELP.chars$LgSUBTLWF <- as.numeric(log10(ELP.chars$SUBTLWF+1))
ELP.chars$OLD <- as.numeric(ELP.chars$OLD)


OSC.targs <- left_join(OSC.targs,ELP.chars,by=c("targ"="Word"))

# find morph family size for everything in OSC.targs
# add family size
require(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
mysheets <- read_excel_allsheets("~/Documents/MorphReadingDev/202003_LexiconProject/English/MorphoLEX_en.xlsx")
OSC.targs$famsize <- NA
for (i in 2:(length(mysheets)-3)){ # for each sheet in morpholex csv - note that first is intro sheet, last three morpheme sheets
  # if any words from oscM18 are in that sheet, add their morphological family size to oscM18
  toadd <- which(OSC.targs$targ %in% mysheets[[i]]$Word)
  if (length(toadd) > 0) {
    for (j in toadd) {
      OSC.targs$famsize[j] <- mysheets[[i]]$ROOT1_FamSize[which(mysheets[[i]]$Word==OSC.targs$targ[j])]
    }
  }
}
OSC.targs$logfamsize <- log10(OSC.targs$famsize+1)

# Ok, now I've calculated OSC for all the targs and elp words that I can - join features. 
oscM18 <- left_join(oscM18,OSC.targs,by=c("Word"="targ"))
rast.pem <- left_join(rast.pem,oscM18,by=c("targ"="Word"))

#################################
########## ELP DATA SECTION #####
#################################

elp.ana <- oscM18 %>% filter(numembed > 1) %>%
  select(c(Word,OSC,OSC.tcent,OSC.ccent,numembed,# BG_Mean,OLD,
           Length,LgSUBTLWF,logfamsize,logRT)) %>%
  filter(numembed > 1) %>% 
  drop_na()

# Are these predictors correlated? 
round(cor(elp.ana[,c(2,3,4,6,7,8)]),4)

baseline.model <- lm(logRT ~ Length + LgSUBTLWF + logfamsize, elp.ana)
summary(baseline.model)
osc.model <- lm(logRT ~ OSC + LgSUBTLWF + Length + logfamsize, elp.ana)
summary(osc.model)
osc.tcent.model <- lm(logRT ~ OSC.tcent + LgSUBTLWF + Length + logfamsize, elp.ana)
summary(osc.tcent.model)
osc.ccent.model <- lm(logRT ~ OSC.ccent + LgSUBTLWF + Length + logfamsize, elp.ana)
summary(osc.ccent.model)

# Does ccent explain variance in tcent model's residuals, or vice versa? 
elp.ana$resids <- residuals(osc.tcent.model)
summary(lm(resids ~ OSC.ccent, elp.ana)) # marginally
elp.ana$ccentresids <- residuals(osc.ccent.model)
summary(lm(ccentresids ~ OSC.tcent,elp.ana)) # very much so! 

# elp.ana %>% mutate(freq = round(LgSUBTLWF)) %>% 
#   ggplot(.,aes(x=OSC.ccent,y=logRT)) + 
#   geom_point() + geom_smooth(method="lm") + facet_grid(freq~.)

ic.table <- AIC(baseline.model,osc.model,osc.tcent.model,osc.ccent.model)
ic.table$BIC <- BIC(baseline.model,osc.model,osc.tcent.model,osc.ccent.model)$BIC
ic.table[order(ic.table$AIC),]
ic.table$model <- c("stem-centered, new values","stem-centered, original","baseline","cluster-centered")

# and last, if one of them is best, does the other explain additional variance? 

#################################
####### PRIMING DATA SECTION ####
#################################

rast.pem$prime.type <- factor(rast.pem$prime.type)
rast.pem$prime.type <- relevel(rast.pem$prime.type,"orthographic")

write.table(rast.pem,"primingdat_OSCcalcandotherprops_compare_updated_CELEXandorth.txt",quote = FALSE,row.names = FALSE)
# rast.pem <- read.table("primingdat_OSCcalcandotherprops_compare_updated_CELEXandorth.txt",header=TRUE)

#########################################################################################
######################### TARGET-FOCUSED ANALYSES #######################################
#########################################################################################

# - how well do these properties predict the target lexical decision times from ELP data? 
# do analysis, first using same base model of just fam size, frequency and length as in Marelli 2015, Marelli 2018 papers. 
targ.ana <- rast.pem %>% 
  select(c(name,prime.type,prime,unrel.RT,targ,study,OSC,OSC.tcent,OSC.ccent,numembed,SUBTLWF,LgSUBTLWF,
           Length,famsize,logfamsize,OLD,logRT)) %>% 
  filter(numembed > 1) %>% 
  drop_na() # 688 pairs, 547 *unique* pairs, 520 unique targets

### **** Does OSC differ significantly between conditions? ****
condition.effect <- lm(OSC ~ prime.type, targ.ana[-which(duplicated(targ.ana$targ)),])
summary(condition.effect)
condition.effect.tcent <- lm(OSC.tcent ~ prime.type, targ.ana[-which(duplicated(targ.ana$targ)),])
summary(condition.effect.tcent)
condition.effect.ccent <- lm(OSC.ccent ~ prime.type, targ.ana[-which(duplicated(targ.ana$targ)),])
summary(condition.effect.ccent)
# Yes, all three OSC metrics differ significantly across the conditions, in the direction you'd expect (greater OSC -> more transparent)

### **** How well does OSC predict target RTs from ELP, relative to transparency condition? ****
### (note that I changed this analysis, it used to be interaction between LgSUBTLWF and main 
### factor (primetype, OSC, etc) and OLD was another covariate)
RT.condition <- lm(logRT ~ prime.type + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(RT.condition)
RT.osc <- lm(logRT ~ OSC + LgSUBTLWF + logfamsize + Length , targ.ana)
summary(RT.osc)
RT.osc.tcent <- lm(logRT ~ OSC.tcent + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(RT.osc.tcent)
RT.osc.ccent <- lm(logRT ~ OSC.ccent + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(RT.osc.ccent)
ic.table <- AIC(RT.condition,RT.osc,RT.osc.tcent,RT.osc.ccent)
ic.table$BIC <- BIC(RT.condition,RT.osc,RT.osc.tcent,RT.osc.ccent)$BIC
ic.table[order(ic.table$AIC),] 
# original, then tcent are best, then ccent, then transparency condition 

### **** How well does OSC predict unrelated RTs, relative to transparency condition? ****
### (made a similar change here) 
targ.ana$lg.unrel.RT <- log10(targ.ana$unrel.RT + 1)
unrel.condition <- lm(lg.unrel.RT ~ prime.type + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(unrel.condition)
unrel.osc <- lm(lg.unrel.RT ~ OSC + LgSUBTLWF + logfamsize + Length , targ.ana)
summary(unrel.osc)
unrel.osc.tcent <- lm(lg.unrel.RT ~ OSC.tcent + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(unrel.osc.tcent)
unrel.osc.ccent <- lm(lg.unrel.RT ~ OSC.ccent + LgSUBTLWF + logfamsize + Length, targ.ana)
summary(unrel.osc.ccent)
ic.table <- AIC(unrel.condition,unrel.osc,unrel.osc.tcent,unrel.osc.ccent)
ic.table$BIC <- BIC(unrel.condition,unrel.osc,unrel.osc.tcent,unrel.osc.ccent)$BIC
ic.table[order(ic.table$AIC),] 
# tcent, then condition, then original, then ccent. Since main point is that stem-centered is better than cluster centered, this is fine. 

#########################################################################################
######################### PRIME-FOCUSED ANALYSES ########################################
#########################################################################################

# how well do these properties predict the priming magnitudes from set of studies? 

# add prime frequency - unless it's already there, then skip this. 
colnames(vocS) <- c("prime","primeSUBTLWF")
rast.pem <- left_join(rast.pem,vocS,by="prime")
rast.pem$primelogSUBTL <- log10(rast.pem$primeSUBTLWF + 1)

# add prime length - unless it's already there, then skip this. 
rast.pem$primelength <- sapply(rast.pem$prime,str_length)

# add affix family size. 
take.aff <- function(x){
  if (is.na(x)) {return(NA)
  } else {
    return(gsub("[^a-z]","",strsplit(x,split="[>]")[[1]][2]))
  }
}
suff.df <- read.table("suffdf.txt")
onestem.onesuff.df <- read.table("onestemonesuffdf.txt")
rast.pem <- left_join(rast.pem,onestem.onesuff.df[,c(2,6)],by=c("prime" = "Word"))
rast.pem$aff <- sapply(rast.pem$MorphoLexSegm,take.aff)

# go through the ones that don't have MorphoLex entries and plug in the affixes
aff.ops <- unique(sapply(onestem.onesuff.df$MorphoLexSegm,take.aff),suff.df$morpheme)
rast.pem$aff.fam.size <- NA
for (i in 1:nrow(rast.pem)) {
  # try to figure out the affix of the cxword without entering manually
  if (is.na(rast.pem$aff[i]) & gsub(rast.pem$targ[i],"", rast.pem$prime[i]) %in% aff.ops){
    rast.pem$aff[i] <- gsub(rast.pem$targ[i],"", rast.pem$prime[i])
  } else if (is.na(rast.pem$aff[i])) {print(paste("no aff for",rast.pem$prime[i]))}
  if (rast.pem$aff[i] %in% suff.df$morpheme) {
    rast.pem$aff.fam.size[i] <- suff.df$family_size[which(suff.df$morpheme==rast.pem$aff[i])]
  }
}


# Add in a transparency metric - how semantically similar are the prime and the target?
# UNLESS transparency is already there, in which case skip this
rast.pem$trans <- NA
for (ri in 1:nrow(rast.pem)) {
  if (rast.pem$prime[ri] %in% glovevecs$Word & rast.pem$targ[ri] %in% glovevecs$Word) {
    rast.pem$trans[ri] <- cos.sim(glovevecs[glovevecs$Word==rast.pem$prime[ri],2:(ncol(glovevecs)-2)],
                                  glovevecs[glovevecs$Word==rast.pem$targ[ri],2:(ncol(glovevecs)-2)])
  } else {print(ifelse(rast.pem$prime[ri] %in% glovevecs$Word,rast.pem$targ[ri],rast.pem$prime[ri]))}
}

setwd("~/Documents/MorphReadingDev/202012_OSCemergence/stimulusGen")
write.table(rast.pem,"primingdat_OSCcalcandotherprops_compare_updated_CELEXandorth.txt",quote = FALSE,row.names = FALSE)
# rast.pem <- read.table("primingdat_OSCcalcandotherprops_compare_updated_CELEXandorth.txt",header=TRUE)

require(lme4)
require(lmerTest)
require(MuMIn)

## FIRST VERSION : INCLUDING orthographic condition: 
prime.ana <- rast.pem %>% 
  select(name,prime.type,PEM,prime,targ,study,OSC,OSC.tcent,OSC.ccent,
         numembed,LgSUBTLWF,primelogSUBTL,Length,primelength,OLD,
         logfamsize,trans) %>% drop_na()
# 676 items, 513 unique targets, 535 unique prime-target pairs

base.mod <- lmer(PEM ~ primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana) # random intercepts for pair or target, inevitably led to singular model fit. 
summary(base.mod)
osc.mod <- lmer(PEM ~ OSC*primelogSUBTL + OSC*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
# osc.mod <- lmer(PEM ~ OSC + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.mod)
osc.tcent.mod <- lmer(PEM ~ OSC.tcent*primelogSUBTL + OSC.tcent*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
# osc.tcent.mod <- lmer(PEM ~ OSC.tcent + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.tcent.mod)
osc.ccent.mod <- lmer(PEM ~ OSC.ccent*primelogSUBTL + OSC.ccent*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
# osc.ccent.mod <- lmer(PEM ~ OSC.ccent + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.ccent.mod)
ic.table <- AIC(base.mod,osc.mod,osc.tcent.mod,osc.ccent.mod)
ic.table$BIC <- BIC(base.mod,osc.mod,osc.tcent.mod,osc.ccent.mod)$BIC
ic.table[order(ic.table$AIC),] # ccent is best, then orig, tcent, base. 

r.squaredGLMM(base.mod)
r.squaredGLMM(osc.mod)
r.squaredGLMM(osc.tcent.mod)
r.squaredGLMM(osc.ccent.mod)

## SECOND VERSION: EXCLUDING ORTHOGRAPHIC CONDITION: 
prime.ana <- rast.pem %>% 
  filter(prime.type != "orthographic") %>% 
  select(name,prime.type,PEM,prime,targ,study,OSC,OSC.tcent,OSC.ccent,
         numembed,LgSUBTLWF,primelogSUBTL,Length,primelength,OLD,
         logfamsize,trans) %>% drop_na() # not using affix family size for now, for analyses WITH orthographic 
# 480 items, 385 unique targets, 397 unique prime-target pairs
base.mod <- lmer(PEM ~ primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana) # random intercepts for pair or target, inevitably led to singular model fit. 
summary(base.mod)
# osc.mod <- lmer(PEM ~ OSC*primelogSUBTL + OSC*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
osc.mod <- lmer(PEM ~ OSC + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.mod)
# osc.tcent.mod <- lmer(PEM ~ OSC.tcent*primelogSUBTL + OSC.tcent*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
osc.tcent.mod <- lmer(PEM ~ OSC.tcent + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.tcent.mod)
# osc.ccent.mod <- lmer(PEM ~ OSC.ccent*primelogSUBTL + OSC.ccent*LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
osc.ccent.mod <- lmer(PEM ~ OSC.ccent + primelogSUBTL + LgSUBTLWF + Length + logfamsize + OLD + trans + (1|targ) + (1|study), prime.ana)
summary(osc.ccent.mod)
ic.table <- AIC(base.mod,osc.mod,osc.tcent.mod,osc.ccent.mod)
ic.table$BIC <- BIC(base.mod,osc.mod,osc.tcent.mod,osc.ccent.mod)$BIC
ic.table[order(ic.table$AIC),] # original, then ccent, then tcent.

# get R^2 values for these models - for WoW poster table comparison. 
# install.packages("MuMIn")
require(MuMIn)
r.squaredGLMM(base.mod)
r.squaredGLMM(osc.mod)
r.squaredGLMM(osc.tcent.mod)
r.squaredGLMM(osc.ccent.mod)

require(ggplot2)
quant.OSC <- quantile(prime.ana$OSC.ccent)
get.quant <- function(x,qs){return(min(which(qs >= x))-1)}
prime.ana$osc.q <- sapply(prime.ana$OSC.ccent,get.quant,qs=quant.OSC)
prime.ana %>% 
  ggplot(.,aes(x=LgSUBTLWF,y=PEM,col=factor(osc.q))) + 
  ggtitle("cluster-centered OSC x targfreq on PEM, with orthographic condition") + 
 # geom_point() + theme_light() + 
  geom_smooth(method="lm",se=FALSE)# (method="lm",se=FALSE,alpha=0.7) # Come back to this - figure out the error. 
ggplot(prime.ana, aes(x=OSC.ccent,y=PEM)) + geom_point() + geom_smooth(method="lm")

