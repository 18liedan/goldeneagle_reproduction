setwd("/Users/amnl/Desktop/mhcwork/forR")
library("dplyr")

data <- read.csv("alleles.csv", header = TRUE)
alleles <- subset(data, select = -c(allele))

data_wild <- read.csv("alleleswild.csv", header=TRUE)
alleles_wild <- subset(data_wild, select = -c(allele))

##MHCshare###
shared <- data.frame(matrix(NA,ncol=16,nrow=1))
pairnames <- c("dewasora","nishimechokai","nishimefuuto","nishimeshiun","tsukikosenshu",
               "tsukikofuuto","komachioume","tatsukoshinano","fuukikonara","tokihimekami",
               "okinasora","okinafukui","hakusankake","fukuiagano","shionagano", "koumeagano")
colnames(shared) <- pairnames

shared$dewasora <- length(intersect(na.omit(alleles$X30031dewa), na.omit(alleles$X30032sora)))
shared$nishimechokai <- length(intersect(na.omit(alleles$X12774nishime), na.omit(alleles$X12775chokai)))
shared$nishimefuuto <- length(intersect(na.omit(alleles$X12774nishime), na.omit(alleles$X31237fuuto)))
shared$nishimeshiun <- length(intersect(na.omit(alleles$X12774nishime), na.omit(alleles$X31071shiun)))
shared$tsukikosenshu <- length(intersect(na.omit(alleles$X31232tsukiko), na.omit(alleles$X12776senshu)))
shared$tsukikofuuto<- length(intersect(na.omit(alleles$X31232tsukiko), na.omit(alleles$X31237fuuto)))
shared$komachioume <- length(intersect(na.omit(alleles$X14474komachi), na.omit(alleles$X14475oume)))
shared$tatsukoshinano <- length(intersect(na.omit(alleles$X31233tatsuko), na.omit(alleles$X31234shinano)))
shared$fuukikonara <- length(intersect(na.omit(alleles$X30949fuuki), na.omit(alleles$X30950konara)))
shared$tokihimekami <- length(intersect(na.omit(alleles$X31060toki), na.omit(alleles$X31061himekami)))
shared$okinasora <- length(intersect(na.omit(alleles$X29330okina), na.omit(alleles$X30032sora)))
shared$okinafukui <- length(intersect(na.omit(alleles$X29330okina), na.omit(alleles$X31914fukui)))
shared$hakusankake <- length(intersect(na.omit(alleles$X31250hakusan), na.omit(alleles$X31252kakehashi)))
shared$fukuiagano <- length(intersect(na.omit(alleles$X31914fukui),na.omit(alleles$X31915agano)))
shared$shionagano <- length(intersect(na.omit(alleles$X31076shion),na.omit(alleles$X31915agano)))
shared$koumeagano <- length(intersect(na.omit(alleles$X31064koume),na.omit(alleles$X31915agano)))

write.csv(shared,"/Users/amnl/Desktop/mhcwork/forR/shared_out.csv")

###Pairwise Difference###

##Trial##
#Vab: number of different alleles between individuals a and b
length(na.omit(alleles$X12774nishime))+length(na.omit(alleles$X12775chokai))-
  2*(length(intersect(na.omit(alleles$X12774nishime), na.omit(alleles$X12775chokai))))

#Fa: total number of alleles carried by individual a
length(na.omit(alleles$X12774nishime))

#Fb: total number of alleles carried by individual b
length(na.omit(alleles$X12775chokai))

#PD = 100*Vab/(Fa+Fb) percent difference of alleles between individual a and b 
100*(length(na.omit(alleles$X12774nishime))+length(na.omit(alleles$X12775chokai))-
       2*(length(intersect(na.omit(alleles$X12774nishime), na.omit(alleles$X12775chokai)))))/
  (length(na.omit(alleles$X12774nishime))+length(na.omit(alleles$X12775chokai)))


##apply to all columns##
#make PD function
PD <- function(a, b) {
  Fa <- length(na.omit(a))
  Fb <- length(na.omit(b))
  shared <- length(intersect(na.omit(a),na.omit(b)))
  Vab <- (Fa+Fb)-(2*shared)
  PD <- 100*Vab/(Fa+Fb)
  PD
}

apd_out <- as.data.frame(matrix(nrow=51,ncol=51))
colnames(apd_out) = c(colnames(alleles))

for (n in 1:51) {
  apd_out[ ,n] <- as.vector((apply(alleles, 2, PD, a = alleles[ ,n])))
}

write.csv(apd_out,"/Users/amnl/Desktop/mhcwork/forR/apd_out.csv")

#apply runs function across column (MARGIN = 2) so specifying one variable is enough? 
#i.e. by default, either a or b are already designated as the value in the columns?
#and then the specified variable (in this case, each column) will be the 
#other variable that is being compared?

#do for wild
apd_out_wild <- as.data.frame(matrix(nrow=20,ncol=20))
colnames(apd_out_wild) = c(colnames(alleles_wild))

for (n in 1:20) {
  apd_out_wild[ ,n] <- as.vector((apply(alleles_wild, 2, PD, a = alleles_wild[ ,n])))
}

write.csv(apd_out_wild,"/Users/amnl/Desktop/mhcwork/forR/apd_out_wild.csv")

###AAdist###
library("dplyr")

##count the number of differences between each site per allele
aa <- data.frame(read.csv("aavsites.csv", colClasses="character", header = TRUE))
aa1 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa2 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa3 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa4 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa5 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa6 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa7 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa8 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa9 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa10 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa12 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa13 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa14 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa15 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa16 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa17 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]
aa19 <- data.frame(matrix(NA,ncol=1,nrow=17))[-1]

for (i in 1:17) { #17 alleles,
  for (j in 1:17)  { #17 variable sites
    aa1[i,j] <- all(aa[1,j] == aa[i,j])
    aa2[i,j] <- all(aa[2,j] == aa[i,j])
    aa3[i,j] <- all(aa[3,j] == aa[i,j])
    aa4[i,j] <- all(aa[4,j] == aa[i,j])
    aa5[i,j] <- all(aa[5,j] == aa[i,j])
    aa6[i,j] <- all(aa[6,j] == aa[i,j])
    aa7[i,j] <- all(aa[7,j] == aa[i,j])
    aa8[i,j] <- all(aa[8,j] == aa[i,j])
    aa9[i,j] <- all(aa[9,j] == aa[i,j])
    aa10[i,j] <- all(aa[10,j] == aa[i,j])
    aa12[i,j] <- all(aa[11,j] == aa[i,j])
    aa13[i,j] <- all(aa[12,j] == aa[i,j])
    aa14[i,j] <- all(aa[13,j] == aa[i,j])
    aa15[i,j] <- all(aa[14,j] == aa[i,j])
    aa16[i,j] <- all(aa[15,j] == aa[i,j])
    aa17[i,j] <- all(aa[16,j] == aa[i,j])
    aa19[i,j] <- all(aa[17,j] == aa[i,j])
  }
}

##make aadist matrix for pbr
aadist_pbr <- data.frame(matrix(NA, ncol=1, nrow=17))[-1]
aadist_pbr$V1 <- 15-rowSums(aa1[1:15])
aadist_pbr$V2 <- 15-rowSums(aa2[1:15])
aadist_pbr$V3 <- 15-rowSums(aa3[1:15])
aadist_pbr$V4 <- 15-rowSums(aa4[1:15])
aadist_pbr$V5 <- 15-rowSums(aa5[1:15])
aadist_pbr$V6 <- 15-rowSums(aa6[1:15])
aadist_pbr$V7 <- 15-rowSums(aa7[1:15])
aadist_pbr$V8 <- 15-rowSums(aa8[1:15])
aadist_pbr$V9 <- 15-rowSums(aa9[1:15])
aadist_pbr$V10 <- 15-rowSums(aa10[1:15])
aadist_pbr$V12 <- 15-rowSums(aa12[1:15])
aadist_pbr$V13 <- 15-rowSums(aa13[1:15])
aadist_pbr$V14 <- 15-rowSums(aa14[1:15])
aadist_pbr$V15 <- 15-rowSums(aa15[1:15])
aadist_pbr$V16 <- 15-rowSums(aa16[1:15])
aadist_pbr$V17 <- 15-rowSums(aa17[1:15])
aadist_pbr$V19 <- 15-rowSums(aa19[1:15])

write.csv(aadist_pbr,"/Users/amnl/Desktop/mhcwork/forR/aadist_pbr.csv")

##make aadist matrix or nonpbr
aadist_non_pbr <- data.frame(matrix(NA, ncol=1, nrow=17))[-1]
aadist_non_pbr$V1 <- 2-rowSums(aa1[16:17])
aadist_non_pbr$V2 <- 2-rowSums(aa2[16:17])
aadist_non_pbr$V3 <- 2-rowSums(aa3[16:17])
aadist_non_pbr$V4 <- 2-rowSums(aa4[16:17])
aadist_non_pbr$V5 <- 2-rowSums(aa5[16:17])
aadist_non_pbr$V6 <- 2-rowSums(aa6[16:17])
aadist_non_pbr$V7 <- 2-rowSums(aa7[16:17])
aadist_non_pbr$V8 <- 2-rowSums(aa8[16:17])
aadist_non_pbr$V9 <- 2-rowSums(aa9[16:17])
aadist_non_pbr$V10 <- 2-rowSums(aa10[16:17])
aadist_non_pbr$V12 <- 2-rowSums(aa12[16:17])
aadist_non_pbr$V13 <- 2-rowSums(aa13[16:17])
aadist_non_pbr$V14 <- 2-rowSums(aa14[16:17])
aadist_non_pbr$V15 <- 2-rowSums(aa15[16:17])
aadist_non_pbr$V16 <- 2-rowSums(aa16[16:17])
aadist_non_pbr$V17 <- 2-rowSums(aa17[16:17])
aadist_non_pbr$V19 <- 2-rowSums(aa19[16:17])

write.csv(aadist_non_pbr,"/Users/amnl/Desktop/mhcwork/forR/aadist_nonpbr.csv")

##aadist per pair
library(data.table)
library(tidyverse)
library(tidytable)

data <- read.csv("alleles.csv", header = TRUE)
alleles <- subset(data, select = -c(allele))

dewa30031 <- as.vector(na.omit(alleles$X30031dewa))
sora30032 <- as.vector(na.omit(alleles$X30032sora))
nishime12774 <- as.vector(na.omit(alleles$X12774nishime))
chokai12775 <- as.vector(na.omit(alleles$X12775chokai))
fuuto31237 <- as.vector(na.omit(alleles$X31237fuuto))
shiun31071 <- as.vector(na.omit(alleles$X31071shiun))
tsukiko31232 <- as.vector(na.omit(alleles$X31232tsukiko))
senshu12776 <- as.vector(na.omit(alleles$X12776senshu))
okina29330  <- as.vector(na.omit(alleles$X29330okina))
komachi14474 <- as.vector(na.omit(alleles$X14474komachi))
oume14475 <- as.vector(na.omit(alleles$X14475oume))
tatsuko31233 <- as.vector(na.omit(alleles$X31233tatsuko))
shinano31234 <- as.vector(na.omit(alleles$X31234shinano))
fuuki30949 <- as.vector(na.omit(alleles$X30949fuuki))
konara30950 <- as.vector(na.omit(alleles$X30950konara))
toki31060 <- as.vector(na.omit(alleles$X31060toki))
himekami31061 <- as.vector(na.omit(alleles$X31061himekami))
hakusan31250 <- as.vector(na.omit(alleles$X31250hakusan))
kakehashi31252 <- as.vector(na.omit(alleles$X31252kakehashi))
agano31915 <- as.vector(na.omit(alleles$X31915agano))
fukui31914 <- as.vector(na.omit(alleles$X31914fukui))
shion31076 <- as.vector(na.omit(alleles$X31076shion))
koume31064 <- as.vector(na.omit(alleles$X31064koume))

dewasora <- as.data.frame(crossing.(dewa30031, sora30032))
dewsor<-paste(dewasora$dewa30031,dewasora$sora30032, sep = ",")
nishimechokai <- as.data.frame(crossing.(nishime12774, chokai12775))
nischo <- paste(nishimechokai$nishime12774,nishimechokai$chokai12775, sep = ",")
nishimefuuto <- as.data.frame(crossing.(nishime12774, fuuto31237))
nisfuu<-paste(nishimefuuto$nishime12774,nishimefuuto$fuuto31237, sep = ",")
nishimeshiun <- as.data.frame(crossing.(nishime12774,shiun31071))
nisshi<-paste(nishimeshiun$nishime12774,nishimeshiun$shiun31071, sep = ",")
tsukikosenshu <- as.data.frame(crossing.(tsukiko31232, senshu12776))
tsusen<-paste(tsukikosenshu$tsukiko31232,tsukikosenshu$senshu12776, sep = ",")
tsukikofuuto <- as.data.frame(crossing.(tsukiko31232, fuuto31237))
tsufuu<-paste(tsukikofuuto$tsukiko31232,tsukikofuuto$fuuto31237, sep = ",")
komachioume <- as.data.frame(crossing.(komachi14474, oume14475))
komoum<-paste(komachioume$komachi14474,komachioume$oume14475, sep = ",")
tatsukoshinano <- as.data.frame(crossing.(tatsuko31233, shinano31234))
tatshi<-paste(tatsukoshinano$tatsuko31233,tatsukoshinano$shinano31234, sep = ",")
fuukikonara <- as.data.frame(crossing.(fuuki30949, konara30950))
fuukon<-paste(fuukikonara$fuuki30949,dewasora$konara30950, sep = ",")
tokihimekami <- as.data.frame(crossing.(toki31060, himekami31061))
tokhim<-paste(tokihimekami$toki31060,tokihimekami$himekami31061, sep = ",")
okinasora <- as.data.frame(crossing.(okina29330, sora30032))
okisor<-paste(okinasora$okina29330, okinasora$sora30032, sep = ",")
okinafukui <- as.data.frame(crossing.(okina29330, fukui31914))
okifuk<-paste(okinafukui$okina29330, okinafukui$fukui31914, sep = ",")
hakusankake <- as.data.frame(crossing.(hakusan31250, kakehashi31252))
hakkak<-paste(hakusankake$hakusan31250, hakusankake$kakehashi31252, sep = ",")
fukuiagano <- as.data.frame(crossing.(fukui31914, agano31915))
fukaga<-paste(fukuiagano$fukui31914, fukuiagano$agano31915, sep=",")
shionagano <- as.data.frame(crossing.(shion31076,agano31915))
shiaga<-paste(shionagano$shion31076, shionagano$agano31915, sep=",")
koumeagano <- as.data.frame(crossing.(koume31064, agano31915))
kouaga<-paste(koumeagano$koume31064, koumeagano$agano31915, sep=",")

addToDf <- function(df, df2){
  nRow <- nrow(df)
  lngth <- nrow(df2)
  if(nRow > lngth){
    df2[(lngth+1):nRow,] <-NA
  }else if(nRow < lngth){
    df[(nRow+1):lngth,] <- NA
  }
  cbind(df,df2)
}

pairalleles <- data.frame(dewasora)
pairalleles <- addToDf(pairalleles, nishimechokai)
pairalleles <- addToDf(pairalleles, nishimefuuto)
pairalleles <- addToDf(pairalleles, nishimeshiun)
pairalleles <- addToDf(pairalleles, tsukikosenshu)
pairalleles <- addToDf(pairalleles, tsukikofuuto)
pairalleles <- addToDf(pairalleles, komachioume)
pairalleles <- addToDf(pairalleles, tatsukoshinano)
pairalleles <- addToDf(pairalleles, fuukikonara)
pairalleles <- addToDf(pairalleles, tokihimekami)
pairalleles <- addToDf(pairalleles, okinasora)
pairalleles <- addToDf(pairalleles, okinafukui)
pairalleles <- addToDf(pairalleles, hakusankake)
pairalleles <- addToDf(pairalleles, fukuiagano)
pairalleles <- addToDf(pairalleles, shionagano)
pairalleles <- addToDf(pairalleles, koumeagano)

write.csv(pairalleles,"/Users/amnl/Desktop/mhcwork/forR/allelecombinations.csv")
