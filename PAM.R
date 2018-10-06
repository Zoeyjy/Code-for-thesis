library(R2WinBUGS)


# Data reading and preparation

data<-read.csv("DataSA.csv", header=TRUE)
head(data)
dim(data)
data$Date<-seq(as.Date("2011/01/01"), as.Date("2012/12/31"), by = "day")
head(data)
summary(data)
data$Date[which(data$CWOD==0)] #"2012-01-03" "2012-01-21"
data$Date[which(data$WNH4==0)] #"2011-10-08" "2011-10-09" "2011-10-10"
data$Date[which(data$TAL==0)] # "2012-12-30"

data$Date[which(data$TCA<0)] #"2011-08-13" "2012-06-11" "2012-07-15" "2012-12-25"
data$Date[which(data$TCA==0)] # no observation

data$Date[which(data$TCU<0)] #"2012-04-29" "2012-07-15"
data$Date[which(data$TCU==0)] # no observation

data$Date[which(data$TNI==0)] #"2012-12-25"
data$Date[which(data$TTI==0)] #"2012-04-07"
data$Date[which(data$TV==0)] #"2012-12-25"

# There are values rounded to zero (i.e values that are reported as zero) 
# or values below detection limit (the value is reported to be zero as the detection method could not prove its presence)
# Then, We apply suggestion from the book "Modeling and Analysis of Compositional Data" page 17 (Vera Pawlowsky-Glahn,Juan José Egozcue,Raimon Tolosana-Delgado)
# and we replace zero with something small, below or around the rounding/detection limit (see also page 20).

# Martin-Fernandez et al 2003 suggested taking 1/2 or 2/3 of the detection limit 
# J. A. Martín-FernándezC. Barceló-VidalV. Pawlowsky-Glahn (2003) Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. Mathematical Geology, 35:253-278

dataNew <- data

summary(data$CWOD[which(data$CWOD>0)]) # 0.100
dataNew$CWOD[which(data$CWOD==0)] <- min(data$CWOD[which(data$CWOD>0)])*0.5

summary(data$WNH4[which(data$WNH4>0)]) # min=0.02
dataNew$WNH4[which(data$WNH4==0)] <- min(data$WNH4[which(data$WNH4>0)])*0.5

summary(data$TAL[which(data$TAL>0)]) # min=0.00158
dataNew$TAL[which(data$TAL==0)] <- min(data$TAL[which(data$TAL>0)])*0.5

summary(data$TCA[which(data$TCA>0)]) # min=0.00159
dataNew$TCA[which(data$TCA<0)] <- min(data$TCA[which(data$TCA>0)])*0.5

summary(data$TCU[which(data$TCU>0)]) # min=0.000260
dataNew$TCU[which(data$TCU<0)] <- min(data$TCU[which(data$TCU>0)])*0.5

summary(data$TNI[which(data$TNI>0)]) # min=0.000030
dataNew$TNI[which(data$TNI==0)] <- min(data$TNI[which(data$TNI>0)])*0.5

summary(data$TTI[which(data$TTI>0)]) # min=0.000810
dataNew$TTI[which(data$TTI==0)] <- min(data$TTI[which(data$TTI>0)])*0.5

summary(data$TV[which(data$TV>0)]) # min=0.00106
dataNew$TV[which(data$TV==0)] <- min(data$TV[which(data$TV>0)])*0.5

summary(dataNew)
# standardization z ((x-m)/s) # possible problem with closure....
#dataZ<-decostand(air.sp, "standardize", na.rm=TRUE)
#head(dataZ)
#dataZ<-round(dataZ, digits=5)
#dataZ<-as.matrix(dataZ)



air.sp<-dataNew[,2:24]
head(air.sp)
summary(air.sp)
air.cor<-round(cor(air.sp, use="pairwise.complete.obs"),3)
air.cor

# log transformation of the data

air.spLog <- round(log(air.sp),5)
head(air.spLog)

#############
### Model ###
#############

T = nrow(air.spLog)
nS = ncol(air.spLog)
m = 10
S<- solve(round(cov(air.spLog, use="pairwise.complete.obs"),5))

summary(wind)
min1 <- min(wind$u, na.rm=TRUE)
max1 <- max(wind$u, na.rm=TRUE)
min2 <- min(wind$v, na.rm=TRUE)
max2 <- max(wind$v, na.rm=TRUE)

#R0=diag(nS) #identity matrix
R0=round(solve(cor(air.spLog, use="pairwise.complete.obs")),5)

datinits <- array(NA, c(T,nS))
datinitsdf = data.frame(datinits)
colnames(datinitsdf) <- colnames(air.spLog)
head(datinitsdf)



medianS<- as.vector(apply(air.spLog, 2, quantile, probs = c(0.5), na.rm=TRUE))
medianS2 <- t(rbind(do.call(cbind, rep(list(medianS), T))))
dim(medianS2)
datinitsdf[is.na(air.spLog)]  <- medianS2[which(is.na(air.spLog))]
dim(datinitsdf) 
head(datinitsdf)

datinitsM<- as.matrix(datinitsdf)

wind$Date[which(is.na(wind$u))]
summary(wind)
windAug<- subset(wind, wind$Date >= as.Date("2011-08-01") & wind$Date <= as.Date("2011-08-31"))
head(windAug)
tail(windAug)
summary(windAug)

wind2<-wind
wind2$u[is.na(wind$u)] <- 1.4249 # this is the median for u of August
wind2$v[is.na(wind$v)] <- 0.8137 # this is the median for v of August

air.spLogM <- as.matrix(air.spLog)


vm=m-1
z1=sample(seq(1,m),T,replace=TRUE)
z2=sample(seq(1,m),T,replace=TRUE)

#############################
##### MODEL  KERNEL SBP #####
#############################
# NOTE: I need define all the inits or it crashes
# see http://esapubs.org/archive/appl/A025/026/Supplement1_WinBUGS_code.R

DataModel<-list("T"=T,"nS"=nS, "m"=m, "species"=air.spLogM, "wu"=wind2$u, "wv"=wind2$v, "R0"=R0, "min1"=min1, "max1"=max1 ,"min2"=min2, "max2"=max2, "lmax"=731, ndf=24)

inits <- list(
list(int=0.1,  Taue=rep(1,nS), species=datinitsM, v=runif(vm, 0,1), z=z1,
     invbw = c(0.49,0.20,1.28,0.63,0.73,0.45,0.29,0.11,0.12),
     knot1 = c( -1.19,3.49,1.13,-0.41,1.49,3.11,4.66,0.88,-3.07),
     knot2 = c(0.01,3.88,-2.89,0.015,1.14,-2.33,-2.31,0.49,-2.39),
     lambda = 3.06, v = c(0.40,0.75,0.22,0.77,0.01,0.03,0.02,0.86,0.25)), 
list(int=0.5,  Taue=rep(2,nS), species=datinitsM, v=runif(vm, 0,1), z=z2,
     invbw = c(1.56,0.03,1.07,0.79,0.96,6.63,1.20,0.33,0.05),
     knot1 = c(-0.38,4.34,4.39,-0.57,4.59,4.52,-2.03,-0.34,-0.57),
     knot2 = c(0.07,-0.20,3.58,0.01,-3.18,2.40,3.14,0.02,-2.24),
     lambda = 0.55, v = c(0.64,0.09,0.26,0.98,0.31,0.83,0.37,0.99,0.48)))

parameters=c("theta","int","K", "z", "v","p", "mu.Sp", "lambda", "invbw")
  
modelsim <- bugs(DataModel, inits, parameters, model.file="MODEL_KSBP.txt",
                   n.chains = 2, n.iter = 25000, n.burnin = 23000, n.thin=2,
                   debug=FALSE, DIC=FALSE, 
                   bugs.directory="C:/winbugs14/WinBUGS14", 
                   working.directory =  "C:/SA_Study/KSBP", bugs.seed=112510)
#n.iter = 10000, n.burnin = 8000
#n.iter = 20000, n.burnin = 18000
#n.iter = 22000, n.burnin = 20000

write.csv(modelsim$sims.matrix,paste0("modpar",".csv"))




# matrix of z
clmat <- modelsim$sims.list$z
dim(clmat)
# score matrix
matM <- matrix(0,731,731)
matMall <- matrix(0,731,731)


for (b  in 1:2000){
    for (a in 1:731){
        rowa <- replicate(731, clmat[b,a])
        matM[a,which(rowa==clmat[b,])] <- 1
    }
    matMall <- matMall + matM
}

# probability matrix (Similarity matrix)
matS <- matMall/2000
# dissimilarity matrix
matDS <- 1-matS


######## Partitioning around matrix
library(cluster)
getmode(modelsim$sims.list$K)#  5 clusters as mode

partition <- pam(matDS, 5, diss=TRUE, keep.diss = FALSE)

cluster <- partition$clustering

length(which(cluster==1))  #268
length(which(cluster==2))  # 91
length(which(cluster==3))  # 193
length(which(cluster==4))  # 178
length(which(cluster==5))  # 1
which(cluster ==5)    #368

# switch cluster 1 and cluster 2
cluster[which(cluster==1)] <- 12
cluster[which(cluster==2)] <- 1
cluster[which(cluster==12)] <- 2

length(which(cluster==1))  #91
length(which(cluster==2))  # 268


mu.sp <- modelsim$sims.list$mu.Sp
# split mu.sp for clusters
mucl1 <- mu.sp[,which(cluster == 1),]  # 2000*268*23
mucl2 <- mu.sp[,which(cluster == 2),]  #
mucl3 <- mu.sp[,which(cluster == 3),]  #
mucl4 <- mu.sp[,which(cluster == 4),]  #
mucl5 <- mu.sp[,which(cluster == 5),]  #

## cluster 1
ave1 <- apply(mucl1,3,rowMeans)
## cluster 2
ave2 <- apply(mucl2,3,rowMeans)
## cluster 3
ave3 <- apply(mucl3,3,rowMeans)
## cluster 4
ave4 <- apply(mucl4,3,rowMeans)

####### boxplot
namecol <- c("EC","OC","CWOD","WNO3","WSO4","WCL","WNH4","TAL","TBA","TCA",	"TCU",	"TFE","TK","TMN","TMO",	"TNA","TNI","TPB","TSB","TTI","TV",	"TZN","TMG")

library(reshape2)
colnames(ave1) <- namecol
new1 <- melt(ave1)
colnames(new1) <- c("iteration","comp","value")

colnames(ave2) <- namecol
new2 <- melt(ave2)
colnames(new2) <- c("iteration","comp","value")

colnames(ave3) <- namecol
new3 <- melt(ave3)
colnames(new3) <- c("iteration","comp","value")

colnames(ave4) <- namecol
new4 <- melt(ave4)
colnames(new4) <- c("iteration","comp","value")

colnames(mucl5) <- namecol
new5 <- melt(mucl5)
colnames(new5) <- c("iteration","comp","value")


library(ggplot2)

pdf("pam25.pdf")

ggplot(new1, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 1")

ggplot(new2, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 2")

ggplot(new3, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 3")

ggplot(new4, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 4")

ggplot(new5, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 5")

dev.off()

############## Back transform

org1 <- new1
org1$value <- exp(org1$value)

org2 <- new2
org2$value <- exp(org2$value)

org3 <- new3
org3$value <- exp(org3$value)

org4 <- new4
org4$value <- exp(org4$value)

org5 <- new5
org5$value <- exp(org5$value)

library(ggplot2)

pdf("pam25_originalscale.pdf")

ggplot(org1, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 1")+
coord_cartesian(ylim = c(0, 10))

ggplot(org2, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 2") +
coord_cartesian(ylim = c(0, 10))

ggplot(org3, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 3")+
coord_cartesian(ylim = c(0, 10))

ggplot(org4, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 4")+
coord_cartesian(ylim = c(0, 10))

ggplot(org5, aes(x=comp, y=value, fill=comp)) +
geom_boxplot(alpha=0.3) +
theme(legend.position="none")+
xlab("PM10 components") + ylab("Posterior estimates")+
labs(title = "Cluster 5")+
coord_cartesian(ylim = c(0, 10))


dev.off()


# quantile(o5[0.20])



########## determine the partition, sulinfo
library(cluster)
unique(modelsim$sims.list$K)  # 5, 4, 6
getmode(modelsim$sims.list$K)#  5 clusters as mode

partition <- pam(matDS, 5, diss=TRUE, keep.diss = FALSE)
partition$silinfo$avg.width    #0.998632

par4 <- pam(matDS, 4, diss=TRUE, keep.diss = FALSE)
par4$silinfo$avg.width   #0.9933249


par6 <- pam(matDS, 6, diss=TRUE, keep.diss = FALSE)
par6$silinfo$avg.width    #0.75513


#### correlation matrix

upper<- air.cor
upper[upper.tri(air.cor,diag=TRUE)]<-""
upper<-as.data.frame(upper)
upper

library(xtable)
xtable(upper)

##### summary statistics
statsum <- as.data.frame(do.call(rbind, lapply(air.sp, summary)))
statsum <- round(statsum, 4)
statsum$range <- do.call(paste, c(statsum[c("Min.", "Max.")], sep = "-"))

statsumfinal <- as.data.frame(cbind(statsum$Mean, statsum$range, statsum$`1st Qu.`,statsum$Median,statsum$`3rd Qu.`))
rownames(statsumfinal) <- namecol
colnames(statsumfinal) <- c("Mean","Range","25th","50th","75th")
xtable(statsumfinal)

# cluster summary statistics  original scale
# ave1, ave2, ave3, ave4, org5
o1 <- exp(ave1)
o2 <- exp(ave2)
o3 <- exp(ave3)
o4 <- exp(ave4)
o5 <- exp(mucl5)


mean(o1[,1])
mean(o1[,1]) - qnorm(0.975)*sd(o1[,1])/sqrt(2000)
mean(o1[,1]) + qnorm(0.975)*sd(o1[,1])/sqrt(2000)
library(Publish)
ci.mean(o1[,1])


mci <- function(x){
    ciupper <- round(apply(x, 2, function(y) mean(y) + qnorm(0.975)*sd(y)/sqrt(2000)),4)
    cilower <- round(apply(x, 2, function(y) mean(y) - qnorm(0.975)*sd(y)/sqrt(2000)),4)
    mu <- round(apply(x, 2,mean),4)
    c1 <- paste(cilower, ciupper, sep = ", ")
    c2 <- paste(replicate(23,"("),c1,replicate(23,")"))
    paste(mu,c2)
}

mci1 <- as.data.frame(mci(o1))
mci2 <- as.data.frame(mci(o2))
mci3 <- as.data.frame(mci(o3))
mci4 <- as.data.frame(mci(o4))
mci5 <- as.data.frame(mci(o5))
confi <- cbind(mci1,mci2,mci3,mci4,mci5)
rownames(confi) <- namecol
colnames(confi) <- c("Cluster 1 (91 days)", "Cluster 2 (268 days)", "Cluster 3 (193 days)", "Cluster 4 (178 days)","Cluster 5 (1 day; day 368)")

xtable(confi)


# heatmap
prob <- modelsim$sims.list$p

cm <- vector()

for (i in 1:10){
    cm <- append(cm,colMeans(prob[,,i]))
}

dayname <-  rep(1:731,10)
clname <- c(rep(1,731), rep(2,731), rep(3,731), rep(4,731), rep(5,731), rep(6,731), rep(7,731), rep(8,731), rep(9,731), rep(10,731))

prob2 <- as.data.frame(cbind(clname,dayname,cm))

colnames(prob2) <- c("Cluster", "Day", "Probability")

library(ggplot2)
pdf("heatmap.pdf")
ggplot(prob2, aes(x=Cluster, y=Day, z=Probability)) + geom_tile(aes(fill = Probability)) +
theme_bw() +
scale_fill_gradient(low="white", high="black")
dev.off()

trial <- prob[,1,]
sum(colMeans(trial))



########## monthly and seasonal distributions
# generate an array of dates  1/1/2011 - 31/12/2012
dates <- seq(as.Date("2011/1/1"), as.Date("2012/12/31"), "days")
d1 <- dates[which(cluster==1)]
d2 <- dates[which(cluster==2)]
d3 <- dates[which(cluster==3)]
d4 <-dates[which(cluster==4)]
d5 <- dates[which(cluster==5)]

vie <- function(x){
    lk <- data.frame(date = x,
    year = as.numeric(format(x, format = "%Y")),
    month = as.numeric(format(x, format = "%m")),
    day = as.numeric(format(x, format = "%d")))
    
    t(as.data.frame(table(lk$month)))
}


vie(d1)
vie(d2)
vie(d3)
vie(d4)
vie(d5)



round(as.numeric(vie(d1)[2,])/sum(as.numeric(vie(d1)[2,])),2)

round(as.numeric(vie(d2)[2,])/sum(as.numeric(vie(d2)[2,])),2)

round(as.numeric(vie(d3)[2,])/sum(as.numeric(vie(d3)[2,])),2)

round(as.numeric(vie(d4)[2,])/sum(as.numeric(vie(d4)[2,])),3)

par(mfrow=c(2,2))
pie(c(23,32,44),labels = c("Winter 23%","Spring 32%","Autumn 44%"),main ="Cluster 1" )
pie(c(23,24,30,24), labels = c("Winter 23%","Spring 24%","Summer 30%", "Autumn 24%"),main ="Cluster 2")
pie(c(25,30,21,24),labels = c("Winter 25%","Spring 30%","Summer 21%", "Autumn 24%"),main ="Cluster 3")
pie(c(28,19,37,17),labels = c("Winter 28%","Spring 19%","Summer 36%", "Autumn 17%"),main ="Cluster 4")

par(mfrow=c(1,1))

library(xtable)
xtable(vie(d1))
xtable(vie(d2))
xtable(vie(d3))
xtable(vie(d4))
xtable(vie(d5))

#########  pie chart for each cluster

#p1 <- as.data.frame(cbind(namecol, colMeans(o1)))
#colnames(p1) <- c("Particle","value")
#library(ggplot2)
#pieplot <- function(h){
#  pp <- as.data.frame(cbind(namecol, colMeans(h)))
#  colnames(pp) <- c("Particle","value")
#  bp<- ggplot(pp, aes(x="", y=value, fill=Particle))+
#    geom_bar(width = 1, stat = "identity")
#  piechart <- bp + coord_polar("y", start=0)
#  piechart + scale_fill_grey() + theme_minimal()
#}


