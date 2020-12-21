## ----unnamed-chunk-1---- ##
#install.packages("reshape2")
#install.packages("unmarked")
#install.packages("fields")
 
 
 
## ----unnamed-chunk-2---- ##
lambda<- 5
 
 
 
## ----unnamed-chunk-3---- ##
kernalsPerBag<- rpois(6,lambda)
 
 
 
## ----unnamed-chunk-4---- ##
ourData<- matrix(c(0,1,0,1,2,0,
    3,2,3,2,2,0,
    2,3,4,4,4,1,
    2,3,1,1,4,1,
    1,4,3,1,5,1,
    1,2,5,1,5,0,
    2,2,3,3,6,0,
    2,4,1,4,5,0,
    0,1,5,0,5,0,
    2,1,1,1,2,0,
    2,2,1,2,6,0,
    2,2,1,2,2,1,
    3,3,3,2,6,1,
    0,1,0,2,2,0,
    2,1,2,0,2,0),ncol=6,nrow=15,byrow=TRUE)
ourData
 
 
 
## ----unnamed-chunk-5---- ##
# TRANSPOSE THE DATA TO HAVE 'VISITS' AS COLUMNS   
# AND SITES AS ROWS 
ourData<-t(ourData) 
head(ourData)
 
 
 
## ----unnamed-chunk-6---- ##
maxCounts<- apply(ourData,1, max)
maxCounts
 
 
 
## ----unnamed-chunk-7---- ##
trueValues<- c(4, 5, 6, 4, 11, 1)
plot(x=trueValues,y=maxCounts)
abline(0,1) # add a 1:1 line
 
 
 
## ----unnamed-chunk-8---- ##
library(unmarked)
data <- unmarkedFramePCount(y = ourData)

# ~DETECTION ~ ABUNDANCE
fit <- pcount(~1 ~ 1, # P THEN LAMBDA
    data=data, 
    K=50) # SET THIS HIGHER THAN YOUR EXPECTED ABUNDANCE
summary(fit)
plogis(coef(fit)[2])

# estimates of N
N_hat<- bup(ranef(fit)) # s4 class
plot(trueValues,N_hat, xlab="True density", ylab="Predicted density")
abline(0,1)# a 1:1 line
 
 
 
## ----message=FALSE---- ##
# Prepare data
library(unmarked)
data <- unmarkedFramePCount(y = y)

# ~DETECTION ~ ABUNDANCE
fit <- pcount(~1 ~ 1, # P THEN LAMBDA
    data=data, 
    K=50) # SET THIS HIGHER THAN YOUR EXPECTED ABUNDANCE
summary(fit)
 
 
 
## ----unnamed-chunk-9---- ##
nsamples<- 50 # i = 1,2,3,...20
beta_0<- 1.386 # UNDERLYING DENSITY
gamma_0<- -0.405 # LOG ODDS CAPTURE PROBABILITY

# TRANSFORM TO REAL VALUES
lambda <- exp(beta_0)# close to 4
lambda
p<- exp(gamma_0)/(1+exp(gamma_0) )
p # close to 0.4
 
 
 load("study-area.Rdata")
## ----unnamed-chunk-10---- ##
# SIMULATE ABUNDANCES 
set.seed(1985)# FOR REPRODUCABILITY; LAST YEAR DLR WAS IN VAN HALEN
sa$N<- rpois(nrow(sa),lambda)
 
 
 
## ----unnamed-chunk-11---- ##
sample_indx<- sample(1:nrow(sa),
    nsamples,replace=FALSE)
sampleSites<- sa[sample_indx,]
 
 
 
## ----unnamed-chunk-12---- ##
# GENERATE CAPTURE HISTORIES
visits<-5 # k = 1,2,3,4,5
# MATRIX TO HOLD VALUES
y<- matrix(0,nsamples,visits) # ROW FOR EAC SAMPLE SITE
for(i in 1:nsamples) # LOOP OVER EACH SAMPLE SITE
	{
    for(k in 1:visits)# LOOP OVER EACH VISIT AT EACH SITE
        {
        y[i,k]<- rbinom(1,sampleSites$N[i],p)#obs count for visit k and site i
        }
	}
 
 
 
## ----unnamed-chunk-13---- ##
head(y)
 data <- unmarkedFramePCount(y = y)

# ~DETECTION ~ ABUNDANCE
fit <- pcount(~1 ~ 1, # P THEN LAMBDA
    data=y, 
    K=50) # SET THIS HIGHER THAN YOUR EXPECTED ABUNDANCE
summary(fit)
 
 
## ----unnamed-chunk-14---- ##
# Density
lambda
# ESTIMATE IS ON LOG SCALE
exp(coef(fit)[1]) # should be close to lambda

# Capture probability
p
# ESTIMATE IS ON LOG ODDS SCALE
exp(coef(fit)[2])/(1+exp(coef(fit)[2])) # should be close p
 
 
 
## ----unnamed-chunk-15---- ##
plot(N~depth,sa,ylab="Abundance",xlab="Depth",las=1)
 
 
 
## ----unnamed-chunk-16---- ##
nsamples<- 40
indx<- sample(1:nrow(sa),nsamples)
sampleSites<- sa[indx,]
 
 
 
## ----unnamed-chunk-17---- ##
# GENERATE CAPTURE HISTORIES
visits<-5
p<- exp(gamma_0+gamma_1*sampleSites$depth)/
    (1+exp(gamma_0+gamma_1*sampleSites$depth))
y<- matrix(0,nsamples,visits)
for(i in 1:nsamples)
	{
	y[i,]<- rbinom(visits,sampleSites$N[i],p[i])
	}
 
 
 
## ----unnamed-chunk-18---- ##
# PREPARE DATA
data <- unmarkedFramePCount(y = y,
    siteCovs=data.frame(depth=sampleSites$depth))
# FIT THE MODEL WITH DEPTH AS A COVARIATE FOR LAMBDA AND P
fit <- pcount(~depth +1 ~depth+1, 
    data=data, 
    K=150)
fit
 
 
 
## ----unnamed-chunk-19---- ##
n_reps<- 10000
N_sim<- matrix(0,nrow=nrow(sa),ncol=n_reps)

for(i in 1:n_reps)
    {
    N_sim[,i]<- rpois(nrow(sa),
        lambda=sa$pred)
    
    }
    totalN<- colSums(N_sim)
    hist(totalN) 
abline(v=sum(sa$N))   
 
 
 
## ----unnamed-chunk-20---- ##
brks<-c(0,2000,2300,10000)#breakpoints
labs<-c("Small (<2000)","Medium (2000-2300)",
    "Large (2300+)")

totalN_b<-cut(x=totalN,
    breaks=brks,
    labels=labs,
    inlude.lowest=TRUE)
table(totalN_b)
table(totalN_b)/n_reps
 
 
 
