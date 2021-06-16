# Set seed
set.seed(112358)

########################
## Necessary packages ##
########################

# Load packages
library(foreign)
library(pROC)
library(igraph)

###############
## Functions ##
###############

# New AUC function
au2 <- function(response,
                predictor,
                direction){
  test <- sum(response == 0,
              na.rm = TRUE) > 0
  test <- test & sum(response == 1,
                     na.rm = TRUE) > 0
  if(test){
    return(auc(response = response,
               predictor = predictor,
               direction = direction))
  } else {
    return(NA)
  }
}

# Create one HIV variable at each wave
oneHIV <- function(data){
  names(data)[names(data) == "labinfected"] <- "hiv"
  hivNA <- is.na(data$hiv)
  data$hiv[hivNA] <- data$selfdiag[hivNA]
  data$selfdiag <- NULL
  data$hiv <- as.numeric(data$hiv == "YES")
  return(data)
}

# Simulate one round
sim <- function(pLambda){
  
  # Number of people
  N <- nrow(pLambda)
  
  # Number of venues in population
  M <- ncol(pLambda)
  
  # Simulate Poisson processes
  Z <- matrix(data = rpois(n = N * M,
                           lambda = as.numeric(pLambda)),
              nrow = N,
              ncol = M)
  nHits <- rowSums(Z)
  K <- sum(nHits)
  
  # Data
  # One column for each person-encounter
  # First column is subject number
  # Second column is time of encounter
  # Third column is venue
  data <- matrix(data = as.numeric(NA),
                 nrow = K,
                 ncol = 3)
  data[ ,1] <- rep(x = 1:N,
                   times = nHits)
  data[ ,2] <- runif(n = K)
  data[ ,3] <- rep(x = rep(x = 1:M,
                           times = N),
                   times = as.numeric(t(Z)))
  data <- data[order(data[ ,2]), ]
  
  # Wide data
  # List
  # One element of list for each venue UNLESS there are 0 or 1 encounters
  # at that venue
  # Note that this means the number of elements in "wide" might be less
  # than M
  wide <- list()
  for(j in 1:M){
    temp <- data[ ,3] == j
    if(sum(temp) %in% c(0,1)){
      wide[[j]] <- NULL
    } else {
      wide[[j]] <- data[temp, ]
      count <- nrow(wide[[j]])
      if(count %% 2 != 0){
        wide[[j]] <- wide[[j]][1:(count - 1), ]
      }
    }
  }
  
  # Return
  return(wide)
}

# Create sexual contact matrix
# Number of sexual encounters between all individuals.
# Note that this means a couple that has sex multiple times gets counted
# all those times.
# In the form of an adjacency matrix (with zero diagonal).
sexConMat <- function(wide,N){
  A2 <- matrix(data = 0,
               nrow = N,
               ncol = N)
  m <- length(wide)
  for(j in 1:m){
    if(!is.null(wide[[j]])){
      n2 <- nrow(wide[[j]]) / 2
      for(i in 1:n2){
        ego <- wide[[j]][2 * i - 1,1]
        alter <- wide[[j]][2 * i,1]
        A2[ego,alter] <- A2[ego,alter] + 1
        A2[alter,ego] <- A2[alter,ego] + 1
      }
    }
  }
  diag(A2) <- 0
  return(A2)
}

# Determine infections
detInf <- function(A2,y,pTrans){
  hiv <- colSums(A2 * y)
  hiv <- rbinom(n = length(hiv),
                size = 1,
                prob = 1 - (1 - pTrans)^hiv)
  hiv[y] <- 1
  return(hiv)
}

# Create affiliation matrix, but delete occurrences when someone
# matched to himself
# Returns matrix in decreasing order of column sums
getAffMat <- function(wide1,N){
  m1 <- length(wide1)
  A1 <- matrix(data = 0,
               nrow = N,
               ncol = m1)
  for(j in 1:m1){
    for(i in 1:N){
      present <- wide1[[j]][ ,1] == i
      K <- length(present)
      odd <- rep(x = c(TRUE,FALSE),
                 times = K / 2)
      even <- rep(x = c(FALSE,TRUE),
                  times = K / 2)
      double <- present[odd] & present[even]
      A1[i,j] <- sum(present) - 2 * sum(double)
    }
  }
  A1 <- A1[ ,order(colSums(A1),
                   decreasing = TRUE)]
  
  return(A1)
}

# Function to calculate Q based on adjacency matrix and HIV statuses
getQ <- function(A,y){
  q <- t(A) %*% y / colSums(A)
  return(q)
}

# Function to get risk
getR <- function(q,z,pTrans){
  terms <- ((1 - pTrans) * q + 1 - q)^z
  r <- 1 - prod(terms)
  return(r)
}

# Function to calculate AUC using new estimator
getNewAUC <- function(A,y,q,pTrans){
  # Calculate risk for people who were HIV-negative at Wave 1
  # and have HIV data at Wave 2
  n <- length(y)
  risk <- rep(x = NA,
              times = n)
  for(i in 1:n){
    z <- A[i, ]
    risk[i] <- getR(q = q,
                    z = z,
                    pTrans = pTrans)
  }
  
  # Test performance of risk estimator
  return(au2(response = y,
             predictor = risk,
             direction = "<"))
}

# Get AUCs
getAllAUC <- function(a1,
                      hivMidSam,
                      negMidSam,
                      hivNewSam,
                      pTrans){
  
  # Estimate prevalence per venue, based only on sample
  q <- getQ(A = a1,
            y = hivMidSam)
  
  # Get AUC from new method
  newAUC <- getNewAUC(A = a1[negMidSam, ],
                      y = hivNewSam[negMidSam],
                      q = q,
                      pTrans = pTrans)
  
  # Determine if logistic regression does as well
  mod4 <- glm(formula = hivMidSam ~ .,
              family = binomial(),
              data = data.frame(cbind(a1,hivMidSam)))
  logAUC <- au2(response = hivNewSam[negMidSam],
                predictor = mod4$fitted.values[negMidSam],
                direction = "<")
  
  # Determine if logistic regression does as well, future knowledge
  mod5 <- glm(formula = hivNewSam ~ .,
              family = binomial(),
              data = data.frame(cbind(a1,hivNewSam)))
  futAUC <- au2(response = hivNewSam[negMidSam],
                predictor = mod5$fitted.values[negMidSam],
                direction = "<")
  
  # Determine if logistic regression does as well
  mod6 <- glm(formula = hivMidSam ~ colSums(t(a1)),
              family = binomial)
  logAUC6 <- au2(response = hivNewSam[negMidSam],
                 predictor = mod6$fitted.values[negMidSam],
                 direction = "<")
  
  # Determine if logistic regression does as well, future knowledge
  mod7 <- glm(formula = hivNewSam ~ colSums(t(a1)),
              family = binomial)
  futAUC7 <- au2(response = hivNewSam[negMidSam],
                 predictor = mod7$fitted.values[negMidSam],
                 direction = "<")
  
  return(c(logAUC,
           logAUC6,
           futAUC7,
           newAUC,
           futAUC))
}

# Big simulation
bigSim <- function(pLambda,
                   y,
                   N,
                   contact = FALSE,
                   pTrans = 1,
                   combine = 3,
                   smallest = 3,
                   biggest = 3,
                   contaminate = 0.5,
                   n = -1,
                   qBasedOnX = FALSE){
  
  # Set sample size
  if(n < 0){
    n <- nrow(pLambda)
  }
  
  # Number of venues in population
  M <- ncol(pLambda)
  
  # Create population
  kPop <- sample(x = nrow(pLambda),
                 size = N,
                 replace = TRUE)
  pLambda <- pLambda[kPop, ]
  yPop <- y[kPop]
  
  # Simulate one round
  wide1 <- sim(pLambda = pLambda)
  
  # Create sexual contact matrix
  A2 <- sexConMat(wide = wide1,
                  N = N)
  
  # Determine infections
  hivMid <- detInf(A2 = A2,
                   y = yPop,
                   pTrans = pTrans)
  
  # Create affiliation matrix
  A1 <- getAffMat(wide1 = wide1,
                  N = N)
  
  # Sample
  kSam <- sort(sample(x = N,
                      size = n))
  
  # Simulate another round
  wide2 <- sim(pLambda = pLambda)
  
  # Base q on X instead of Z?
  if(qBasedOnX){
    A1 <- getAffMat(wide1 = wide2,
                    N = N)
  }
  
  # Sample
  hivMidSam <- hivMid[kSam]
  a1 <- A1[kSam, ]
  
  # Coarse
  d <- ncol(a1) %/% combine
  indicator <- ncol(a1) %% combine > 0
  a1Combine <- matrix(data = 0,
                      nrow = nrow(a1),
                      ncol = d + indicator)
  for(i in 1:d){
    a1Combine[ ,i] <- rowSums(a1[ ,((i - 1) * combine + 1):(i * combine)])
  }
  if(indicator){
    columns <- (d * combine + 1):ncol(a1)
    if(length(columns) > 1){
      a1Combine[ ,(d + 1)] <- rowSums(a1[ ,columns])
    } else {
      a1Combine[ ,(d + 1)] <- a1[ ,columns]
    }
  }
  rm(d,
     indicator)
  
  # Smallest
  a1Smallest <- a1[ ,1:(ncol(a1) - smallest)]
  
  # Biggest
  a1Biggest <- a1[ ,(biggest + 1):ncol(a1)]
  
  # Contaminate
  a1Contam <- a1
  for(i in 1:nrow(a1Contam)){
    total <- sum(a1Contam[i, ])
    if(total > 0){
      newVec <- round(a1Contam[i, ] * (1 - contaminate))
      a1Contam[i, ] <- newVec + rmultinom(n = 1,
                                          size = total - sum(newVec),
                                          prob = rep(x = 1,
                                                     times = length(newVec)))
    }
  }
  rm(total,
     newVec)
  
  # Calculate quartiles for number of partners
  A2TF <- A2 >= 1
  numPar <- colSums(A2TF)
  quar <- unname(quantile(numPar[kSam]))[2:4]
  
  # Clean up
  rm(A2TF,
     numPar)
  
  # Create sexual contact matrix
  A2 <- sexConMat(wide = wide2,
                  N = N)
  
  # Determine infections
  hivNew <- detInf(A2 = A2,
                   y = hivMid,
                   pTrans = pTrans)
  
  # Who was HIV-negative in sample?
  hivNewSam <- hivNew[kSam]
  negMidSam <- (hivMidSam == 0)
  
  # Determine number of new infections between Wave 1 and Wave 2
  numNewInf <- sum(hivNewSam[negMidSam])
  
  # Determine number of person-days at-risk between Wave 1 and Wave 2
  pdar <- (sum(negMidSam) - numNewInf) * 6 * 31 + numNewInf * 3 * 31
  
  # Determine number of new infections per 100 person-years at-risk
  rate <- numNewInf / pdar * 365.25
  
  # Vector to return
  vec <- c(getAllAUC(a1 = a1,
                     hivMidSam = hivMidSam,
                     negMidSam = negMidSam,
                     hivNewSam = hivNewSam,
                     pTrans = pTrans),
           getAllAUC(a1 = a1Combine,
                     hivMidSam = hivMidSam,
                     negMidSam = negMidSam,
                     hivNewSam = hivNewSam,
                     pTrans = pTrans),
           getAllAUC(a1 = a1Smallest,
                     hivMidSam = hivMidSam,
                     negMidSam = negMidSam,
                     hivNewSam = hivNewSam,
                     pTrans = pTrans),
           getAllAUC(a1 = a1Biggest,
                     hivMidSam = hivMidSam,
                     negMidSam = negMidSam,
                     hivNewSam = hivNewSam,
                     pTrans = pTrans),
           getAllAUC(a1 = a1Contam,
                     hivMidSam = hivMidSam,
                     negMidSam = negMidSam,
                     hivNewSam = hivNewSam,
                     pTrans = pTrans),
           quar,
           rate)
  names(vec) <- c("Logistic 1 AUC",
                  "Logistic 1 AUC Total",
                  "Logistic 2 AUC Total",
                  "New Method AUC",
                  "Logistic 2 AUC",
                  "Logistic 1 AUC",
                  "Logistic 1 AUC Total",
                  "Logistic 2 AUC Total",
                  "New Method AUC",
                  "Logistic 2 AUC",
                  "Logistic 1 AUC",
                  "Logistic 1 AUC Total",
                  "Logistic 2 AUC Total",
                  "New Method AUC",
                  "Logistic 2 AUC",
                  "Logistic 1 AUC",
                  "Logistic 1 AUC Total",
                  "Logistic 2 AUC Total",
                  "New Method AUC",
                  "Logistic 2 AUC",
                  "Logistic 1 AUC",
                  "Logistic 1 AUC Total",
                  "Logistic 2 AUC Total",
                  "New Method AUC",
                  "Logistic 2 AUC",
                  "Num Partners Q1",
                  "Num Partners Q2",
                  "Num Partners Q3",
                  "Infection Rate")
  
  if(contact){
    # Return sample sexual contact network
    return(A2[kSam,kSam])
  } else {
    # Return AUC
    return(vec)
  }
  
}

#####################
## Manipulate data ##
#####################

# Load data
data1 <- read.dta(file = "Wave 1/EGO.dta")
data2 <- read.dta(file = "Wave 2/EGO.dta")
data3 <- read.dta(file = "Wave 1/SEXUAL.dta")

# Keep useful variables
keep1 <- c("su_id",
           "idate",
           "agecalc",
           "seed",
           "labinfected",
           "selfdiag")
keep2 <- c("su_id",
           "idate",
           "labinfected",
           "selfdiag")
keep3 <- c("su_id",
           "numpartners",
           "lastsexdt",
           "numsex",
           "partnersex",
           "howmeet",
           "howmeetapp")
data1 <- data1[ ,keep1]
data2 <- data2[ ,keep2]
data3 <- data3[ ,keep3]

# Clean up
rm(keep1,
   keep2,
   keep3)

# Number of participants at Wave 1
nrow(data1)

# Number of sex partners at Wave 1
nrow(data3)

# Number of participants at Wave 2
nrow(data2)

# Range of dates for Wave 1
range(data1$idate)

# Range of dates for Wave 2
range(data2$idate)

# Number of seeds
sum(data1$seed)

# Range of ages
range(data1$agecalc)

# Were there any new participants at Wave 2?
sum(!(data2$su_id %in% data1$su_id))

# How many participants are missing lab HIV test at each wave?
sum(is.na(data1$labinfected))
sum(is.na(data1$selfdiag))
sum(is.na(data1$labinfected) & is.na(data1$selfdiag))
sum(is.na(data2$labinfected))
sum(is.na(data2$selfdiag))
sum(is.na(data2$labinfected) & is.na(data2$selfdiag))

# Create one HIV variable at each wave
data1 <- oneHIV(data1)
data2 <- oneHIV(data2)

# Any missing date?
sum(is.na(data1$idate))
sum(is.na(data2$idate))

# Merge the two waves
data <- merge(x = data1,
              y = data2,
              by = "su_id",
              suffixes = c("1","2"),
              all = TRUE)

# Clean up
rm(data1,
   data2)

# Delete people with missing HIV variable at Wave 1
sum(is.na(data$hiv1))
data <- data[!is.na(data$hiv1), ]
data3 <- data3[data3$su_id %in% data$su_id, ]

# Number of sex partners in the data set
nrow(data3)

# How many participants had sex partner info?
length(unique(data3$su_id))
nrow(data) - length(unique(data3$su_id))

# Add idate1 to data3
data3 <- merge(x = data3,
               y = data.frame(su_id = data$su_id,
                              idate1 = data$idate1),
               by = "su_id")

# Remove partners from greater than six months ago
data3 <- data3[as.numeric(data3$idate1 - data3$lastsexdt) <= 186, ]

# How many participants had sex partner info?
length(unique(data3$su_id))

# How many sex partners?
nrow(data3)

# Clean up
data3$idate1 <- NULL
data3$lastsexdt <- NULL

# # Remove partners who are not male or whose gender is missing
# data3 <- data3[!is.na(data3$partnersex), ]
# data3 <- data3[data3$partnersex == "Male", ]
# 
# # How many participants had sex partner info?
# length(unique(data3$su_id))
# 
# # How many sex partners?
# nrow(data3)

# Remove partnersex variable
data3$partnersex <- NULL

# Remove partners with missing howmeet
data3 <- data3[!is.na(data3$howmeet), ]

# How many participants had sex partner info?
length(unique(data3$su_id))

# How many sex partners are in the data?
nrow(data3)

# Who met at least one sex partner "through somebody else"
# or "knew each other previously"?
throughKnew <- c("THROUGH SOMEBODY ELSE YOU BOTH KNEW",
                 "KNEW EACH OTHER PREVIOUSLY (VOLUNTEERED)")
idsThroughKnew <- sort(unique(data3$su_id[data3$howmeet %in% throughKnew]))

# Remove partners with howmeet "through somebody else"
# or "knew each other previously"
data3 <- data3[data3$howmeet != "THROUGH SOMEBODY ELSE YOU BOTH KNEW", ]
data3 <- data3[data3$howmeet != "KNEW EACH OTHER PREVIOUSLY (VOLUNTEERED)", ]

# How many participants had sex partner info?
length(unique(data3$su_id))

# How many sex partners are in the data?
nrow(data3)

# Specify website, mobile app, or phone service
levels(data3$howmeet) <- c(levels(data3$howmeet),
                           levels(data3$howmeetapp))
isOnline <- data3$howmeet == "PHONE OR INTERNET"
data3$howmeet[isOnline] <- data3$howmeetapp[isOnline]

# Again, remove partners with missing howmeet
data3 <- data3[!is.na(data3$howmeet), ]

# Delete howmeetapp variable
data3$howmeetapp <- NULL

# How many participants had sex partner info?
length(unique(data3$su_id))

# How many sex partners?
nrow(data3)

# Is anyone missing numsex?
sum(is.na(data3$numsex))

# Remove partners with missing numsex
data3 <- data3[!is.na(data3$numsex), ]

# How many participants had sex partner info?
length(unique(data3$su_id))

# How many sex partners?
nrow(data3)

# How many of these participants have HIV data at Wave 2?
sum(unique(data3$su_id) %in% data$su_id[!is.na(data$hiv2)])

# How many partners?
nrow(data3)

# Venues?
summary(data3$howmeet)
sort(table(data3$howmeet))

# Transform levels to numeric values
data3$howmeet <- as.numeric(data3$howmeet)

# How many "venues" are there?
length(unique(data3$howmeet))

# Add "numpartners" to "data"
temp <- unique(data3[ ,1:2])
data <- merge(x = data,
              y = temp,
              by = "su_id")
rm(temp)
data3$numpartners <- NULL

# Summary statistics for the number of partners
summary(data$numpartners)

# Count number of rows per id in "data3"
count <- as.data.frame(table(data3$su_id))
names(count) <- c("su_id",
                  "nrows")

# Summary statistics for the number of partners asked about
summary(count$nrows)

# Merge count with data
data <- merge(x = data,
              y = count,
              by = "su_id")
data$augment <- data$numpartners / data$nrows
sum(data$augment < 1)
data$numpartners <- NULL
data$nrows <- NULL

# Clean up
rm(count)

## Change numsex to 1 for everyone
#data3$numsex <- 1

# Combine sexual encounters through same location
data3 <- data3[order(data3$su_id,
                     data3$howmeet), ]
i <- 1
while(i < nrow(data3)){
  test1 <- c(data3$su_id[i],
             data3$howmeet[i])
  test2 <- c(data3$su_id[i + 1],
             data3$howmeet[i + 1])
  if(sum(test1 == test2) == 2){
    data3$numsex[i] <- data3$numsex[i] + data3$numsex[i + 1]
    data3 <- data3[-(i + 1), ]
  } else {
    i <- i + 1
  }
}

# Clean up
rm(i,
   test1,
   test2)

# Reshape data
data3 <- reshape(data = data3,
                 v.names = "numsex",
                 timevar = "howmeet",
                 idvar = "su_id",
                 direction = "wide")
data3[is.na(data3)] <- 0

# Merge data
data <- merge(x = data,
              y = data3,
              by = "su_id")
rm(data3)

# Augment by number of partners
data[ ,9:ncol(data)] <- data[ ,9:ncol(data)] * data$augment
data$augment <- NULL

# Sample size
sum(!is.na(data$idate1))
sum(!is.na(data$hiv1))
sum(!is.na(data$idate2))
sum(!is.na(data$hiv2))

# One participant had same interview date 2 as interview date 1
same <- data$su_id[which(data$idate1 == data$idate2)]
same
# Not clearly in one wave or the other
# Note overlap in range
data$idate1[data$su_id == same]
range(data$idate1)
range(data$idate2,
      na.rm = TRUE)
test <- data[!is.na(data$idate2), ]
test <- test[test$idate2 <= "2014-07-23", ]
nrow(test)
sum(test$hiv1 != test$hiv2,
    na.rm = TRUE)
sum(test$hiv1 == 0 & test$hiv2 == 1,
    na.rm = TRUE)
# No change in HIV status
data[data$su_id == same,c("hiv1","hiv2")]
# Delete Wave 2 data
data[data$su_id == same,c("idate2","hiv2")] <- NA
# Number of people with HIV status at Wave 2
sum(!is.na(data$hiv2))

# Clean up
rm(same,
   test)

# Who is hiv negative at baseline?
neg <- data$hiv1 == 0
sum(neg)

# How many people were HIV-negative at baseline and have HIV data at Wave 2?
hasHIV2 <- !is.na(data$hiv2)
negAndHasHIV2 <- neg & hasHIV2
sum(negAndHasHIV2)

# Elapsed time between two interview dates
data$elapsed <- as.numeric(data$idate2 - data$idate1)
#hist(data$elapsed)
summary(data$elapsed)
head(sort(data$elapsed))
tail(sort(data$elapsed))
sd(x = data$elapsed,
   na.rm = TRUE)

# Calculate A
A <- as.matrix(data[ ,grepl(pattern = "numsex",
                            x = colnames(data))])

# Degree distribution for individuals
summary(colSums(t(A)))

# Estimate prevalence per venue, based only on sample
q <- getQ(A = A,
          y = data$hiv1)

# Get AUC for new method under different probabilities of transmission
getNewAUC(A = A[negAndHasHIV2, ],
          y = data$hiv2[negAndHasHIV2],
          q = q,
          pTrans = 0.62 / 100)
getNewAUC(A = A[negAndHasHIV2, ],
          y = data$hiv2[negAndHasHIV2],
          q = q,
          pTrans = 1.1 / 100)
getNewAUC(A = A[negAndHasHIV2, ],
          y = data$hiv2[negAndHasHIV2],
          q = q,
          pTrans = 1.43 / 100)

# Logistic regression, using Wave 1 as outcome
log1 <- glm(formula = V1 ~ .,
            family = binomial,
            data = as.data.frame(cbind(data$hiv1,
                                       A)))
au2(response = data$hiv2[negAndHasHIV2],
    predictor = log1$fitted.values[negAndHasHIV2],
    direction = "<")

# Logistic regression, using Wave 2 as outcome
log2 <- glm(formula = V1 ~ .,
            family = binomial,
            data = as.data.frame(cbind(data$hiv2[hasHIV2],
                                       A[hasHIV2, ])))
au2(response = data$hiv2[negAndHasHIV2],
    predictor = log2$fitted.values[neg[hasHIV2]],
    direction = "<")

# Logistic regression, using Wave 1 as outcome
log6 <- glm(formula = data$hiv1 ~ colSums(t(A)),
            family = binomial)
au2(response = data$hiv2[negAndHasHIV2],
    predictor = log6$fitted.values[negAndHasHIV2],
    direction = "<")

# Logistic regression, using Wave 2 as outcome
log7 <- glm(formula = data$hiv2[hasHIV2] ~ colSums(t(A[hasHIV2, ])),
            family = binomial)
au2(response = data$hiv2[negAndHasHIV2],
    predictor = log7$fitted.values[neg[hasHIV2]],
    direction = "<")

# Create tempData
tempData <- data[negAndHasHIV2, ]

# How many new infections were there?
numNewInf <- sum(tempData$hiv2 - tempData$hiv1)

# How many person-days at risk were there?
pdar <- sum(tempData$elapsed[tempData$hiv2 == 0],
            1 / 2 * tempData$elapsed[tempData$hiv2 == 1])

# New infections per person-years at risk
infectionRate <- numNewInf / pdar * 365.25
infectionRate

#####################
## Subset analysis ##
#####################

# Compare AUCs for only those people who didn't report
# any sex partners met through other people or that they already knew
negHIV2NotThroughKnew <- negAndHasHIV2 & !(data$su_id %in% idsThroughKnew)
getNewAUC(A = A[negHIV2NotThroughKnew, ],
          y = data$hiv2[negHIV2NotThroughKnew],
          q = q,
          pTrans = 0.62 / 100)
getNewAUC(A = A[negHIV2NotThroughKnew, ],
          y = data$hiv2[negHIV2NotThroughKnew],
          q = q,
          pTrans = 1.43 / 100)
au2(response = data$hiv2[negHIV2NotThroughKnew],
    predictor = log1$fitted.values[negHIV2NotThroughKnew],
    direction = "<")
au2(response = data$hiv2[negHIV2NotThroughKnew],
    predictor = log2$fitted.values[negHIV2NotThroughKnew[hasHIV2]],
    direction = "<")

####################
## Alternate data ##
####################

# Venue-to-venue graph
vMat <- matrix(data = 0,
               nrow = ncol(A),
               ncol = ncol(A))
for(i in 1:ncol(A)){
  for(j in 1:ncol(A)){
    vMat[i,j] <- sum(A[ ,i] > 0 & A[ ,j] > 0)
  }
}
diag(vMat) <- 0
v <- graph_from_adjacency_matrix(adjmatrix = vMat,
                                 mode = "undirected",
                                 weighted = TRUE)
plot(x = v,
     vertex.label = NA,
     layout = layout_with_graphopt,
     edge.width = E(v)$weight)

# "Largest" venues
colSums(A)
order(colSums(A),
      decreasing = TRUE)
temp <- colSums(A)
names(temp) <- 1:length(temp)
barplot(temp)
rm(temp)

# Create matrix B
B <- cbind(rbind(A,
                 A),
           matrix(data = 0,
                  nrow = nrow(A) * 2,
                  ncol = 10))
rowsToChange <- (nrow(A) + 1):(nrow(A) * 2)
colsTo <- (ncol(A) + 1):ncol(B)
colsFrom <- order(colSums(A),
                  decreasing = TRUE)[1:10]
B[rowsToChange,colsTo] <- B[rowsToChange,colsFrom]
B[rowsToChange,colsFrom] <- 0

# Alternate venue-to-venue matrix
uMat <- matrix(data = 0,
               nrow = ncol(B),
               ncol = ncol(B))
for(i in 1:ncol(B)){
  for(j in 1:ncol(B)){
    uMat[i,j] <- sum(B[ ,i] > 0 & B[ ,j] > 0)
  }
}
diag(uMat) <- 0
u <- graph_from_adjacency_matrix(adjmatrix = uMat,
                                 mode = "undirected",
                                 weighted = TRUE)
plot(x = u,
     vertex.label = NA,
     layout = layout_with_graphopt,
     edge.width = E(u)$weight)

# HIV status
altHIV1 <- rbinom(n = length(data$hiv1),
                  size = 1,
                  prob = data$hiv1 * 0.25)
mean(data$hiv1)
mean(altHIV1)

################
## Simulation ##
################

# Set N
N <- nrow(A) * 5

# Set number of replications
nRep <- 1e3

# Run simulation
result1 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = A,
                                     y = data$hiv1,
                                     N = N,
                                     pTrans = 0.62 / 100)))

# Run simulation
result2 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = A,
                                     y = data$hiv1,
                                     N = N,
                                     pTrans = 1.43 / 100)))

# # These two values of pi are too small. They result in zero transmissions.
# # Run simulation
# result <- t(replicate(n = nRep,
#                       expr = bigSim(pLambda = A,
#                                     y = data$hiv1,
#                                     N = N,
#                                     pTrans = 0.11 / 100)))
# 
# # Results
# summary(result)
#apply(X = result,
#      MARGIN = 2,
#      FUN = sd)
# 
# # Run simulation
# result <- t(replicate(n = nRep,
#                       expr = bigSim(pLambda = A,
#                                     y = data$hiv1,
#                                     N = N,
#                                     pTrans = 0.25 / 100)))
# 
# # Results
# summary(result)
#apply(X = result,
#      MARGIN = 2,
#      FUN = sd)

# Run simulation
result9 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = A,
                                     y = data$hiv1,
                                     N = N,
                                     pTrans = 0.5 / 100)))

# Run simulation
result10 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 0.75 / 100)))

# Run simulation
result11 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 1 / 100)))

# Run simulation
result12 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 1.25 / 100)))

# Run simulation
result13 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 1.5 / 100)))

# Run simulation
result14 <- t(replicate(n = nRep,
                      expr = bigSim(pLambda = A,
                                    y = data$hiv1,
                                    N = N,
                                    pTrans = 1.75 / 100)))

# Run simulation
result15 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 2 / 100)))

# Run simulation
result16 <- t(replicate(n = nRep,
                        expr = bigSim(pLambda = A,
                                      y = data$hiv1,
                                      N = N,
                                      pTrans = 1.1 / 100)))

# Run simulation
result3 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = B,
                                     y = c(data$hiv1,
                                           altHIV1),
                                     N = nrow(B) * 5,
                                     pTrans = 0.62 / 100)))

# Run simulation
result4 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = B,
                                     y = c(data$hiv1,
                                           altHIV1),
                                     N = nrow(B) * 5,
                                     pTrans = 1.1 / 100)))

# Run simulation
result5 <- t(replicate(n = nRep,
                       expr = bigSim(pLambda = B,
                                     y = c(data$hiv1,
                                           altHIV1),
                                     N = nrow(B) * 5,
                                     pTrans = 1.43 / 100)))

#############
## Figures ##
#############

# Boxplot
theseCols <- c(1,5,2:4,
               6,10,7:9,
               11,15,12:14,
               16,20,17:19,
               21,25,22:24)
at <- c(1:5,
        7:11,
        13:17,
        19:23,
        25:29)
comp <- rep(x = c(4,9,14,19,24),
            each = 5)
pdf(file = "draft13a.pdf",
    width = 8.5,
    height = 11)
op <- par()
par(mfrow = c(3,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result1[ ,theseCols[i]],result1[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result1[ ,theseCols],
        ylim = c(0,1.02),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "A",
     cex = 2)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result16[ ,theseCols[i]],result16[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result16[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "B",
     cex = 2)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result2[ ,theseCols[i]],result2[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result2[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "C",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
par(op)
dev.off()

# Values of pi
piValues <- c(0.62,
              0.75,
              1,
              1.1,
              1.25,
              1.43,
              1.5,
              1.75,
              2) / 100

# Infection rates
rates <- cbind(result1[ ,29],
               result10[ ,29],
               result11[ ,29],
               result16[ ,29],
               result12[ ,29],
               result2[ ,29],
               result13[ ,29],
               result14[ ,29],
               result15[ ,29])
meanRates <- colMeans(rates)
seRates <- apply(X = rates,
                 MARGIN = 2,
                 FUN = sd) / sqrt(nRep)

# Line plot
pdf(file = "draft13f.pdf",
    width = 11.3,
    height = 7)
plot(x = piValues,
     y = meanRates,
     ylim = range(c(meanRates - 1.96 * seRates,meanRates + 1.96 * seRates)),
     pch = 19,
     xlab = expression(pi),
     ylab = "Mean Infection Rate +/- 1.96 SE",
     type = "o",
     cex.lab = 1.3,
     cex.axis = 1.3)
arrows(x0 = piValues,
       y0 = meanRates - 1.96 * seRates,
       x1 = piValues,
       y1 = meanRates + 1.96 * seRates,
       length = 0.05,
       angle = 90,
       code = 3)
points(x = piValues[c(1,4,6)],
       y = meanRates[c(1,4,6)],
       col = "red",
       pch = "O",
       cex = 5)
dev.off()

# Venue-to-venue graph
pdf(file = "draft13h.pdf",
    width = 11.3,
    height = 7)
op <- par()
par(mfrow = c(1,2))
plot(x = v,
     vertex.label = NA,
     layout = layout_with_graphopt,
     edge.width = E(v)$weight)
text(x = -1.2,
     y = -1.2,
     labels = "A",
     cex = 2)
plot(x = u,
     vertex.label = NA,
     layout = layout_with_graphopt,
     edge.width = E(u)$weight)
text(x = -1.2,
     y = -1.2,
     labels = "B",
     cex = 2)
par(op)
dev.off()

# Boxplot
pdf(file = "draft13j.pdf",
    width = 8.5,
    height = 11)
op <- par()
par(mfrow = c(3,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result3[ ,theseCols[i]],result3[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result3[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "A",
     cex = 2)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result4[ ,theseCols[i]],result4[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result4[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "B",
     cex = 2)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result5[ ,theseCols[i]],result5[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result5[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "C",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
par(op)
dev.off()

# Boxplot (separated)
theseCols <- c(1,5,2:4,
               6,10,7:9,
               11,15,12:14,
               16,20,17:19,
               21,25,22:24)
at <- c(1:5,
        7:11,
        13:17,
        19:23,
        25:29)
comp <- rep(x = c(4,9,14,19,24),
            each = 5)
pdf(file = "draft13aa.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result1[ ,theseCols[i]],result1[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result1[ ,theseCols],
        ylim = c(0,1.02),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "A",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()
pdf(file = "draft13ab.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result16[ ,theseCols[i]],result16[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result16[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "B",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()
pdf(file = "draft13ac.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result2[ ,theseCols[i]],result2[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result2[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "C",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()

# Boxplot (separated)
pdf(file = "draft13ja.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result3[ ,theseCols[i]],result3[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result3[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "A",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()
pdf(file = "draft13jb.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result4[ ,theseCols[i]],result4[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result4[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "B",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()
pdf(file = "draft13jc.pdf",
    width = 11,
    height = 6.8)
name <- rep(x = 1:5,
            times = 5)
for(i in 1:length(at)){
  if(t.test(result5[ ,theseCols[i]],result5[ ,comp[i]])$p.value < 0.05){
    name[i] <- paste0(name[i],
                      "*")
  }
}
boxplot(result5[ ,theseCols],
        ylim = c(0,1),
        at = at,
        names = name,
        col = rep(x = c("gray20",
                        "gray40",
                        "gray60",
                        "gray80",
                        "gray100"),
                  each = 5))
abline(h = seq(from = 0,
               to = 1,
               by = 0.1),
       lty = "dotted",
       col = "lightgray")
text(x = 1,
     y = 0.1,
     labels = "C",
     cex = 2)
legend(x = "bottom",
       title = "Scenario",
       fill = c("gray20",
                "gray40",
                "gray60",
                "gray80",
                "gray100"),
       legend = c("Perfect",
                  "Coarse",
                  "Smallest",
                  "Largest",
                  "Contamination"),
       horiz = TRUE,
       bg = "white")
title(xlab = "Method",
      ylab = "AUC",
      outer = TRUE,
      line = 3)
dev.off()

save.image("draft13.RData")


