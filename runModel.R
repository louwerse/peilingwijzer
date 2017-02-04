# Peilingwijzer replication script ----------------------------------------
# Author: Tom Louwerse, Leiden University


# Load packages -----------------------------------------------------------
# Script requires JAGS (4.2) and the following packages to be installed

# Read libraries
library(rjags)
library(coda)
library(doSNOW)
library(foreach)


# Load sourcedata ---------------------------------------------------------
load("sourcedata.RData")


# Run model ---------------------------------------------------------------

# Function for running model with multiple house effects
runmodel_mhe = function(dat, 
                        margin=rep(0.0005,nrow(dat)), 
                        outcome=0,
                        n.adapt=1e5,
                        n.iter=5e5,
                        thin=1e2,
                        electionDate=as.Date("2012-09-12"),
                        bugfile="Modellen/current_model.bug", 
                        start=NA) {

  # Recode houseperiod, so that minimum is always 1
  dat$houseperiod <- as.integer(dat$houseperiod) - 
    min(as.integer(dat$houseperiod)) + 1
  
  # Start transition
  startTransition = 2
  if(!is.na(start)) {
    startTransition <- as.numeric(julian(start, origin=as.Date("2010-01-01")) - 
                                    julian(electionDate, origin=as.Date("2010-01-01")) + 1)
  }
  
  # Prepare data for JAGS
  foo <- list(y=dat$percentage,
              vari=dat$var,
              date=as.numeric(dat$datenum - 
                                julian(electionDate, 
                                       origin=as.Date("2010-01-01")) + 1),
              houseperiod=dat$houseperiod,
              org=as.integer(as.factor(as.character(dat$bedrijf))),
              margin=margin,
              NHOUSE=length(unique(dat$bedrijf)),
              NPOLLS=length(dat$percentage),
              NPERIODS=length(julian(electionDate, 
                                     origin=
                                       as.Date("2010-01-01")):max(dat$datenum)),
              NHOUSEPERIODS=length(unique(dat$houseperiod)),
              STARTTRANSITION=startTransition,
              alpha=c(outcome, 
                      rep(NA,length(julian(electionDate, 
                                           origin=as.Date("2010-01-01")):max(dat$datenum))-1))
  )
  
  # Define weights for pollsters
  house_names <- sort(unique(dat$bedrijf))
  house_weights <- matrix(1, foo$NHOUSE, foo$NHOUSEPERIODS)
  
  # Exclude I&O from house effect average for houseperiod 1, because it did not poll
  if("I&O Research" %in% house_names & foo$NHOUSEPERIODS >= 5) {
    io_which <- which(house_names == "I&O Research")
    house_weights[io_which,1] <- 0 # set weight for I&O pre-2014 to 0
  }
  
  # Exclude LISS from house effect average for houseperiods 1-4, 
  # because it did not poll
  if("LISS" %in% house_names) {
    liss_which <- which(house_names == "LISS")
    
    # set weights for periods 1-4
    # This is somewhat complicated, because smaller parties start in 
    # houseperiod 2, 3, or 4 (but houseperiods are always numbered from 1)
    if(foo$NHOUSEPERIODS == 5) house_weights[liss_which,1] <- 0 
    if(foo$NHOUSEPERIODS >= 4) house_weights[liss_which,foo$NHOUSEPERIODS-3] <- 0 
    if(foo$NHOUSEPERIODS >= 3) house_weights[liss_which,foo$NHOUSEPERIODS-2] <- 0 
    if(foo$NHOUSEPERIODS >= 2) house_weights[liss_which,foo$NHOUSEPERIODS-1] <- 0 
  }
  
  # Create matrix for house weights
  house_weights <- house_weights / matrix(colSums(house_weights), 
                                          foo$NHOUSE, 
                                          foo$NHOUSEPERIODS, 
                                          byrow=TRUE)
  foo$house_weights <- house_weights
  
  # Parties that were not polled before a certain date, are set to 0 beforehand
  if(!is.na(start)) {
    setToZero <- julian(electionDate, origin=as.Date("2010-01-01")):
      julian(start, origin=as.Date("2010-01-01")) -
      julian(electionDate, origin=as.Date("2010-01-01")) +1
    foo$alpha[setToZero]  <- 0
  }
  
  inits <- list(y2=dat$percentage)
  
  # Initialize the model
  model = rjags::jags.model(bugfile, foo, n.adapt=n.adapt, inits=inits)
  
  # Run the model
  est = rjags::coda.samples(model, c("alpha", "sigma", "house","pie"), n.iter=n.iter, thin=thin)
  return(est)
}



# Run model for each party ------------------------------------------------
cat("Start analysis\n")

# Initialize parallel processing
workerCount=min(7, as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))-1) 
w <- makeCluster(workerCount)
registerDoSNOW(w)

outputs <- foreach(i=1:length(sourcedata)) %dopar%
  runmodel_mhe(dat=sourcedata[[i]]$dat, 
               sourcedata[[i]]$dat$margin,
               outcome=sourcedata[[i]]$outcome, 
               n.adapt=1e4,n.iter=2e5,thin=80, 
               electionDate=as.Date("2012-09-12"),
               start=sourcedata[[i]]$start,
               bugfile=sourcedata[[i]]$model)
save(outputs, file="outputs.RData")

stopCluster(w)



# Basic results -----------------------------------------------------------
# This provides point estimates and 95% credible intervals
# for percentage party support on most recent day

partylist <- sapply(sourcedata, function(q) q$partyname)

getData <- function(X) {
  alpha.pos <- grep("^alpha",colnames(X[[1]]))
  alpha.last.pos <- grep(paste("^alpha\\[", max(alpha.pos),"\\]", sep=""),
                         colnames(X[[1]]))
  x <- X[[1]][,alpha.last.pos]
  return(c(mean(x),quantile(x,c(0.025,.975))))
}

alldata <- foreach(i=1:length(outputs), .combine=rbind) %do% {
  getData(outputs[[i]])
}

rownames(alldata) <- partylist
colnames(alldata) <- c("Percentage", "Margin_Low", "Margin_High")

print(round(alldata * 100, 1))


