#' Generate continuous usable data set from a
#'
#' @param treatment A column containing the assigned treatments so far
#' @param values A column containing the values from the treatments
#' @param data A data set containing the above columns
#' @param burnin The number of burn-in values desired for each treatment group
#' @param seed A set seed to help with reproducibility
#' @param nu The clinically meaningful significant difference to be found
#' @param trtlvls A vector containing the labels of each respective group
#' @param weight The weight of Sample Size Re-Estimation and Random Adaptive Re-Allocation
#' @param alpha The Type 1 error rate
#' @param power The desired power of the study
#'
#' @return A data set containing the previous data set with the attached next allocations.
#' @export
#'
#' @examples IMAD_cont_dataset()
IMAD_cont_dataset <- function(treatment, values, data, burnin = 10, seed, nu, trtlvls = c("A", "B"), weight = 0.5, alpha = 0.05, power = 0.8){
  #Set warning for lack of seed and set given seed.
  if(is.na(seed)){
    warning("Please set your seed. Keep track of this seed for all future randomization schemes.")
  } else {
    set.seed(seed)
  }

  #Get initial value
  if(is.na(treatment)){
    #Get initial burnin vector;
    burnininit <- as.factor(sample(c(rep(trtlvls[1], burnin), rep(trtlvls[2], burnin))))
    UpdatedTrt <- burnininit[1]
    newdf <- as.data.frame(UpdatedTrt)
    newdf$Values <- rep(NA, 1)
    return(newdf)
  }

  #Set type 2 Error
  beta <- 1 - power
  #Set the treatment vector
  trtvec <- data[, treatment]
  #Identify the factors of the treatment vector
  trtfac <- as.factor(trtvec)
  #Identify if any treatment levels are missing
  for(i in 1:length(trtvec)){
    if(is.na(trtvec[i])){
      stop("At least one treatment group has a missing assigment")
    }
  }

  #Identify if there are not yet 2 treatment levels
  if(length(trtlvls < 2)){
    trtlvls = trtlvls
  } else {
    #Get the treatment names of the treatment vector
    trtlvls <- levels(trtfac)
  }

  #Identify if there are more than 2 treatment levels
  if(length(trtlvls) > 2){
    warning("There are more than 2 treatment levels. Make sure you only have two treatment levels.")
  }

  #Reset seed to get same burninin
  #set.seed(seed)
  burnininit <- as.factor(sample(c(rep(trtlvls[1], burnin), rep(trtlvls[2], burnin))))

  #Isolate the two vectors affiliated with each group
  trtvec1 <- data[trtfac == trtlvls[1], values]
  trtvec2 <- data[trtfac == trtlvls[2], values]

  #Get randomization if burnin not yet met
  if(length(trtvec) < 2 * burnin){

    neededburnin <- burnininit[1:length(trtvec)]
    burnina <- ifelse(neededburnin == trtlvls[1], 1, 0)
    burninb <- ifelse(neededburnin == trtlvls[2], 1, 0)
    veca <- ifelse(trtfac == trtlvls[1], 1, 0)
    vecb <- ifelse(trtfac == trtlvls[2], 1, 0)

    if(sum(burnina == veca) == length(neededburnin) & sum(burninb == vecb) == length(neededburnin)){
      UpdatedTrt <- c(trtfac, burnininit[length(trtvec) + 1])
      newdf <- as.data.frame(UpdatedTrt)
      newdf$Values <- c(data[, values], NA)
    }

    else{
      warning("Check seed: Set seed gives allocation scheme different than previously allotted.
                We recommend using the same seed for reproducibility.
                An allocation will be given based on current dataset without considering seed.")

      #Set number of values to randomize
      grp1lft <- burnin - length(trtvec1)
      grp2lft <- burnin - length(trtvec2)
      #Get vectors of each group needed for each group
      grp1vals <- rep(trtlvls[1], grp1lft)
      grp2vals <- rep(trtlvls[2], grp2lft)
      #Generate randomization of remaining values to be randomized
      randsmp <- as.factor(sample(c(grp1vals, grp2vals)))
      #Isolate first random sample
      randsmp1 <- randsmp[1]
      #Update dataset with newest treatment group
      UpdatedTrt <- c(trtfac, randsmp1)
      #Get updated dataset of proper length
      newdf <- as.data.frame(UpdatedTrt)
      #Bring in previous values
      newdf$Values <- c(data[, values], rep(NA, length(randsmp1)))
    }
  }

  #Get randomization if burnin has been met
  if(length(trtvec) >= (2 * burnin)){

    updatetreat <- burnininit

    if (length(trtvec) > (2 * burnin)){
      for(i in (2 * burnin):(length(trtvec) - 1)){
        df <- data[1:i,]
        trtfacrep <- as.factor(df[, treatment])
        newtrtvec1 <- df[trtfacrep == trtlvls[1], values]
        newtrtvec2 <- df[trtfacrep == trtlvls[2], values]
        nextval <- IMAD_next_val_cont(newtrtvec1, newtrtvec2, trtlvls, nu)
        allogroup <- as.factor(nextval$allocation_group)
        updatetreat <- c(updatetreat, allogroup)
      }
    }

    totveca <- ifelse(updatetreat == trtlvls[1], 1, 0)
    totvecb <- ifelse(updatetreat == trtlvls[2], 1, 0)

    veca <- ifelse(trtfac == trtlvls[1], 1, 0)
    vecb <- ifelse(trtfac == trtlvls[2], 1, 0)

    if(sum(veca == totveca) == length(trtvec) & sum(vecb == totvecb) == length(trtvec)){
      #Get next value given current vectors given burnin has been met
      nextval <- IMAD_next_val_cont(trtvec1, trtvec2, trtlvls, nu)
      #Assign the allocation group
      allogroup <- as.factor(nextval$allocation_group)
      #Update dataset with newest treatment group
      UpdatedTrt <- c(trtfac, allogroup)
      #Get updated dataset of proper length
      newdf <- as.data.frame(UpdatedTrt)
      #Bring in previous values
      newdf$Values <- c(data[, values], rep(NA, length(allogroup)))
    }


    else{
      warning("Check seed: Set seed gives allocation scheme different than previously allotted.
                We recommend using the same seed for reproducibility.
                An allocation will be given based on current dataset without considering seed.")

      #Get next value given current vectors given burnin has been met
      nextval <- IMAD_next_val_cont(trtvec1, trtvec2, trtlvls, nu)
      #Assign the allocation group
      allogroup <- as.factor(nextval$allocation_group)
      #Update dataset with newest treatment group
      UpdatedTrt <- c(trtfac, allogroup)
      #Get updated dataset of proper length
      newdf <- as.data.frame(UpdatedTrt)
      #Bring in previous values
      newdf$Values <- c(data[, values], rep(NA, length(allogroup)))
    }

  }
  return(newdf)
}
