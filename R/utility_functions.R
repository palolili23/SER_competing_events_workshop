nonParametricCumHaz <-
  function(weightVector,
           inputdata,
           grp,
           outcomeProstate = TRUE) {
    outputHazards <- rep(NA, length.out = length(cutTimes))
    counter <- 1
    for (i in cutTimes) {
      if (outcomeProstate) {
        indices <-
          inputdata$dtime == i &
          inputdata$rx == grp &
          inputdata$eventCens == 0 & inputdata$otherDeath == 0
        eventIndicator <- indices & inputdata$prostateDeath == 1
      } else{
        indices <-
          inputdata$dtime == i & inputdata$rx == grp & inputdata$eventCens == 0
        eventIndicator <- indices & inputdata$otherDeath == 1
      }
      outputHazards[counter] <-
        sum(weightVector[eventIndicator]) / sum(weightVector[indices])
      counter <- counter + 1
    }
    return(outputHazards)
  }

nonParametricCumInc <- function(hazard1, hazard2, competing = FALSE) {
  inc <- rep(NA, length.out = length(cutTimes))
  cumulativeSurvival <- c(1, cumprod((1 - hazard1) * (1 - hazard2)))
  counter <- 1
  for (i in 1:length(cutTimes)) {
    if (!competing) {
      inc[i] <- hazard1[i] * (1 - hazard2[i]) * cumulativeSurvival[i]
    } else{
      inc[i] <- hazard1[i] * cumulativeSurvival[i]
    }
  }
  cumInc <- cumsum(inc)
  return(cumInc)
}

### Separable effects code

calculateCumInc <-
  function(inputData,
           timepts = cutTimes,
           competing = FALSE) {
    cumulativeIncidence <-
      matrix(NA, ncol = length(unique(inputData$patno)), nrow = length(cutTimes))
    # Insert the event probabilities at the first time interval:
    cumulativeIncidence[1,] <-
      inputData[inputData$dtime == 0,]$hazardP * inputData[inputData$dtime ==
                                                             0,]$hazardO
    # Create a matrix compatible with 'cumulativeIncidence' with survival probabilities at each time
    survivalProb <-
      t(aggregate(s ~ patno, data = inputData, FUN = cumprod)$s) #split the long data into subsets per subject
    for (i in 2:length(cutTimes)) {
      subInputDataP <-
        inputData[inputData$dtime == (i - 1),]$hazardP #OBS: dtime starts at 0
      subInputDataO <-
        (1 - inputData[inputData$dtime == (i - 1),]$hazardO) #OBS: dtime starts at 0
      if (!competing)
        cumulativeIncidence[i,] <-
        subInputDataO *  subInputDataP * survivalProb[(i - 1),] # OBS: survivalProb entry i is time point i-1
      else
        cumulativeIncidence[i,] <-
        (1 - subInputDataO) * survivalProb[(i - 1),]
    }
    meanCumulativeIncidence <-
      rowMeans(apply(cumulativeIncidence, MARGIN = 2, cumsum)) #rowMeans(cumulativeIncidence)
    return(meanCumulativeIncidence)
  }

discrete_cuminc_prost <- function(weight_vector, inputdata, grp=0, outcome_y=TRUE,follow_up=1:max_time){
  event_vec <- rep(NA, length.out=length(follow_up))
  counter <- 1 
  # count number of individuals in grp (that is, we cound those who were present at baseline)
  n_grp <- sum(inputdata$dtime==0 & inputdata$rx==grp)
  for(i in follow_up){
    if(outcome_y){
      indices <- inputdata$dtime==i & inputdata$rx == grp & inputdata$eventCens==0 & inputdata$otherDeath==0 
      eventIndicator <- indices & inputdata$prostateDeath==1 
    }else{
      indices <- inputdata$dtime==i & inputdata$rx == grp & inputdata$eventCens==0  
      eventIndicator <- indices &  inputdata$otherDeath==1
    }
    event_vec[counter] <- sum(weight_vector[eventIndicator]) / n_grp
    counter <- counter+1
  }
  output_cuminc <- cumsum(event_vec)
  return(output_cuminc)
}
