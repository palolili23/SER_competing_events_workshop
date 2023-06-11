# Causal inference with competing events

A competing risk event is any event that ensures the outcome of interest cannot subsequently occur. For example, in a study where prostate cancer death is the primary outcome, a fatal stroke is a competing event because an individual cannot die of cancer once they have died of stroke. When competing events are present, many possible definitions of a causal effect may be considered. Choosing a causal effect of practical interest requires understanding the interpretation of different counterfactual contrasts and the assumptions needed to identify them using the study data and subject matter knowledge. This workshop will introduce participants to a counterfactual framework for causal inference in the face of competing events. Participants will learn how to articulate and interpret different types of causal effects when competing events are present, and approaches to estimating them under transparent assumptions with the aid of causal diagrams. In part I, we cover counterfactual contrasts of popular parameters from the competing risks literature, including contrasts of cause- specific and subdistribution hazards, and cause-specific cumulative incidences and their relation to total and controlled direct effects from the mediation literature. In part II, we introduce the separable effects, new causal effect definitions that may be of particular clinical relevance in competing events settings. Theoretical concepts will be illustrated via practical examples and R code provided.

## Workshop Details

### In Person:
June 13, 2023
1:00 pm – 5:00 pm

### Virtual:
July 6, 2023
10:00 am – 2:00pm MT

### Instructors:

- Jessica Young
- Mats Stensrud
- L. Paloma Rojas-Saunero

## Repository content

- `R` folder: Includes R scripts and data 

- [Guided scripts](https://palolili23.github.io/SER_competing_events_workshop/R/index.html)

## Disclaimer:

This workshop was based on the `R code` from [Young JG et al. Statistics in Medicine 2020](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471) and [Stensrud MJ et al. Journal of the American Statistical Association 2022](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1765783). Since both papers use the same data source and have similar procedures of data cleaning, we have unified several preliminary steps, specially in the `Data cleaning and setup` section, for better readability and comprehension in terms of the workshop. However, this means that results will not be identical to the results of these manuscripts.

Some of the differences that arise in this version, compared to the code by Stensrud et al. are:

- DES is A=1 and placebo A=0

- `cutTimes` is 60, instead of 101

- `prostRed$hgBinary <- prostRed$hg < 12` instead of `prostRed$hg < 10` 

- Since we include the following code:

```
longProstRed$prostateDeath[longProstRed$eventCens == 1] <- NA

longProstRed$otherDeath[longProstRed$eventCens == 1] <- NA

longProstRed$prostateDeath[longProstRed$otherDeath == 1] <- NA

```

`ipw_d <- cum_pred_O_1 / cum_pred_O_0` was used instead of `ipw_d[!t0] <- cum_pred_O_1[!t0] / cum_pred_O_0[!t0]`

- Several variables names have been modified compared to the code by Young et al. and by Stensrud et al.

- Hand calculations for the interaction terms between the exposure and time of follow-up were removed and indicators were included in the models. 
