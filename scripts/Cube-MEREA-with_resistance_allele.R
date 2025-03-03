###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   MEDEA
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   August 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: MEDEA (Maternal Effect Dominant Embryonic Arrest)
#'
#' This function creates an inheritance cube to model a MEDEA drive system. This
#' system was first discovered in flour beetles. It biases inheritance by expressing
#' a maternal toxin such that offspring die unless they express a zygotic antidote. \cr
#' This drive has 3 alleles at 1 locus:
#'  * ZZ: Wild-type male
#'  * ZW: Wild-type female
#'  * M: MEREA allele
#'  * R: Resistance allele
#'
#' @param rM Breakdown of MEDEA allele, no homing/toxin/antidote, M -> R conversion
#' @param rW De novo resistance generation, W -> R conversion
#' @param Teff Efficacy of the toxin
#' @param eta Genotype-specific mating fitness
#' @param phi Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF Genotype-specific female pupatory success
#' @param xiM Genotype-specific male pupatory success
#' @param s Genotype-specific fractional reduction(increase) in fertility
#'
#' @return Named list containing the inheritance cube, transition matrix, genotypes, wild-type allele,
#' and all genotype-specific parameters.
#' @export
cubeMEREA <- function(rM = 0, Teff = 1.0, eta = NULL, phi = NULL,
                      omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## safety checks
  if(any(c(rM,rW,Teff)<0) || any(c(rM,rW,Teff)>1)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }
  
  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('ZW', 'ZZ', 'MZ', 'MM', 'MW', 'RW', 'RZ', 'MR', 'RR') #changed genotypes to merea alleles #this is every possible genotype
  size <- length(gtype) #because I use it several times
tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix
  
  ## fill tMatrix with probabilities 
  #( 'ZW'(female), 'ZZ'(male), 'MZ'(male), 'MM'(male), 'MW'(female), 'RW'(female), 'RZ'(male), 'MR'(male), 'RR'(male)) #Merea is Z-linked #M and Z are both variants of the Z chromosome
  #Use explicit syntax
  #explicit pros: can check line by line #explicit cons: risk of typos and will not be assigned to correct genotype #therefore implicit is made to reduce typos and bugs
  #later: will need to specify sex-specific alleles
  
  #remember: female on the left and male on the right
  tMatrix['ZW','ZZ', c('ZW', 'ZZ', )] <- c(1, 1)/2
    #100% wildtype #example of explicit assigning syntax #mathematical outcomes come after function arrow and are separated by commas
    
  tMatrix['ZW','MZ', c('MZ', 'ZZ', 'MW', 'ZW', 'RZ', 'RW')] <- c(1-rM, 1, 1-rM, 1, rM, rM)/4
  tMatrix['ZW','MM', c('MZ', 'MW', 'RZ', 'RW')] <- c((1-rM)^2, (1-rM)^2, 1-((1-rM)^2), 1-((1-rM)^2))/2
  tMatrix['ZW','MR', c('MZ', 'RZ', 'MW', 'RW')] <- c(1-rM, (1-rM)+rM*2, 1-rM, (1-rM)+rM*2)/4
  tMatrix['ZW','RZ', c('ZZ', 'RZ', 'ZW', 'RW')] <- c(1, 1, 1, 1)/4
  tMatrix['ZW','RR', c('RZ', 'RW')] <- c(1, 1)/2
  
  tMatrix['MW','ZZ', c('MZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2
  tMatrix['MW','MZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['MW','MM', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MW','MR', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['MW','RZ', c('MZ', 'MR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4
  tMatrix['MW','RR', c('MR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2
  
  tMatrix['RW','ZZ', c('RZ', 'ZW',)] <- c(1, 1)/2
  tMatrix['RW','MZ', c('MR', 'RZ', 'MW', 'ZW', 'RR', 'RW')] <- c((1-rM), 1, (1-rM), 1, rM, rM)/4
  tMatrix['RW','MM', c('MR', 'MW', 'RR', 'RW')] <- c((1-rM), (1-rM), rM, rM)/2
  tMatrix['RW','MR', c('RW', 'MR', 'RR', 'MW')] <- c((1+rM), (1-rM), (1+rM), (1-rM))/4
  tMatrix['RW','RZ', c('RR', 'RZ', 'RW', 'ZW')] <- c(1, 1, 1, 1)/4
  tMatrix['RW','RR', c('RR', 'RW')] <- c(1, 1)/2
  
  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}
  
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors
  
  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))
  
  ## fill mother/offspring specific death, then muliply by efficacy of toxin and antidote
  for(slice in 1:size){
    viabilityMask[c('WM', 'MR'),slice, ] <- matrix( c( 1-Teff, 1, 1-Teff,1, 1, 1), nrow = 2, ncol = size, byrow = TRUE )
    viabilityMask['MM',slice, ] <- c( 1-2*Teff+Teff^2, 1, 1-2*Teff+Teff^2, 1, 1, 1)
  }
  
  # for loop may be unnecessary
  
  # All offspring of MW females are dead regardless of father
  viabilityMask['MW',1:9,1:9] <- matrix( rep(1-Teff,size), nrow = 1, ncol = size, byrow = TRUE )
  # MM offspring always survive, even if born to MW females
  viabilityMask['MW', 1:9, 'MM'] <- matrix(rep(1, size), nrow = 1, ncol = size, byrow = TRUE)
  
  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  
  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "WW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "MM"
  ))
}