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
#'  * W: Wild-type allele
#'  * M: MEREA allele
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
cubeMEREA <- function(rM = 0, rW = 0, Teff = 1.0, eta = NULL, phi = NULL,
                      omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## safety checks
  if(any(c(rM,rW,Teff)<0) || any(c(rM,rW,Teff)>1)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('ZW', 'ZZ', 'MZ', 'MM', 'MW') #changed genotypes to merea alleles #this is every possible genotype #took out resistant alleles
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix
  
  ## fill tMatrix with probabilities 
  #( 'ZW'(female), 'ZZ'(male), 'MZ'(male), 'MM'(male), 'MW'(female)) #Merea is Z-linked #M and Z are both variants of the Z chromosome
  #Use explicit syntax
  #explicit pros: can check line by line #explicit cons: risk of typos and will not be assigned to correct genotype #therefore implicit is made to reduce typos and bugs
  #later: will need to specify sex-specific alleles
  
  #remember: female on the left and male on the right
  tMatrix['ZW','ZZ', c('ZW', 'ZZ')] <- #100% wildtype #example of explicit assigning syntax #mathemetical outcomes come after function arrow and are separated by commas
  
  tMatrix['ZW','MZ', c('MZ', 'ZZ', 'MW', 'ZW')] <-
  tMatrix['ZW','MM', c('MZ', 'MW')] <-
  
  tMatrix['MW','ZZ', c('MZ', 'ZW')] <-
  tMatrix['MW','MZ', c('MW', 'MM', 'MZ')] <-
  tMatrix['MW','MM', c('MM', 'MW')] <-
  
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