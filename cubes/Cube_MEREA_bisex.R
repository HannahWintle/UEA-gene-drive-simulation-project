###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   MEREA
#   Hannah Wintle, Philip Leftwich
#   February 2025
#   cnm20syu@uea.ac.uk
#
#   Inheritance Cube: MEREA (Maternal Effect RECESSIVE Embryonic Arrest)
#
#   This function creates an inheritance cube to model a MEREA drive system.
#   MEREA functions via a maternal-effect toxin that is deposited into the eggs, 
#   causing embryonic lethality unless the offspring inherit a zygotic antidote. 
#   The antidote is only effective when inherited from both parents, making the 
#   system recessive.
#
#   This drive operates with a single allele (MEREA) at one locus.
#
#   Male Genotype Key:
#    * ZZ  : Wild-type male
#    * MZ  : Male with one MEREA allele
#    * MM  : Male with two MEREA alleles
#    * RZ  : Male with resistance allele (no MEREA)
#    * RM  : Male with resistance + MEREA allele
#    * RR  : Male with two resistance alleles (wild-type behaviour)
#
#   Female Genotype Key:
#    * ZW  : Wild-type female
#    * MW  : Female with one MEREA allele
#    * RW  : Female with resistance allele (no MEREA)
#
#' @param rM    Breakdown of MEREA allele, no homing/toxin/antidote, M -> R conversion
#' @param Teff  Efficacy of the toxin
#' @param eta   Genotype-specific mating fitness
#' @param phi   Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF   Genotype-specific female pupatory success
#' @param xiM   Genotype-specific male pupatory success
#' @param s     Genotype-specific fractional reduction (increase) in fertility
#'
#' @return Named list containing the inheritance cube, transition matrix, 
#' genotypes, wild-type allele, and all genotype-specific parameters.
#' @export
#
###############################################################################

cubeMEREA <- function(rM = 0, Teff = 1.0, eta = NULL, phi = NULL,
                      omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## safety checks
  if(any(c(rM, Teff)<0) || any(c(rM, Teff)>1)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }
  
  ## define matrices
  ## Matrix Dimensions Key: [female Genotype, male Genotype, offspring Genotype]
  gtype <- c('ZW', 'ZZ', 'MZ', 'MM', 'MW', 'RW', 'RZ', 'MR', 'RR') #changed genotypes to merea alleles #this is every possible genotype
  size <- length(gtype) #because I use it several times
tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix
  
  ## fill tMatrix with probabilities 
  #( 'ZW'(female), 'ZZ'(male), 'MZ'(male), 'MM'(male), 'MW'(female), 'RW'(female), 'RZ'(male), 'MR'(male), 'RR'(male)) #Merea is Z-linked #M and Z are both variants of the Z chromosome
  #Use explicit syntax
  #explicit pros: can check line by line #explicit cons: risk of typos and will not be assigned to correct genotype #therefore implicit is made to reduce typos and bugs
  #later: will need to specify sex-specific alleles
  
  #remember: female on the left and male on the right
  tMatrix['ZW','ZZ', c('ZW', 'ZZ')] <- c(1, 1)/2
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
  
  tMatrix['RW','ZZ', c('RZ', 'ZW')] <- c(1, 1)/2
  tMatrix['RW','MZ', c('MR', 'RZ', 'MW', 'ZW', 'RR', 'RW')] <- c((1-rM), 1, (1-rM), 1, rM, rM)/4
  tMatrix['RW','MM', c('MR', 'MW', 'RR', 'RW')] <- c((1-rM), (1-rM), rM, rM)/2
  tMatrix['RW','MR', c('RW', 'MR', 'RR', 'MW')] <- c((1+rM), (1-rM), (1+rM), (1-rM))/4
  tMatrix['RW','RZ', c('RR', 'RZ', 'RW', 'ZW')] <- c(1, 1, 1, 1)/4
  tMatrix['RW','RR', c('RR', 'RW')] <- c(1, 1)/2
  
  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
 # boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  #for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}
  
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors
  
  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))
  
  # Kill all offspring of MW females by default
  viabilityMask['MW', 1:size, 1:size] <- 0
  
  # Allow survival only for specific genotypes
  allowed_offspring <- c('MM', 'MR', 'RZ', 'RW', 'RR')
  for (offspring in allowed_offspring) {
    viabilityMask['MW', 1:size, offspring] <- 1
  }
  
  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  
  ##sex modifiers
  #0 is male and 1 is female
  sexes <- c(
    "ZW" = 1,  # Female
    "ZZ" = 0,  # Male
    "MZ" = 0,  # Male
    "MM" = 0,  # Male
    "MW" = 1,  # Female
    "RW" = 1,  # Female
    "RZ" = 0,  # Male
    "MR" = 0,  # Male
    "RR" = 0   # Male
  )
  
  # Define phi (sex ratio at emergence)
  phi <- sexes
  
  ## put everything into a labeled list to return
  #need to change wildtype and release type
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c("ZZ", "ZW"),
    eta = modifiers$eta,
    phi = phi, #assign sex ratios
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = c("MM", "MW")
  ))
}
