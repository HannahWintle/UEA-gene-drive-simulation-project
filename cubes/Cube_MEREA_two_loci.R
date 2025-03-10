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
#   March 2025
#   cnm20syu@uea.ac.uk
#
#   Inheritance Cube: MEREA (Maternal Effect RECESSIVE Embryonic Arrest)
#
#   This drive contains two different transgenes on the MEREA allele. Exu-Cas9 
#   (acting as the toxin) remains the same, but the rescue genes are split:
#    * Ma: Contains rescuegenepromoter-tTAV
#    * Mb: Contains teto-recoded.rescuegene
#
#   Rescue Conditions:
#   - Males must inherit both copies of each transgene to survive (MaMb genotype).
#   - Instead of a simple homozygous rescue system (as in MEREA), males must be 
#     trans-homozygous.
#   - Males with MaMa or MbMb **are not viable**.
#   - Only MaMb males **survive**.
#
#   Male Genotype Key:
#    * ZZ      : Wild-type male
#    * MaZ     : Male with only Ma transgene 
#    * MbZ     : Male with only Mb transgene 
#    * MaMa    : Male with two Ma alleles
#    * MbMb    : Male with two Mb alleles
#    * MaMb    : Male with both Ma and Mb
#    * RZ      : Male with resistance allele
#    * RMa     : Male with resistance + Ma
#    * RMb     : Male with resistance + Mb
#    * RR      : Male with two resistance alleles
#   
#   Female Genotype Key:
#    * ZW   : Wild-type female
#    * MaW  : Female with only Ma transgene
#    * MbW  : Female with only Mb transgene
#    * RW   : Female with resistance allele
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
#' @return Named list containing the inheritance cube, transition matrix, genotypes, 
#' wild-type allele, and all genotype-specific parameters.
#' @export
#
###############################################################################

cubeMEREA_2L <- function(rM = 0, Teff = 1.0, eta = NULL, phi = NULL,
                                           omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## Safety checks
  if(any(c(rM, Teff) < 0) || any(c(rM, Teff) > 1)){
    stop("Parameters are rates. 0 <= x <= 1")
  }
  
  ## Define genotype list
  gtype <- c('ZW', 'ZZ', 'MaZ', 'MbZ', 'MaMb', 'MaMa', 'MbMb', 'RZ', 'RMa', 'RMb', 'RR', 'MaW', 'MbW', 'RW')
  size <- length(gtype)
  
  ## Initialize transition matrix
  tMatrix <- array(data = 0, dim = c(size, size, size), dimnames = list(gtype, gtype, gtype))
  
  ## Fill tMatrix with inheritance probabilities
  tMatrix['ZW', 'ZZ', c('ZW', 'ZZ')] <- c(1, 1) / 2  # Wild-type cross
  
  tMatrix['ZW', 'MaZ', c('MaZ', 'ZZ', 'MaW', 'ZW')] <- c(1-rM, 1, 1-rM, 1) / 4
  tMatrix['ZW', 'MbZ', c('MbZ', 'ZZ', 'MbW', 'ZW')] <- c(1-rM, 1, 1-rM, 1) / 4
  tMatrix['ZW', 'MaMa', c('MaZ', 'MbZ', 'MaW', 'MbW')] <- c((1-rM), (1-rM), (1-rM), (1-rM)) / 4
  tMatrix['ZW', 'MbMb', c('MaZ', 'MbZ', 'MaW', 'MbW')] <- c((1-rM), (1-rM), (1-rM), (1-rM)) / 4
  tMatrix['ZW', 'MaMb', c('MaZ', 'MbZ', 'MaW', 'MbW')] <- c((1-rM), (1-rM), (1-rM), (1-rM)) / 4
  tMatrix['ZW', 'RZ', c('ZZ', 'RZ', 'ZW', 'RW')] <- c(1, 1, 1, 1) / 4
  tMatrix['ZW', 'RMa', c('ZZ', 'RZ', 'ZW', 'RW')] <- c(1, 1, 1, 1) / 4
  tMatrix['ZW', 'RMb', c('ZZ', 'RZ', 'ZW', 'RW')] <- c(1, 1, 1, 1) / 4
  tMatrix['ZW', 'RR', c('RZ', 'RW')] <- c(1, 1) / 2
  
  tMatrix['MaW','ZZ', c('MZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2
  tMatrix['MaW','MaZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['MaW','MbZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['MaW','MaMa', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MaW','MbMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MaW','MaMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MaW','RZ', c('MZ', 'MR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4
  tMatrix['MaW','RMa', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['MaW','RMb', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['MaW','RR', c('MR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2
  
  tMatrix['MbW','ZZ', c('MZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2
  tMatrix['MbW','MaZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['MbW','MbZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['MbW','MaMa', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MbW','MbMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MbW','MaMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['MbW','RZ', c('MZ', 'MR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4
  tMatrix['MbW','RMa', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['MbW','RMb', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['MbW','RR', c('MR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2
  
  tMatrix['RW','ZZ', c('MZ', 'ZW', 'RZ')] <- c(1-rM, 1, rM)/2
  tMatrix['RW','MaZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['RW','MbZ', c('MW', 'MM', 'MZ', 'ZW', 'RZ', 'RW', 'RR', 'MR')] <- c(1-rM, (1-rM)^2, 1-rM, 1, rM, rM, rM^2, 2*(1-rM)*rM)/4
  tMatrix['RW','MaMa', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['RW','MbMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['RW','MaMb', c('MM', 'MW', 'RR', 'RW', 'MR')] <- c((1-rM)^2, (1-rM), rM^2, rM, 2*(1-rM)*rM)/2
  tMatrix['RW','RZ', c('MZ', 'MR', 'ZW', 'RW', 'RZ', 'RR')] <- c((1-rM), (1-rM), 1, 1, rM, rM)/4
  tMatrix['RW','RMa', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['RW','RMb', c('MM', 'MR', 'MW', 'RW', 'RR')] <- c((1-rM)^2, (1-rM)*(1+2*rM), (1-rM), (1+rM), rM*(1+rM))/4
  tMatrix['RW','RR', c('MR', 'RW', 'RR')] <- c((1-rM), 1, rM)/2
  
  ## Define viability mask
  viabilityMask <- array(data = 1, dim = c(size, size, size), dimnames = list(gtype, gtype, gtype))
  viabilityMask[c('MaMa', 'MbMb', 'RZ', 'RMa', 'RMb')] <- 0  # Non-viable combinations
  
  ## Genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  
  ## Define sex modifiers
  # 0 is male and 1 is female
  sexes <- c(
    "ZW" = 1, 
    "ZZ" = 0,
    "MaZ" = 0, 
    "MbZ" = 0, 
    "MaMb" = 0,
    "MaMa" = 0, 
    "MbMb" = 0, 
    "RZ" = 0, 
    "RMa" = 0, 
    "RMb" = 0, 
    "RR" = 0,
    "MaW" = 1, 
    "MbW" = 1, 
    "RW" = 1
  )
  phi <- sexes
  
  ## Return cube structure
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c("ZZ", "ZW"),
    eta = modifiers$eta,
    phi = phi,  # Assign sex ratios
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = c("MaMb", "MaMa", "MbMb", "MaW", "MbW")
  ))
}
