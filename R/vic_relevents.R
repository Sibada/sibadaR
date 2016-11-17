
#' Create the hydraulic parameters of soil from soil property parameters for VIC model.
#'
#' @description
#'
#' @param soil_props A dataframe containing the soil hydraulic parameters, including
#'                   volume fraction of gravity (fraction, not percentage), weight fraction
#'                   of sand, silt and clay, reference bulk density, and weight fraction of
#'                   organic matter, total 6 property parameters. The first 6 columns are
#'                   parameters of the first soil layer, and the next 6 columns are of the
#'                   second soil layer, and so on.
#' @param nlayer Num of layers of each gridcell.
#'
#' @return A dataframe containing the soil hydraulic parameters for each layer of each
#'         gridcell. Including EXPT (Exponent), Ksat (Saturated hydrologic conductivity),
#'         BUBLE (Bubbling pressure), BULKDN (Bulk density), WcrFT (Fractional soil
#'         moisture content at the critical point) and WpFT (Fractional soil moisture
#'         content at the wilting point).
#'
#' @details This function is using the Saxton's formulars. If the soil property parametersis
#'          HWSD, The property parameters are corresponding to T_GRAVEL, T_SAND, T_SILT, T_CLAY,
#'          T_REF_BULK_DENSITY, and T_OC.
#'
#' @references Saxton, K. E, Rawls, W. J. Soil Water Characteristic Estimates by Texture
#'             and Organic Matter for Hydrologic Solutions[J]. Soil Science Society of
#'             America Journal, 2006, 70(5):1569-1578.
#'
#' @export
soil_convert <- function(soil_props, nlayer=3) {
  if(ncol(x) < nlayer * 7) {
    stop("Soil properties is insufficient.")
  }

  soil_hydro_params <-  data.frame(matrix(0, nrow=nrow(soil_props), ncol=6*nlayer))

  for(l in 1:nlayer) {

    S <- soil_props[ , 2+(l-1)*6]
    Si <- soil_props[ , 3+(l-1)*6]
    C <- soil_props[ , 4+(l-1)*6]
    OM <- soil_props[ , 6+(l-1)*6]
    rhoDF <- soil_props[ , 5+(l-1)*6]
    Rv <- soil_props[ , 1+(l-1)*6]

    fb <- S + Si + C
    S <- S/fb
    C <- C/fb

    theta33t <- -0.251*S+0.195*C+0.011*OM+0.006*(S*OM)-0.027*(C*OM)+0.452*S*C+0.299
    theta33 <- theta33t + (1.283*theta33t*theta33t-0.374*theta33t-0.015)
    theta1500t <- -0.024*S+0.487*C+0.006*OM+0.005*S*OM-0.013*C*OM+0.068*S*C+0.031
    theta1500 <- theta1500t + (0.14*theta1500t-0.02)

    thetas_33t <- 0.278*S+0.034*C+0.022*OM-0.018*S*OM-0.027*C*OM-0.584*S*C+0.078
    thetas_33 <- thetas_33t + (0.636*thetas_33t-0.107)

    thetas <- theta33 + thetas_33-0.097*S+0.043
    rhoN <- (1-thetas)*2.685

    DF <- rhoDF/rhoN
    thetasDF <- 1-(rhoDF/2.685)
    theta33DF <- theta33-0.2*(thetas - thetasDF)
    thetas_33DF <- thetasDF - theta33DF

    phiet <- -21.674*S-27.932*C-81.975*thetas_33DF+71.121*S*thetas_33DF+8.294*C*thetas_33DF+14.05*S*C+27.161
    phie <- phiet + 0.02*phiet*phiet-0.113*phiet-0.7

    B <- (log(1500) - log(33))/(log(theta33DF)-log(theta1500))
    A <- exp(log(33) +B*log(theta33DF))
    phai1500_33 <- A*thetas**(-B)
    phi33_e <- 33-(thetas-theta33)*(33-phie)/(thetas-theta33)

    lamda <- 1/B
    Ks <- 1930*(thetas_33DF)**(3-lamda)

    phie[phie<0] <- min(phie[phie>0])
    phie_cm <- phie * 10.19368

    ######################################################
    EXPT <- 3+2/lamda
    Ksat <- Ks * 24  # Convert from mm/h to mm/day
    BUBLE <- phie_cm
    BULKDN <- rhoDF * 1000
    WcrFT <- theta33DF/thetasDF
    WpFT <- theta1500/thetasDF

    WcrFT[WcrFT > 1] <- 1.0
    WpFT[WpFT > 1] <- 1.0

    soil_hydro_params[, l] <- EXPT
    soil_hydro_params[, l + nlayer] <- Ksat
    soil_hydro_params[, l + nlayer*2] <- BUBLE
    soil_hydro_params[, l + nlayer*3] <- BULKDN
    soil_hydro_params[, l + nlayer*4] <- WcrFT
    soil_hydro_params[, l + nlayer*5] <- WpFT

    names(soil_hydro_params)[l] <- sprintf("EXPT_%d", l)
    names(soil_hydro_params)[l + nlayer] <- sprintf("Ksat_%d", l)
    names(soil_hydro_params)[l + nlayer*2] <- sprintf("BUBLE_%d", l)
    names(soil_hydro_params)[l + nlayer*3] <- sprintf("BULKDN_%d", l)
    names(soil_hydro_params)[l + nlayer*4] <- sprintf("WcrFT_%d", l)
    names(soil_hydro_params)[l + nlayer*5] <- sprintf("WpFT_%d", l)
  }

  return(soil_hydro_params)
}

#create_soil_params <- function(coords, elev, soil_props, anprec, T_avg=NA, quarz=NA) {

#}
