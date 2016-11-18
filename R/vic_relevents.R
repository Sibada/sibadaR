
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
    names(soil_hydro_params)[l + nlayer*2] <- sprintf("BUBLE%d", l)
    names(soil_hydro_params)[l + nlayer*3] <- sprintf("BULKDN%d", l)
    names(soil_hydro_params)[l + nlayer*4] <- sprintf("WcrFT%d", l)
    names(soil_hydro_params)[l + nlayer*5] <- sprintf("WpFT%d", l)
  }

  return(soil_hydro_params)
}

#' Create the soil parameters for VIC model.
#'
#' @description Create the soil parameters for VIC model by offering coordinates, elevations
#'              and other datas. Each row of each data must be corresponding to a gridcell.
#'
#' @param coords A data frame containing ongitude and latitude coordinate of the gridcells.
#'               First column must be the longitudes while the second column must be the
#'               latitudes.
#' @param elev A vector containing the elevation of each gridcell.
#'
#' @param soil_props A dataframe like what the function soil_convert needs.
#'
#' @param anprec A vector containing the annual average precipitation of each gridcell.
#'
#' @param avg_T Average soil temperature of each gridcell. Would create from the global
#'              dataset by interpolation if not offered.
#'
#' @param organic Organic parameters. Including the fraction of organic matter, its bulk
#'                density and its partial density for each soil layer.
#'
#' @param Javg Average temperature in July at each gridcell.
#'
#' @param quarzs Quarz content of soil of each gridcell. Would create from the global
#'               dataset by interpolation if not offered.
#'
#' @param nlayer Num of soil layers.
#'
#' @return A data frame containing the soil parameters needed by the VIC model. Can
#'         be write to file use write.table for use.
#'
#'
#' @export
create_soil_params <- function(coords, elev, soil_props, anprec, nlayer=3, avg_T=NA, quarzs=NA, organic=NA, Javg=NA) {
  if(!(nrow(coords) == nrow(elev) & nrow(coords) == nrow(soil_props) & nrow(coords) == anprec)) {
    stop("Input data has different rows.")
  }
  ncell <- nrow(coords)

  lng <- coords[ ,1]
  lat <- coords[ ,2]

  soil_hydraulic <- soil_convert(soil_props, nlayer=nlayer)
  soil_hydraulic <- round(soil_hydraulic, 3)
  offgmt <- lng / 15

  # Create avg_T parameter by interaporation from global soil parameters dataset.
  if(is.na(avg_T)){
    avgT_o <- data.frame(x=VIC_global_soil$LNG, y=VIC_global_soil$LAT, z=VIC_global_soil$AVG_T)
    coordinates(avgT_o) <- ~x + y
    grids <- data.frame(x=lng, y=lat)
    coordinates(grids) <- ~x + y
    gridded(grids) <- TRUE

    avg_T <- data.frame(idw(z~1, avgT_o, grids, nmax=4, maxdist=0.7071068, debug.level=0))[ , 3]
  }

  # Create quarz parameters by interaporation from global soil parameters dataset.
  if(is.na(quarzs)){
    quarz_os <- data.frame(VIC_global_soil$QUARZ1, VIC_global_soil$QUARZ2, VIC_global_soil$QUARZ3)
    quarzs <- data.frame(QUARZ1=rep(0, ncell))
    for(l in 1:nlayer){
      if(l > 3) {
        quarz_o <- quarz_os[ , 3]
      } else {
        quarz_o <- quarz_os[ , l]
      }
      quarz_o <- data.frame(x=VIC_global_soil$LNG, y=VIC_global_soil$LAT, z=quarz_o)
      coordinates(quarz_o) <- ~x + y
      grids <- data.frame(x=lng, y=lat)
      coordinates(grids) <- ~x + y
      gridded(grids) <- TRUE

      quarz <- data.frame(idw(z~1, quarz_o, grids, nmax=4, maxdist=0.7071068, debug.level=0))[ , 3]
      quarzs[sprintf('QUARZ%d', l)] <- quarz
    }
  } else {
    if(ncol(quarzs) != nlayer) {
      stop("Column of quarzs is not equal to layers.")
    }

    names(quarzs) <- paste("QUARZ", 1:nlayer, sep="")
  }

  # Organic part.
  if(!is.na(organic)) {
    if(ncol(organic) != nlayer*3){
      stop("Column of organic is not equal.")
    }
    names(organic) <- c(paste("ORGANIC", 1:nlayer, sep=""),
                        paste("BULKDN_ORG", 1:nlayer, sep=""),
                        paste("PARTDN_ORG", 1:nlayer, sep=""))
  }

  #############################################################################
  # Fill the table of soil parameters.
  #############################################################################
  soil_params <- data.frame(RUN=rep(1, ncell))
  names(soil_params) <- "#RUN"
  soil_params$GRID <- 1:ncell
  soil_params$LAT <- lat
  soil_params$LNG <- lng

  soil_params$INFILT <- 0.25
  soil_params$Ds <- 0.1
  soil_params$Ds_MAX <- 16
  soil_params$Ws <- 0.8

  soil_params$C <- 2
  soil_params <- cbind(soil_params, soil_hydraulic[ , 1:(nlayer*2)])

  for(l in 1:nlayer) {
    soil_params[sprintf('PHI_%d', l)] <- -999
  }

  for(l in 1:nlayer) {
    soil_params[sprintf('MOIST_%d', l)] <- 0.36
  }

  soil_params$ELEV <- round(elev, 1)

  soil_params$DEPTH_1 <- 0.1
  for(l in 2:nlayer) {
    soil_params[sprintf('DEPTH_%d', l)] <- 0.3
  }

  soil_params$AVG_T <- avg_T
  soil_params$DP <- 4.0

  soil_params <- cbind(soil_params, soil_hydraulic[ , (nlayer*2+1):(nlayer*3)])
  soil_params <- cbind(soil_params, quarzs)
  soil_params <- cbind(soil_params, soil_hydraulic[ , (nlayer*3+1):(nlayer*4)])

  for(l in 1:nlayer) {
    soil_params[sprintf('PARKDN%d', l)] <- 2650
  }

  if(!is.na(organic)) {
    soil_params <- cbind(soil_params, organic)
  }

  soil_params$OFF_GMT <- offgmt
  soil_params <- cbind(soil_params, soil_hydraulic[ , (nlayer*4+1):(nlayer*6)])
  soil_params$Z0_SOIL <- 0.010
  soil_params$Z0_SNOW <- 0.001
  soil_params$PRCP <- anprec

  for(l in 1:nlayer) {
    soil_params[sprintf('RESM%d', l)] <- 0.02
  }

  soil_params$FS_ACTV <- 1

  if(!is.na(Javg)) {
    soil_params$JULY_TAVG <- Javg
  }

  return(soil_params)
}

#' Create the veg params file for VIC model.
#'
#' @description Create the vegetration parameter file (generally named "veg_params.txt") with
#'              the output stat data of each single vegetration tile and the gridcell id which
#'              it distributed in. The vegetration data is usually from the dataset of UMD.
#'
#' @param pre_veg_params A vector of integer, the rare two bit means the id ofa kind of
#'                       vegetration conver, and the front bits means the id of the gridcell.
#'                       It's generally generated by Arcgis using the raster calculater.
#' @param out_file Path of the output veg_param file.
#' @param aban A fraction, when a kind of vegeraction cover fraction is little then it, it will
#'             not be considerated in.
#'
#' @export
create_veg_params <- function(pre_veg_params, out_file, aban=0.03) {
  veg_stat <- pre_veg_params

  # Root zone data
  root_pm <- c('0.10 0.05 1.00 0.45 5.00 0.50',
              '0.10 0.05 1.00 0.45 5.00 0.50',
              '0.10 0.05 1.00 0.45 5.00 0.50',
              '0.10 0.05 1.00 0.45 5.00 0.50',
              '0.10 0.05 1.00 0.45 5.00 0.50',
              '0.10 0.10 1.00 0.65 1.00 0.25',
              '0.10 0.10 1.00 0.65 1.00 0.25',
              '0.10 0.10 1.00 0.65 0.50 0.25',
              '0.10 0.10 1.00 0.65 0.50 0.25',
              '0.10 0.10 1.00 0.70 0.50 0.20',
              '0.10 0.10 0.75 0.60 0.50 0.30')

  id <- round(veg_stat/100)
  veg_type <- veg_stat - id*100

  veg_stat <- data.frame(id,veg_type)
  veg_stat <- veg_stat[order(veg_stat$id),]

  vn <- 0
  tl <- nrow(veg_stat)

  type_stat <- rep(0,11)

  out_text <- c()
  r <- 1

  for(i in 1:tl){
    current_id <- veg_stat$id[i]
    tp <- veg_stat$veg_type[i]

    vn <- vn + 1

    if(tp %in% c(1:11))
      type_stat[tp] <- type_stat[tp] + 1

    if(i == tl || veg_stat$id[i + 1] != current_id){
      type_stat <- type_stat / vn

      type_sum <- length(type_stat[type_stat > aban])
      out_text[r] <- paste(current_id,type_sum, sep='\t')
      r <- r + 1

      for(t in 1:11){
        if(type_stat[t] > aban){
          out_text[r] <- paste(t,round(type_stat[t],5),root_pm[t], sep='\t')
          r <- r + 1
        }
      }
      type_stat <- rep(0,11)
      vn <- 0
    }
  }

  write.table(out_text, out_file, row.names=F, col.names=F, quote=F)
}

