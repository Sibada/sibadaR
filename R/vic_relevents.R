
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
  if(ncol(soil_props) < nlayer * 6) {
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
#' @import gstat
#' @export
create_soil_params <- function(coords, elev, soil_props, anprec, nlayer=3, avg_T=NA, quarzs=NA, organic=NA, Javg=NA) {

  if(!(nrow(coords) == length(elev) & nrow(coords) == nrow(soil_props) & nrow(coords) == length(anprec))) {
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

#' Create the flow direction of each gridcell for VIC routing model.
#'
#' @description Create the flow direction data of each VIC gridcell for the routing model
#'              of VIC by providing a high resolution flow accumulation data (about 90m
#'              usually caculated from SRTM DEM by ArcGIS) and the domain file of VIC model.
#'              WARNING: This is a simple heuristic method and more fit for small gridcell.
#'              It's usually needs further manual calibration by set_direc and plot_flow_direc.
#'
#' @param flow_file An netCDF4 file of the flow accumulation data from high resolution DEM
#'                  (about 90m). Must containing all the area.
#' @param domain_file Domain file (netCDF4 format) of VIC model.
#' @param arc_code If the direction code is ArcInfo type. If FALSE, it would set to 1 to 8.
#'
#' @param diag_tsh Diag threshould. Determine if flow to the diag gridcell.
#'
#' @return A list including the meta parameters (num of rows and columns, size of the cells,
#'         x and y corner) and the matrix of the flow direction data.
#'
#' @export
create_flow_direction <- function(flow_file, domain_file, arc_code=TRUE, diag_tsh=0.382) {
  if(arc_code) {
    dir_code <- c(64, 128, 1, 2, 4, 8, 16, 32)
  } else {
    dir_code <- 1:8
  }

  grids_nc <- nc_open(domain_file)
  grids_var <- ncvar_get(grids_nc, grids_nc$var[[1]])
  glons <- grids_nc$dim$lon$vals
  glats <- grids_nc$dim$lat$vals
  nc_close(grids_nc)

  nrows <- length(glons)
  ncols <- length(glats)
  csizex <- abs(mean(glons[2:nx]-glons[1:(nx-1)]))
  csizey <- abs(mean(glats[2:ny]-glats[1:(ny-1)]))

  flow_nc <- nc_open(flow_file)

  direc <- grids_var
  direc[!is.na(direc)] <- 0
  for(col in 1:ncols) {
    for(row in 1:nrows) {
      if(is.na(direc[row, col]) | direc[row, col] > 0) next

      glon <- glons[row]
      glat <- glats[col]

      mf <- get_border(flow_nc, glon, glat, csizex, csizey)
      mfd <- mf[1]
      mfs <- mf[2]

      nextrow <- row
      nextcol <- col
      if(mfd == 1) {
        nextcol <- col - 1
      } else if(mfd == 2) {
        nextrow <- row + 1
      } else if(mfd == 3) {
        nextcol <- col + 1
      } else {
        nextrow <- row - 1
      }

      if(mfd == 1) {
        old <- dir_code[8]
        osd <- dir_code[1]
        ord <- dir_code[2]
      } else if(mfd == 2) {
        old <- dir_code[2]
        osd <- dir_code[3]
        ord <- dir_code[4]
      } else if(mfd == 3) {
        old <- dir_code[4]
        osd <- dir_code[5]
        ord <- dir_code[6]
      } else if(mfd == 4) {
        old <- dir_code[6]
        osd <- dir_code[7]
        ord <- dir_code[8]
      }

      if(mfd == 1 & nextcol < 1 | mfd == 2 & nextrow > nrows |
         mfd == 3 & nextcol > ncols | mfd == 4 & nextrow < 1) {
        direc[row, col] <- osd
        next
      }

      nglon <- glons[nextrow]
      nglat <- glats[nextcol]

      nmf <- get_border(flow_nc, nglon, nglat, csizex, csizey)
      nmfd <- nmf[1]
      nmfs <- nmf[2]

      nmfd <- (nmfd + 4 - mfd) %% 4 + 1

      if(nmfd == 3) {
        direc[row, col] <- osd
        direc[nextrow, nextcol] <- osd
      } else if(nmfd == 4 & mfs < diag_tsh & nmfs < diag_tsh) {
        direc[row, col] <- old
      } else if(nmfd == 2 & mfs > 1-diag_tsh & nmfs > 1-diag_tsh) {
        direc[row, col] <- ord
      } else {
        direc[row, col] <- osd
      }
    }
  }
  direc <- t(direc)
  nc_close(flow_nc)
  arcgrid <- list('ncols'=nrows,
                  'nrows'=ncols,  # The row and col is inversed of nc file.
                  'xllcorner'=min(glons)-csizex/2,
                  'yllcorner'=min(glats)-csizey/2,
                  'cellsize'=(csizex+csizey)/2,
                  'grid'=direc)
  return(arcgrid)
}

# Assistant function.
get_border <- function(nc, x, y, xsize, ysize) {
  flons <- nc$dim$x$vals
  flats <- nc$dim$y$vals
  l <- x - xsize/2
  r <- x + xsize/2
  t <- y + ysize/2
  b <- y - ysize/2
  cols <- which(flons > l & flons < r)
  rows <- which(flats > b & flats < t)
  ncol <- length(cols)
  nrow <- length(rows)
  lf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[1]), c(1, nrow))
  rf <- ncvar_get(nc, nc$var[[1]], c(cols[ncol], rows[1]), c(1, nrow))
  tf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[1]), c(ncol, 1))
  bf <- ncvar_get(nc, nc$var[[1]], c(cols[1], rows[120]), c(ncol, 1))
  lf[is.na(lf)] <- 0
  rf[is.na(rf)] <- 0
  tf[is.na(tf)] <- 0
  bf[is.na(bf)] <- 0

  lf <- rev(lf)
  bf <- rev(bf)
  maxps <- c(which.max(tf), which.max(rf), which.max(bf), which.max(lf))
  maxs <- c(max(tf), max(rf), max(bf), max(lf))
  maxside <- which.max(maxs)

  if(maxside == 1 | maxside == 3) {
    maxp <- maxps[maxside]/ncol
  } else {
    maxp <- maxps[maxside]/nrow
  }
  return(c(maxside, maxp))
}

#' Convert the grid data (save as matrix) to points.
#'
#' @description Conver the grid data to points, showed as a table, with the columns of
#'              coordinates for each point and a column for their value.
#'
#' @param grid A matrix of the gridded data.
#' @param x X coordinates for each column of the grid.
#' @param y Y coordinates for each column of the grid.
#' @param csize Size of each gridcell.
#' @param xcor X coordinate of the southwest corner of the grid.
#' @param xcor Y coordinate of the southwest corner of the grid.
#'
#' @return A table of the points, including x and y coordinate and the values, converted
#'         from the grid.
#' @export
grid2points <- function(grid, x=NULL, y=NULL, csize=NULL, xcor=NULL, ycor=NULL, NA_value=NULL, y_rev=FALSE) {
  if((is.null(x) | is.null(y)) & (is.null(csize) | is.null(xcor) | is.null(ycor))) {
    stop("Must provide x, y or csize, xcor, ycor")
  }

  if(is.null(x) | is.null(y)) {
    x <- (1:ncol(grid)) * csize - csize/2 + xcor
    y <- (nrow(grid):1) * csize - csize/2 + ycor
  }
  if (y_rev) y <- rev(y)

  xs <- c()
  ys <- c()
  vs <- c()
  if(!is.null(NA_value)) grid[grid == NA_value] = NA
  for(row in 1:nrow(grid)) {
    for(col in 1:ncol(grid)) {
      if(is.na(grid[row, col])) next

      xs <- append(xs, x[col])
      ys <- append(ys, y[row])
      vs <- append(vs, grid[row, col])
      }
  }
  points <- data.frame(xs, ys, vs)
  names(points) <- c('x', 'y', 'value')
  return(points)
}


#' Convert the point data (a table including columns of the coordinates and value) to
#' grid.
#'
#' @description Conver the point data, usually a table with columns including coordinates
#'              and values to grids.
#'
#' @param grid A matrix of the gridded data.
#' @param x Which column store the x cordinate.
#' @param y Which column store the y cordinate.
#' @param val Which column store the point value.
#' @param out_file The output file path. If provide, it will write as an ArcInfo file
#'                 without return.
#'
#' @return A list including the meta parameters (num of rows and columns, size of the cells,
#'         x and y corner) and the matrix of the gridcell value.
#' @export
points2grid <- function(points, x=NULL, y=NULL, val=NULL, out_file=NULL) {
  acc <- 6

  if(is.null(x)) {
    xs <- points[ , 1]
  } else {
    xs <- points[ , x]
  }
  if(is.null(y)) {
    ys <- points[ , 2]
  } else {
    ys <- points[ , y]
  }
  print(val)
  if(is.null(val)) {
    v <- points[ , 3]
  } else {
    v <- points[ , val]
  }

  ux <- sort(unique(round(xs, acc)))
  uy <- sort(unique(round(ys, acc)))

  ncol <- length(ux)
  nrow <- length(uy)

  itvx <- unique(round(ux[2:ncol] - ux[1:(ncol-1)], acc))
  itvy <- unique(round(uy[2:nrow] - uy[1:(nrow-1)], acc))
  if(length(itvx) > 1 | length(itvy) > 1 |
     itvx[1] != itvy[1]) {
    stop("Intervals of points are not equal. Cannot be gridded.")
  }
  cellsize <- itvx

  xcor <- min(ux) - cellsize/2
  ycor <- min(uy) - cellsize/2

  cols <- (xs - xcor + cellsize/2)/cellsize
  rows <- (ys - ycor + cellsize/2)/cellsize
  rows <- max(rows) - rows + 1

  grid <- matrix(nrow=nrow, ncol=ncol)
  for(p in 1:nrow(points)) {
    grid[rows[p], cols[p]] <- v[p]
  }
  if(is.null(out_file)){
    arcgrid <- list('ncols'=ncol,
                    'nrows'=nrow,
                    'xllcorner'=xcor,
                    'yllcorner'=ycor,
                    'cellsize'=cellsize,
                    'grid'=grid)
    return(arcgrid)
  } else {
    meta <- c(paste('ncols         ', ncol),
              paste('nrows         ', nrow),
              paste('xllcorner     ', xcor),
              paste('yllcorner     ', ycor),
              paste('cellsize      ', cellsize),
              paste('NODATA_value  ', '-9999'))
    writeLines(meta, out_file)
    write.table(grid, out_file, append=TRUE, row.names=FALSE, col.names=FALSE, na="-9999")
  }
}

# river_shp='~/IMERG/predata/dense_river.shp'

# flow_file <- '~/IMERG/predata/flow_bj.nc'
# domain_file <- '~/IMERG/predata/grids.nc'

# direc_grid=create_flow_direction(flow_file, domain_file, diag_tsh=0.3)

#' Plot the flow direction of the flow direction data.
#'
#' @description Plot the flow direction of a flow direction data (created by
#'              create_flow_direction()) for visullization and calibration. Can
#'              provide a shp file of the river network to plot helping the
#'              mannual calibration.
#'
#' @param direc Flow direction data created by create_flow_direction().
#' @param arc_code If use the ArcInfo direction code.
#' @param river_shp Path of the river shp file.
#'
#' @export
plot_flow_direc <- function(direc, arc_code=TRUE, river_shp=NULL, row_reverse=FALSE) {
  if(arc_code) {
    dir_code <- c(64, 128, 1, 2, 4, 8, 16, 32)
  } else {
    dir_code <- 1:8
  }

  if(!is.data.frame(direc)){
    csize <- direc$cellsize
    xcor <- direc$xllcorner
    ycor <- direc$yllcorner
    if(row_reverse) {
      direc$grid <- direc$grid[nrow(direc$grid):1, ]
    }
    direc <- grid2points(direc$grid, xcor=xcor, ycor=ycor, csize=csize)
  }

  xticks <- sort(unique(direc[ , 1]))
  yticks <- sort(unique(direc[ , 2]))
  xlabs <- round((xticks - min(xticks))/(xticks[2]-xticks[1]) + 1)
  ylabs <- round((yticks - min(yticks))/(yticks[2]-yticks[1]) + 1)

  plot(direc[ , 1:2], cex=0, xaxt='n', yaxt='n', xlab=NA, ylab=NA, asp=1)
  axis(1, xticks, xlabs)
  axis(3, xticks, xlabs)
  axis(2, yticks, ylabs)
  axis(4, yticks, ylabs)
  abline(h=yticks, col='grey')
  abline(v=xticks, col='grey')

  if(!is.null(river_shp)) {
    rs <- readShapeLines(river_shp)
    lines(rs, col='blue')
  }
  points(direc[ , 1:2], pch=20, cex=0.8)

  xcz <- sort(unique(direc[ ,1]))
  xcz <- mean(xcz[2:length(xcz)] - xcz[1:(length(xcz))-1])
  ycz <- sort(unique(direc[ ,2]))
  ycz <- mean(ycz[2:length(ycz)] - ycz[1:(length(ycz))-1])

  for(p in 1:nrow(direc)) {
    dx <- 0
    dy <- 0
    x <- direc[p, 1]
    y <- direc[p, 2]
    d <- direc[p, 3]
    if(d == dir_code[1]) {
      dy <- ycz
    } else if(d == dir_code[2]) {
      dy <- ycz
      dx <- xcz
    } else if(d == dir_code[3]) {
      dx <- xcz
    } else if(d == dir_code[4]) {
      dy <- -ycz
      dx <- xcz
    } else if(d == dir_code[5]) {
      dy <- -ycz
    } else if(d == dir_code[6]) {
      dy <- -ycz
      dx <- -xcz
    } else if(d == dir_code[7]) {
      dx <- -xcz
    } else if(d == dir_code[8]) {
      dx <- -xcz
      dy <- ycz
    }
    # arrows(c(x, x+dx), c(y, y+dy))
    Arrows(x, y, x+dx, y+dy, arr.length = 0.15, arr.adj=0.15)
  }
}

#' Revise the flow direction of the flow direction data.
#'
#' @description Revise the flow direction of the flow direction data created by
#'              create_flow_direction().
#'
#' @param direc_grid Flow direction data created by create_flow_direction().
#' @param row Row of the gridcell to be revised.
#' @param col Column of the gridcell to be revised.
#' @param direc New flow direction of the gridcell. For more convinient calibration,
#'              the direction code is sat as the keypad. It means that 1 is southwest,
#'              2 is south, 3 is southeast, and so on. set keypad_dir to FALSE can
#'              revise the direction code directly.
#' @param arc_code If use the ArcInfo flow direction code.
#' @param keypad_dir If use keypad flow direction code.
#'
#' @return The flow direction data after revise.
#'
#' @export
set_direc <- function(direc_grid, row, col, direc, arc_code=TRUE, keypad_dir=TRUE) {
  if(!is.list(direc_grid) | is.null(direc_grid$nrow | is.null(direc_grid$grid)))
    stop("direc_grid is not a ArcInfo grid.")
  if(arc_code) {
    dir_map <- c(8, 4, 2, 16, NA, 1, 32, 64, 128)
  } else {
    dir_map <- c(6, 5, 4, 7, NA, 3, 8, 1, 2)
  }
  nrow <- direc_grid$nrows
  if(keypad_dir) di <- dir_map[direc] else di <- direc
  direc_grid$grid[nrow-row+1, col] <- di
  return(direc_grid)
}

#' Write ArcInfo ASCII grid data.
#'
#' @description Write the grid data to the file as a ArcInfo type ASCII file by
#'              offering a list of grid data like the output of points2grid()
#'              or create_flow_direction(), or a matrix and some necessary
#'              informations like cell size, x and y corner of the grid data and
#'              so on.
#'
#' @param grid A list of grid data or a matrix. When it's a matrix, xcor, ycor and
#'             csize are necessary.
#' @param out_file Output file path.
#' @param xcor X corner of the grid. Necessary when grid is a matrix.
#' @param ycor Y corner of the grid. Necessary when grid is a matrix.
#' @param csize Cell size of the grid. Necessary when grid is a matrix.
#' @param NA_value A value to reprent the NA value.
#'
#' @export
write_arc_grid <- function(grid, out_file, xcor=NULL, ycor=NULL, csize=NULL, NA_value=-9999) {
  if(is.list(grid)) {
    if(is.null(grid$xllcorner) | is.null(grid$yllcorner) | is.null(grid$cellsize)) {
      stop("grid should have xllcorner, yllconer and cellsize.")
    }
    ncols=ncol(grid)
    nrows=nrow(grid)
    xllcorner <- grid$xllcorner
    yllcorner <- grid$yllcorner
    cellsize <- grid$cellsize
    griddata <- grid$grid
  } else {
    if(!is.matrix(grid)) {
      stop("grid should be a list or a matirx.")
    }
    if(is.null(xcor) | is.null(ycor) | is.null(csize)) {
      stop("xcor, ycor and csize must provided when grid is a matrix.")
    }
    ncols=ncol(grid)
    nrows=nrow(grid)
    xllcorner <- xcor
    yllcorner <- ycor
    cellsize <- csize
    griddata <- grid
  }

  griddata <- griddata[nrows:1, ]
  meta <- c(paste('ncols         ', round(ncols,  6)),
            paste('nrows         ', round(nrows,  6)),
            paste('xllcorner     ', round(xllcorner, 6)),
            paste('yllcorner     ', round(yllcorner, 6)),
            paste('cellsize      ', round(cellsize, 6)),
            paste('NODATA_value  ', NA_value))
  writeLines(meta, out_file)
  write.table(griddata, out_file, append=TRUE, row.names=FALSE, col.names=FALSE, na="-9999")
}

#' Read ArcInfo ASCII grid data.
#'
#' @description Read the grid data from ArcInfo type ASCII file.
#'
#' @param grid_file A ArcInfo ASCII grid file.
#' @return A grid data in list type, including nrows, ncols, xllcorner, yllcorner and cellsize.
#'
#' @export
read_arc_grid <- function(grid_file) {
  params <- strsplit(readLines(grid_file, 6), '\\s+')
  xcor <- as.numeric(params[[3]][2])
  ycor <- as.numeric(params[[4]][2])
  csize <- as.numeric(params[[5]][2])
  NA_value <- as.numeric(params[[6]][2])

  grid <- read.table(grid_file, skip=6)
  ncols <- ncol(grid)
  nrows <- nrow(grid)

  grid <- grid[nrows:1, ]
  grid[grid == NA_value] <- NA

  in_grid <- list(ncols=ncols,
                  nrows=nrows,
                  xllcorner=xcor,
                  yllcorner=ycor,
                  cellsize=csize,
                  grid=grid)

  return(in_grid)
}

#' Read the rows and columns of grids of a basin from a rout data file.
#'
#' @description Read the rows and columns of gridcells of a basin from a rout data file
#'              created by VIC Hime.
#'
#' @param data_path File path of the rou data file.
#' @return A data frame contains the rows and columns of the gridcells in the basin.
#'
#' @export
read_basin <- function(data_path) {
  basin <- fromJSON(readLines(data_path))$basin
  basin <- data.frame(t(data.frame(basin)))
  rownames(basin) <- 1:nrow(basin)
  colnames(basin) <- c('col', 'row')
  return(basin)
}

# direc_grid = set_direc(direc_grid, 19, 11, 1)
# direc_grid = set_direc(direc_grid, 14, 3, 3)
# direc_grid = set_direc(direc_grid, 17, 8, 9)
# plot_flow_direc(direc_grid, river_shp=river_shp)

#write_arc_grid(direc_grid,'~/IMERG/direc.txt')
