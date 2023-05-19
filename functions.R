# FUNCTIONS USED FOR BUILDING RISK MODEL

# function to crop rasters
get_crop <- function(CropLayer, Raster) {
  return(terra::crop(Raster, CropLayer, snap = "out"))
}

# function to generate zonal statistics
get_zonal <- function(ZonalLayer, Raster, Stat) {
  return(as_tibble(exact_extract(Raster, ZonalLayer, Stat)))
}

# function to create adjacency matrices
get_adjacency <- function(SpatialUnits, Name, FileLocation) {

  # get spatial neighbours
  nb <- poly2nb(SpatialUnits, snap = 2) # assume anything within 2m is a neighbour

  # create adjacency matrix
  nb2INLA(paste(FileLocation, Name, ".adj", sep = ""), nb)

  # return graph
  return(inla.read.graph(paste(FileLocation, Name, ".adj", sep = "")))
}

# function to format data for zero-inflated binomial model
get_zib_format <- function(Response, Ntrials, Covariates) {

  # set up response variables
  Resp1 <- as.vector(ifelse(Response > 0, 1, 0))
  Resp2 <- as.vector(Response[which(Response > 0),])
  Resp <- as_tibble(cbind(as.matrix(c(Resp1, rep(NA, length(Resp2)))), as.matrix(c(rep(NA, length(Resp1)), Resp2))))
  names(Resp) <- c("P", "N")

  # set up number of trials
  Trials <- as_tibble(rbind(as.matrix(rep(1, nrow(Response))), as.matrix(Ntrials[which(Response > 0),])))
  names(Trials) <- c("Ntrials")

  # set up intercepts
  Int1 <- rbind(as.matrix(rep(1, nrow(Response))), as.matrix(rep(0, length(which(Response > 0)))))
  Int2 <- rbind(as.matrix(rep(0, nrow(Response))), as.matrix(rep(1, length(which(Response > 0)))))
  Int <-as_tibble(cbind(Int1, Int2))
  names(Int) <- c("IntP", "IntN")

  # set up covariates
  Cov1 <- as_tibble(rbind(as.matrix(Covariates), matrix(rep(NA, length(which(Response > 0)) * ncol(Covariates)), nrow = length(which(Response > 0)), ncol = ncol(Covariates))))
  Cov2 <- as_tibble(rbind(matrix(rep(NA, nrow(Covariates) * ncol(Covariates)), nrow = nrow(Covariates), ncol = ncol(Covariates)), as.matrix(Covariates[which(Response > 0),])))
  names(Cov1) <- paste(names(Covariates), "P", sep = "")
  names(Cov2) <- paste(names(Covariates), "N", sep = "")
  Cov <- bind_cols(Cov1, Cov2)

  # return formatted data
  return(bind_cols(Resp, Trials, Int, Cov))
}