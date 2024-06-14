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

fit_model <- function(KMR, ClearType, SpatUnits, RespData, CovsCD, SA1sPoly) {
# function to fit model
# KMR = text field for KMR - one of "CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC"
# ClearType = clearing type - one of 1 = agriculture, 2 = infrastructure, 3 = forestry
# SpatUnits = spatial units as list of Spatvector objects (properties) for each KMR
# RespData = response data consisting listr of clearing and woody vegetation data for each KMR
# CovsCD = covariates for each spatial unit as a list for each KMR
# SA1Poly = SA1s spatial representation as a list of Spatvector objectd for each KMR

  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area, -Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())

  # get response data and ensure clearing is < woody vegetation
  Response <- RespData[[KMR]] %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>%
             mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody)

  # set up data for INLA models

  # response - how many cells cleared over time period
  if (ClearType == 1) {
    R <- Response %>% dplyr::select(YAg) %>% as.matrix()
  } else if (ClearType == 2) {
    R <- Response %>% dplyr::select(YIn) %>% as.matrix()
  } else if (ClearType == 3) {
    R <- Response %>% dplyr::select(YFo) %>% as.matrix()
  }

  # number of trials - how many cells woody at start of time period
  NT <- Response %>% dplyr::select(N) %>% as.matrix()

  # probability of clearing model

  # format data for model fitting
  RP <- R[which(NT > 0)] # only fit to data for properties with woody vegetation in 2011
  RP <- as.vector(ifelse(RP > 0, 1, 0)) # recode to binary cleared/not cleared
  NTP <- rep(1, length(RP)) # set number of trials to 1 for all properties in 2011
  ResponseP <- as_tibble(cbind(RP, NTP))
  names(ResponseP) <- c("P", "Ntrials")
  CP <- Covs[which(NT > 0), ] # only fit to data for properties with woody vegetation
  CP <- CP %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
  DataP <- bind_cols(ResponseP, CP)

  # get adjacency matrix for SA1s containing properties with forest cover
  SA1IDs <- CP %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  SA1sPolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
  Adj <- SA1sPolyKMR %>% get_adjacency(paste0("modelP_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")

  # fit clearing versus no clearing model
  formula <- as.formula(paste0(paste("P", paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultP <- inla(formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)

  # proportion cleared|clearing model

  # format data for model fitting
  RN <- R[which(R > 0)] # only fit to data from properties with some clearing
  NTN <- NT[which(R > 0)] # only fit to data from properties with some clearing
  ResponseN <- as_tibble(cbind(as.matrix(RN), as.matrix(NTN))) %>% mutate(Prop = RN / NTN) %>% mutate(Prop = ifelse(Prop < 1, Prop, 0.999)) # calculate the proportion cleared (when proportion = 1 set to 0.999)
  names(ResponseN) <- c("N", "Ntrials", "Prop")
  CN <- Covs[which(R > 0), ] # only fit to data from properties with some clearing
  CN <- CN %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
  DataN <- bind_cols(ResponseN, CN)

  # get adjacency matrix for SA1s containing properties with some clearing
  SA1IDs <- CN %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  SA1PolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
  Adj <- SA1sPolyKMR %>% get_adjacency(paste0("modelN_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")

  # fit proportion cleared|clearing model
  formula <- as.formula(paste0(paste("Prop", paste(names(CN %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultN <- inla(formula, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)

  # return models
  return(list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly))
}

predict_model <- function(Model) {
# function to create spatial predictions from model
# Model = the fitted model (an output from fit_model)

  # define objects
  PModel = Model$PModel
  NModel = Model$NModel
  KMR = Model$KMR
  ClearType = Model$ClearType
  SpatUnits = Model$SpatUnits
  RespData = Model$RespData
  CovsCD = Model$CovsCD
  SA1sPoly = Model$SA1sPoly
  rm(Model)

  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area, -Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())

  # get response data and ensure clearing is < woody vegetation
  Response <- RespData[[KMR]] %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>% mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody)

  # set up data for INLA models

  # response - how many cells cleared over time period
  if (ClearType == 1) {
    R <- Response %>% dplyr::select(YAg) %>% as.matrix()
  } else if (ClearType == 2) {
    R <- Response %>% dplyr::select(YIn) %>% as.matrix()
  } else if (ClearType == 3) {
    R <- Response %>% dplyr::select(YFo) %>% as.matrix()
  }

  # number of trials - how many cells woody at start of time period
  NT <- Response %>% dplyr::select(N) %>% as.matrix()

  # probability of clearing model

  # format data for model fitting
  RP <- R[which(NT > 0)] # only fit to data for properties with woody vegetation in 2011
  RP <- as.vector(ifelse(RP > 0, 1, 0)) # recode to binary cleared/not cleared
  NTP <- rep(1, length(RP)) # set number of trials to 1 for all properties in 2011
  ResponseP <- as_tibble(cbind(RP, NTP))
  names(ResponseP) <- c("P", "Ntrials")
  CP <- Covs[which(NT > 0), ] # only fit to data for properties with woody vegetation
  CP <- CP %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
  DataP <- bind_cols(ResponseP, CP)

  # get adjacency matrix for SA1s containing properties with forest cover
  SA1IDs <- CP %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  SA1sPolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
  Adj <- SA1sPolyKMR %>% get_adjacency(paste0("modelP_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")

  # format data for model predictions
  # only include properties that have woody vegetation in 2011
  RPPred <- R[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
  RPPred <- as.vector(ifelse(RPPred > 0, 1, 0)) # recode to binary cleared/not cleared
  NTPPred <- rep(1, length(RPPred)) # set number of trials to 1 for all properties
  ResponsePPred <- as_tibble(cbind(RPPred, NTPPred))
  names(ResponsePPred) <- c("P", "Ntrials")
  CPPred <- Covs %>% left_join(CP %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
  CPPred <- CPPred[which(NT > 0), ] # only make predictions for properties with woody vegetation in 2011
  DataPPred <- bind_cols(ResponsePPred, CPPred)
  DataPPred <- DataPPred %>% mutate(P = NA) # set whether cleared or not to NA so as to make predictions

  # combine fitting and prediction data
  DataP <- bind_rows(DataPPred, DataP)

  # fit clearing versus no clearing model
  formula <- as.formula(paste0(paste("P", paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultP <- inla(formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)

  # proportion cleared|clearing model

  # format data for model fitting
  RN <- R[which(R > 0)] # only fit to data from properties with some clearing
  NTN <- NT[which(R > 0)] # only fit to data from properties with some clearing
  ResponseN <- as_tibble(cbind(as.matrix(RN), as.matrix(NTN))) %>% mutate(Prop = RN / NTN) %>% mutate(Prop = ifelse(Prop < 1, Prop, 0.999)) # calculate the proportion cleared (when proportion = 1 set to 0.999)
  names(ResponseN) <- c("N", "Ntrials", "Prop")
  CN <- Covs[which(R > 0), ] # only fit to data from properties with some clearing
  CN <- CN %>% mutate(SA1ID = as.integer(factor(SA1))) # recode indices for SA1s for random-effect
  DataN <- bind_cols(ResponseN, CN)

  # get adjacency matrix for SA1s containing properties with some clearing
  SA1IDs <- CN %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  SA1PolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
  Adj <- SA1sPolyKMR %>% get_adjacency(paste0("modelN_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")

  # format data for model predictions
  # this time we only include properties that have woody vegetation in 2011
  RNPred <- R[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
  NTNPred <- NT[which(NT > 0)] # only make predictions for properties with woody vegetation in 2011
  ResponseNPred <- as_tibble(cbind(as.matrix(RNPred), as.matrix(NTNPred))) %>% mutate(Prop = NA) # set Prop to NA so we get predictions for these properties
  names(ResponseNPred) <- c("N", "Ntrials", "Prop")
  CNPred <- Covs %>% left_join(CN %>% distinct(SA1, SA1ID), join_by(SA1 == SA1)) # recode indices for SA1s for random-effect (same IDs as for training data)
  CNPred <- CNPred[which(NT > 0), ] # only make predictions for properties with woody vegetation in 2011
  DataNPred <- bind_cols(ResponseNPred, CNPred)

  # combine fitting and prediction data
  DataN <- bind_rows(DataNPred, DataN)

  # fit proportion cleared|clearing model
  formula <- as.formula(paste0(paste("Prop", paste(names(CN %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultN <- inla(formula, data = DataN, family = "beta", control.inla = list(control.vb = list(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)

  # spatialise predictions
  Layer <- SpatUnits[[KMR]] %>% bind_cols(R = as.vector(R), NT = as.vector(NT), SUID = Covs$SUID) %>% dplyr::filter(NT > 0) %>% mutate(ActualProp =  ifelse(NT > 0, R / NT, NA)) # only for properties with woody cover
  PredictionsP <- as_tibble(ResultP$summary.fitted.values$mean[1:nrow(DataPPred)]) %>% bind_cols(SUID = as.vector(Covs$SUID[which(NT > 0)]))
  names(PredictionsP) <- c("PredP", "SUID")
  PredictionsN <- as_tibble(ResultN$summary.fitted.values$mean[1:nrow(DataNPred)]) %>% bind_cols(SUID = as.vector(Covs$SUID[which(NT > 0)]))
  names(PredictionsN) <- c("PredN", "SUID")
  PredictionsCombined <- PredictionsP %>% left_join(PredictionsN, by = join_by(SUID == SUID)) %>% mutate(PredAll = PredP * PredN)
  Layer <- Layer %>% left_join(PredictionsCombined, by = join_by(SUID == SUID)) %>% dplyr::select(-Shape_Length, -Shape_Area)

  return(list(Layer = Layer, PModel = PModel, NModel = NModel))
}