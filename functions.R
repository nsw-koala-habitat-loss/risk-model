# FUNCTIONS USED FOR BUILDING RISK MODEL

# function to crop rasters----
get_crop <- function(CropLayer, Raster) {
  return(terra::crop(Raster, CropLayer, snap = "out"))
}

# function to generate zonal statistics----
get_zonal <- function(ZonalLayer, Raster, Stat) {
  return(as_tibble(exact_extract(Raster, ZonalLayer, Stat)))
}

# function to generate zonal statistics allowing larger cells in memory----
get_zonal2 <- function(ZonalLayer, Raster, Stat) {
  return(as_tibble(exact_extract(Raster, ZonalLayer, Stat, max_cells_in_memory = 3.5e+08 )))
}

# function to create adjacency matrices----
get_adjacency <- function(SpatialUnits, Name, FileLocation) {
  
  # get spatial neighbours
  nb <- poly2nb(SpatialUnits, snap = 2) # assume anything within 2m is a neighbour
  
  # create adjacency matrix
  nb2INLA(paste(FileLocation, Name, ".adj", sep = ""), nb)
  
  # return graph
  return(inla.read.graph(paste(FileLocation, Name, ".adj", sep = "")))
}

# function to format data for zero-inflated binomial model----
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

# Function to fit model ----
fit_model <- function(KMR, ClearType, SpatUnits, RespData, CovsCD, SA1sPoly) {
 
  # KMR = text field for KMR - one of "CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC"
  # ClearType = clearing type - one of 1 = agriculture, 2 = infrastructure, 3 = forestry
  # SpatUnits = spatial units as list of Spatvector objects (properties) for each KMR
  # RespData = response data consisting listr of clearing and woody vegetation data for each KMR
  # CovsCD = covariates for each spatial unit as a list for each KMR
  # SA1Poly = SA1s spatial representation as a list of Spatvector objectd for each KMR
  
  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())
  
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
  ResultP <- inla(formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
  
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
  ResultN <- inla(formula, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
  
  # return models
  return(list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly))
}

# Wrapper function for INLA to limit run time----
## To prevent "idling". This should be use with IN:A_with_Retry Function
INLA_with_TimeLimit <- function(TimeLimit,...){
  setTimeLimit(cpu = TimeLimit, elapsed = TimeLimit, transient = TRUE)
  inla(...)
} 

# Wrapper to handle INLA errors and retry----
INLA_with_Retry <- function(N_retry=3, Initial_Tlimit = 1000,...){
  retry <- as.integer(0)
  result <- NULL
  while(retry<=N_retry){
    result <- tryCatch({
      INLA_with_TimeLimit(TimeLimit = (Initial_Tlimit + (retry*500)), ...)}, 
      error = function(e){
        cat("Error in INLA model fitting. Retrying...", retry+1, "\n", sep = "")
        Sys.sleep(1)
        NULL
      })
    if(!is.null(result)){break} 
    retry <- retry + 1
    if(retry > N_retry){
      stop("Error in INLA model fitting. Exceeded max number of retries \n")
      break 
    }
  }
  return(result)
}

# function to fit model with additional control for debugging----
fit_model2 <- function(KMR, ClearType, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD, SA1sPoly = SA1s, Explanatory = "All", Verbose = TRUE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL) {
  # ## Same as fit_model but with additional control for Explanatory variables and Verbose, allow for debugging
  # KMR = text field for KMR - one of "CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC"
  # ClearType = clearing type - one of 1 = agriculture, 2 = infrastructure, 3 = forestry
  # SpatUnits = spatial units as list of Spatvector objects (properties) for each KMR
  # RespData = response data consisting listr of clearing and woody vegetation data for each KMR
  # CovsCD = covariates for each spatial unit as a list for each KMR
  # SA1Poly = SA1s spatial representation as a list of Spatvector objectd for each KMR
  # Explanatory = covariates to include in model - either "All" or a vector of covariate names
  # Verbose = logical to print progress
  # N_retry = number of retries for INLA model fitting (passed to INLA_with_TimeLimit)
  # Initial_Tlimit = initial time limit for INLA model fitting (passed to INLA_with_Retry)
  # OutputDir = directory to save model output (optional)
  
  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())
  
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
  ExplV <- if(Explanatory == "All") {paste(names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID)), collapse=" + ")} else {paste(Explanatory)}
  
  # get adjacency matrix for SA1s containing properties with forest cover
  SA1IDs <- CP %>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  SA1sPolyKMR <- SA1sPoly[[KMR]] %>% left_join(SA1IDs, join_by(SA1_CODE21 == SA1), keep = TRUE) %>% filter(!is.na(SA1ID)) %>% arrange(SA1ID)
  Adj <- SA1sPolyKMR %>% get_adjacency(paste0("modelP_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")
  
  # fit clearing versus no clearing model
  formula <- as.formula(paste0(paste("P", ExplV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)

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
  formula <- as.formula(paste0(paste("Prop", ExplV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = Adj, scale.model = TRUE)"))
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  
  # return models
  if(!is.null(OutputDir)){
    Model <- list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly)
    output_FPath <- file.path(OutputDir, paste0("Model_", KMR, "_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, ".qs"))
    qsave(Model, file = output_FPath)
  }
  return(list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly))
}

# function to create spatial predictions from model----
predict_model <- function(Model) {
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
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())
  
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
  ResultP <- inla(formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
  
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
  ResultN <- inla(formula, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = TRUE)
  
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


# function for model selection----
Select_model <- function(KMR = "CC", ClearType = 1, SpatUnits = SUs_Ag, RespData = ZStats_Woody_Ag, CovsCD = ZStats_Covs_Ag, SA1sPoly = SA1s, Direction = "Forward Complete", Verbose = FALSE, N_retry=3, Initial_Tlimit = 1000, OutputDir = NULL) {
  # KMR = text field for KMR - one of "CC", "CST", "DRP", "FW", "NC", "NT", "NS", "R", "SC"
  # ClearType = clearing type - one of 1 = agriculture, 2 = infrastructure, 3 = forestry
  # SpatUnits = spatial units as list of Spatvector objects (properties) for each KMR
  # RespData = response data consisting listr of clearing and woody vegetation data for each KMR
  # CovsCD = covariates for each spatial unit as a list for each KMR
  # SA1Poly = SA1s spatial representation as a list of Spatvector objectd for each KMR
  # Direction = direction of model selection - one of "Forward", "F", "Backward", "B", "Forward Complete", "FC", "Backward Complete", "BC"
  # Verbose = logical to print progress
  # N_retry = number of retries for INLA model fitting (passed to INLA_with_TimeLimit)
  # Initial_Tlimit = initial time limit for INLA model fitting (passed to INLA_with_Retry)
  # OutputDir = directory to save model output (optional)
  
  tic("Total Time")
  
  # Filter input data to KMR only (save memory)
  SpatUnits <- SpatUnits[[KMR]] %>%   list() %>% setNames(KMR)
  RespData <- RespData[[KMR]] %>% list() %>% setNames(KMR)
  CovsCD <- CovsCD[[KMR]] %>% list() %>% setNames(KMR)
  SA1sPoly <- SA1sPoly[[KMR]] %>% list() %>% setNames(KMR)
  
  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())
  
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
  
  
  ## Format INLA Parameters ----
  ### Probability of clearing model ----
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
  AdjP <- SA1sPolyKMR %>% get_adjacency(paste0("modelP_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")
  
  ### Proportion cleared|clearing model ----
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
  AdjN <- SA1sPolyKMR %>% get_adjacency(paste0("modelN_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, KMR, "_Adj_SA1s"), "output/neighbours/")
  
  
  if(Direction == "Forward" | Direction == "F"){
    
    #'######################################
    ## Forward stepwise model selection ----
    #'######################################  
    tic("Forward Model Selection")
    
    # Get a list of Full explanatory variables and explanatory variables yet to be tested
    left_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
    Sel_explV_ls <- list()
    Best_DIC_ls <- list() ##**##
    DIC_ls <- dDIC_ls <-  list() 
    
    ### Null Model----
    tic("Null Model")
    # fit clearing versus no clearing model 
    ForP_H0 <- as.formula(paste0(paste("P", 1, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
    ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H0, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # fit proportion cleared|clearing model
    ForN_H0 <- as.formula(paste0(paste("Prop", 1, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
    ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H0, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # Calculate DIC for the null model
    Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
    Best_DIC <- Current_DIC
    Current_dDIC <- Current_DIC + (2*0)
    Best_dDIC <- Current_dDIC
    
    Best_DIC_ls[1] <- Best_DIC ##**##
    names(Best_DIC_ls)[1] <- "H0" ##**##
    DIC_ls[1] <- Current_DIC 
    dDIC_ls[1] <- Current_dDIC
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- "H0" 
    
    toc(log = TRUE) # null model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", "1", "\n", sep = "")
    
    while(length(left_explV_ls) > 0){
      #### Step forward... ----
      tic(paste0("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      
      # Reset best explanatory variable
      Best_explV <- NULL
      
      for(x in 1: length(left_explV_ls)){
        
        # Update the list of selected explanatory variables
        explV <- paste(c(Sel_explV_ls, left_explV_ls[x]), collapse=" + ")
        
        # fit clearing versus no clearing model 
        ForP_H1 <- as.formula(paste0(paste("P", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
        tic(paste0("ExplV ", x, " of ", length(left_explV_ls)))
        # cat("Running: ", deparse(ForP_H1), "\n", sep = "")
        ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        
        # fit proportion cleared|clearing model
        ForN_H1 <- as.formula(paste0(paste("Prop", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
        # cat(deparse(ForN_H1), "\n", sep = "")
        ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        toc(log = TRUE) # current model
        
        # Calculate DIC for the current model and compare it with the best model
        Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
        cat("DIC = ", format(round(Current_DIC, 2)), "   ExplV: ", explV, "\n", sep = "")
        DIC_ls[length(DIC_ls)+1] <- Current_DIC 
        Current_dDIC <- Current_DIC + (2*length(explV))
        dDIC_ls[length(dDIC_ls)+1] <- Current_dDIC
        names(DIC_ls)[length(DIC_ls)] <- names(dDIC_ls)[length(dDIC_ls)] <- paste(explV) 
        
        # Update the best DIC if the DIC of the current model is better
        if(Current_dDIC < Best_dDIC){
          Best_DIC <- Current_DIC
          Best_dDIC <- Current_dDIC
          Best_explV <- left_explV_ls[x]
        }
      }
      
      # Update the list of selected explanatory variables if there is a new Best explanatory variable & the DIC difference is greater than 2
      if(!is.null(Best_explV)){
        
        Sel_explV_ls <- c(Sel_explV_ls, Best_explV)
        left_explV_ls <- setdiff(Full_explV_ls, Sel_explV_ls)
        
        Best_DIC_ls[length(Best_DIC_ls)+1] <- Best_DIC ##**##
        names(Best_DIC_ls)[length(Best_DIC_ls)] <- paste(unlist(Sel_explV_ls), collapse=" + ") ##**##
        # cat("\014")    
        toc(log = TRUE) # Step
        cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "")
        cat("Explanatory variable not tested yet: ", paste(left_explV_ls, collapse = ", "), "\n\n", sep = "")
        
      }else if (length(DIC_ls)==1){
        Sel_explV_ls <- 1
        warning("Null model selected!")
        break
      } else {
        break
      }
      
    }
    toc(log = TRUE) # Forward Model Selection
  }
  
  else if(Direction == "Backward" | Direction  == "B"){
    
    #'##############################
    ## Backward model selection ----
    #'##############################
    
    tic("Backward Model Selection")
    
    # Get a list of Full explanatory variables and start by selecting all
    Sel_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
    Best_DIC_ls <- list() ##**##
    DIC_ls <- dDIC_ls <-  list() 
    worst_explV_ls <- list()
    
    ### Full Model ----
    tic("Full Model")
    # fit clearing versus no clearing model 
    ForP_H0 <- as.formula(paste0(paste("P", paste(Full_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
    ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H0, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # fit proportion cleared|clearing model
    ForN_H0 <- as.formula(paste0(paste("Prop", paste(Full_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
    ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H0, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # Calculate DIC for the full model
    Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
    Best_DIC <- Current_DIC
    Current_dDIC <- Current_DIC + (2*length(Full_explV_ls))
    Best_dDIC <- Current_dDIC
    
    Best_DIC_ls[1] <- Best_DIC ##**##
    names(Best_DIC_ls)[1] <- paste(unlist(Full_explV_ls), collapse=" + ") ##**##
    DIC_ls[1] <- Current_DIC
    dDIC_ls[1] <- Current_DIC + (2*length(Full_explV_ls))
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- paste(unlist(Full_explV_ls), collapse=" + ")
    
    toc(log = TRUE) # Full Model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", paste(Full_explV_ls, collapse = ", "), "\n", sep = "")
    
    
    while(length(Sel_explV_ls) > 0 ){
      #### Step backward... ----
      tic(paste0("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      Worst_explV <- NULL
      
      for(x in 1: length(Sel_explV_ls)){
        
        explV <- paste(Sel_explV_ls[-x], collapse=" + ")
        
        # fit clearing versus no clearing model 
        ForP_H1 <- as.formula(paste0(paste("P", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
        tic(paste0("ExplV ", x, " of ", length(Sel_explV_ls)))
        # cat("Running: ", deparse(ForP_H1), "\n", sep = "")
        ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        
        # fit proportion cleared|clearing model
        ForN_H1 <- as.formula(paste0(paste("Prop", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
        # cat(deparse(ForN_H1), "\n", sep = "")
        ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        toc(log = TRUE) # Current Model
        
        # Calculate DIC for the current model and compare it with the best model
        Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
        cat("DIC = ", format(round(Current_DIC, 2)), "   ExplV: ", explV, "\n", sep = "")
        DIC_ls[length(DIC_ls)+1] <- Current_DIC 
        Current_dDIC <- Current_DIC+ if(explV!=""){2*length(unlist(strsplit(explV, "\\+")))} else {0}
        dDIC_ls[length(dDIC_ls)+1] <- Current_dDIC
        names(DIC_ls)[length(DIC_ls)] <- names(dDIC_ls)[length(dDIC_ls)] <- if(explV!="") {paste(explV)} else {"H0"} 
        
        # Update the best DIC if the DIC of the current model is better
        if(Current_dDIC < Best_dDIC){
          Best_DIC <- Current_DIC
          Best_dDIC <- Current_dDIC
          Worst_explV <- Sel_explV_ls[x]
        }
      }
      
      # Update the list of selected explanatory variables if there is a new Best explanatory variable & the DIC difference is greater than 2
      if(!is.null(Worst_explV)){
        Sel_explV_ls <- setdiff(Sel_explV_ls, Worst_explV)
        worst_explV_ls[length(worst_explV_ls)+1] <- Worst_explV
        
        Best_DIC_ls[length(Best_DIC_ls)+1] <-  Best_DIC ##**##
        names(Best_DIC_ls)[length(Best_DIC_ls)] <- paste(unlist(Sel_explV_ls), collapse=" + ")
        # cat("\014")    
        toc(log = TRUE) # Step
        cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "")
        cat("Explanatory variable removed: ", paste(worst_explV_ls, collapse = ", "), "\n\n", sep = "")
        
      }else if (length(Sel_explV_ls )==0){
        Sel_explV_ls <- 1
        warning("Null model selected!")
        break
      } else {
        break
      }
    }
    toc(log = TRUE) # Backward Model Selection
  } 
  
  else if (Direction == "Forward complete" | Direction == "FC"){
    #'######################################
    ## Forward complete model selection ----
    #'######################################
    
    tic("Forward Complete Model Selection")
    # Get a list of Full explanatory variables and explanatory variables yet to be tested
    left_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
    Sel_explV_ls <- list()
    Best_DIC_ls <- list() ##**##
    DIC_ls <- dDIC_ls <-  list() 
    
    ## Null Model ----
    tic("Null Model")
    # fit clearing versus no clearing model 
    ForP_H0 <- as.formula(paste0(paste("P", 1, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
    ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H0, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # fit proportion cleared|clearing model
    ForN_H0 <- as.formula(paste0(paste("Prop", 1, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
    ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H0, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # Calculate DIC for the null model
    Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
    Best_DIC <- ResultP$dic$dic + ResultN$dic$dic ##**##
    Best_DIC_ls[1] <- Best_DIC ##**##
    names(Best_DIC_ls)[1] <- "H0" ##**##
    DIC_ls[1] <- Current_DIC 
    dDIC_ls[1] <- Current_DIC + (2*0)
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- "H0" 
    
    toc(log = TRUE) # null model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", "1", "\n", sep = "")
    
    while(length(left_explV_ls) > 0){
      #### Step forward... ----
      tic(paste0("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      Best_explV <- NULL 
      Best_DIC <- Inf 
      
      for(x in 1: length(left_explV_ls)){
        
        explV <- paste(c(Sel_explV_ls, left_explV_ls[x]), collapse=" + ")
        
        # fit clearing versus no clearing model 
        ForP_H1 <- as.formula(paste0(paste("P", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
        tic(paste0("ExplV ", x, " of ", length(left_explV_ls)))
        # cat("Running: ", deparse(ForP_H1), "\n", sep = "")
        ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        
        # fit proportion cleared|clearing model
        ForN_H1 <- as.formula(paste0(paste("Prop", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
        # cat(deparse(ForN_H1), "\n", sep = "")
        ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        toc(log = TRUE) # current model
        
        # Calculate DIC for the current model and compare it with the best model
        Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
        cat("DIC = ", format(round(Current_DIC, 2)), "   ExplV: ", explV, "\n", sep = "")
        DIC_ls[length(DIC_ls)+1] <- Current_DIC 
        dDIC_ls[length(dDIC_ls)+1] <- Current_DIC + 2*length(unlist(strsplit(explV, "\\+")))
        names(DIC_ls)[length(DIC_ls)] <- names(dDIC_ls)[length(dDIC_ls)] <- paste(explV) 
        
        # Update the best DIC if the DIC of the current model is better
        if(Current_DIC < Best_DIC){
          Best_DIC <- Current_DIC
          Best_explV <- left_explV_ls[x]
        }
      }
      
      Sel_explV_ls <- c(Sel_explV_ls, Best_explV)
      left_explV_ls <- setdiff(Full_explV_ls, Sel_explV_ls)
      
      Best_DIC_ls[length(Best_DIC_ls)+1] <- Best_DIC
      names(Best_DIC_ls)[length(Best_DIC_ls)] <- paste(unlist(Sel_explV_ls), collapse=" + ")
      # cat("\014")    
      toc(log = TRUE) # Step
      cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "") 
      
      cat("Explanatory variable not tested yet: ", paste(left_explV_ls, collapse = ", "), "\n\n", sep = "")
    }
    
    Best_explvs <- names(dDIC_ls)[which.min(dDIC_ls)]
    if(Best_explvs != "H0"){
      Sel_explV_ls <- unlist(strsplit(Best_explvs, " \\+ "))
    }else{
      Sel_explV_ls <- 1
      warning("Null model selected!")
      }
    toc(log = TRUE) # Forward Complete Model Selection
    FC_DIC_ls <- DIC_ls ##**##
    FC_dDIC_ls <- dDIC_ls ##**##
  }
  
  else if (Direction == "Backward Complete" | Direction == "BC"){
    #'######################################
    ## Backward complete model selection ----
    #'######################################
    
    tic("Backward complete Model Selection")
    
    # Get a list of Full explanatory variables and start by selecting all
    Sel_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
    Best_DIC_ls <- list() ##**##
    DIC_ls <- dDIC_ls <-  list() 
    worst_explV_ls <- list()
    
    ### Full Model ----
    tic("Full Model")
    # fit clearing versus no clearing model 
    ForP_H0 <- as.formula(paste0(paste("P", paste(Full_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
    ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H0, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # fit proportion cleared|clearing model
    ForN_H0 <- as.formula(paste0(paste("Prop", paste(Full_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
    ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H0, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
    
    # Calculate DIC for the full model
    Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
    Best_DIC <- ResultP$dic$dic + ResultN$dic$dic ##**##
    Best_DIC_ls[1] <- Best_DIC ##**##
    names(Best_DIC_ls)[1] <- paste(unlist(Full_explV_ls), collapse=" + ") ##**##
    DIC_ls[1] <- Current_DIC
    dDIC_ls[1] <- Current_DIC + (2*length(Full_explV_ls))
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- paste(unlist(Full_explV_ls), collapse=" + ")
    
    toc(log = TRUE) # Full Model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", paste(Full_explV_ls, collapse = ", "), "\n", sep = "")
    
    while(length(Sel_explV_ls) > 0 ){
      #### Step backward... ----
      tic(paste0("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      Worst_explV <- NULL
      Best_DIC <- Inf
      
      for(x in 1: length(Sel_explV_ls)){
        
        explV <- paste(Sel_explV_ls[-x], collapse=" + ")
        
        # fit clearing versus no clearing model 
        ForP_H1 <- as.formula(paste0(paste("P", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
        tic(paste0("ExplV ", x, " of ", length(Sel_explV_ls)))
        # cat("Running: ", deparse(ForP_H1), "\n", sep = "")
        ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        
        # fit proportion cleared|clearing model
        ForN_H1 <- as.formula(paste0(paste("Prop", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
        # cat(deparse(ForN_H1), "\n", sep = "")
        ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        toc(log = TRUE) # Current Model
        
        # Calculate DIC for the current model and compare it with the best model
        Current_DIC <- ResultP$dic$dic + ResultN$dic$dic
        cat("DIC = ", format(round(Current_DIC, 2)), "   ExplV: ", explV, "\n", sep = "")
        DIC_ls[length(DIC_ls)+1] <- Current_DIC 
        dDIC_ls[length(dDIC_ls)+1] <- Current_DIC + if(explV!=""){2*length(unlist(strsplit(explV, "\\+")))} else {0}
        names(DIC_ls)[length(DIC_ls)] <- names(dDIC_ls)[length(dDIC_ls)] <- if(explV!="") {paste(explV)} else {"H0"} 
        
        # Update the best DIC if the DIC of the current model is better
        if(Current_DIC < Best_DIC){
          Best_DIC <- Current_DIC
          Worst_explV <- Sel_explV_ls[x]
        }
      }
      
      Sel_explV_ls <- setdiff(Sel_explV_ls, Worst_explV)
      worst_explV_ls[length(worst_explV_ls)+1] <- Worst_explV
      
      Best_DIC_ls[length(Best_DIC_ls)+1] <-  Best_DIC ##**##
      names(Best_DIC_ls)[length(Best_DIC_ls)] <- paste(unlist(Sel_explV_ls), collapse=" + ") ##**##
      # cat("\014")    
      toc(log = TRUE) # Step
      cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "")
      cat("Explanatory variable removed: ", paste(worst_explV_ls, collapse = ", "), "\n\n", sep = "")
    
    }
    
    Best_explvs <- names(dDIC_ls)[which.min(dDIC_ls)]
    
    if(Best_explvs != "H0"){
      Sel_explV_ls <- unlist(strsplit(Best_explvs, " \\+ "))
    }else{
      Sel_explV_ls <- 1
      warning("Null model selected!")
    }
    toc(log = TRUE) # Forward Complete Model Selection
  } else{
    stop("Invalid model selection method\n Please select one of 'Forward', 'F', 'Backward', 'B', 'Complete Forward', 'CF' or 'Complete Backward', 'CB'")
  }

  
  ## Fit final model with selected explanatory variables ----
  # fit clearing versus no clearing model 
  ForP_H1 <- as.formula(paste0(paste("P", paste(Sel_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
  cat("Running Best model: ", deparse(ForP_H1), "\n")
  ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  
  # fit proportion cleared|clearing model
  ForN_H1 <- as.formula(paste0(paste("Prop", paste(Sel_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
  cat(deparse(ForN_H1), "\n")
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  toc(log = TRUE) # Total Time
  
  # filter CovCD to only include selected explanatory variables (For Predict function)
  CovsCD <- CovsCD[[KMR]] %>% dplyr::select(all_of(Sel_explV_ls)) %>%  list() %>% setNames(KMR)
  
  # return models
  if(!is.null(OutputDir)){
    Model <- list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly, DIC_ls = DIC_ls)
    output_FPath <- file.path(OutputDir, paste0("SelModel_", KMR, "_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, Direction, ".qs"))
    qsave(Model, file = output_FPath)
  }
  return(list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly, DIC_ls = DIC_ls))
}