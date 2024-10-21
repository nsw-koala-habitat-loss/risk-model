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
      message("Error in INLA model fitting. Exceeded max number of retries \n")
      result <- list(dic = list(dic = NULL))
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
  ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE,config = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)

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
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE, config = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  
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

# function to create spatial predictions from model----
predict_model2 <- function(Model, Verbose = FALSE, N_retry=10, Initial_Tlimit = 1000) {
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
  ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, formula, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  
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

# function to create spatial predictions by sampling from posterior----
predict_model3 <- function(model, N = 1000, RandEff = "SA1ID", verbose = TRUE){
  # model = the fitted model (an output from fit_model)
  # N = number of samples to draw from the posterior
  # RandEff = the name of the random effect
  
  # define objects
  PModel = model$PModel
  NModel = model$NModel
  KMR = model$KMR
  ClearType = model$ClearType
  SpatUnits = model$SpatUnits
  RespData = model$RespData
  CovsCD = model$CovsCD
  SA1sPoly = model$SA1sPoly
  rm(model)
  
  # get attribute table of spatial units and join covariates
  Covs <- SpatUnits[[KMR]] %>% st_drop_geometry() %>% as_tibble() %>% dplyr::select(-KMR, -Shape_Length, -Shape_Area) %>% bind_cols(CovsCD[[KMR]]) %>% mutate(SUID = 1:n())
  
  # get response data and ensure clearing is < woody vegetation
  Response <- RespData[[KMR]] %>% mutate(YAg = sum.aloss, YIn = sum.iloss, YFo = sum.floss, N = sum.woody) %>% mutate(N = ifelse(N < YAg + YIn + YFo, YAg + YIn + YFo, N)) %>% dplyr::select(-sum.aloss, -sum.iloss, -sum.floss, -sum.woody , -sum.Khab)
  
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
  
  # format data for prediction
  CPred <- Covs[which(Response$N>0),] %>% dplyr::select(-SA1, -SUID) # Get Covariate data only for properties with woody vegetation
  CPred_MM <- model.matrix(~., data = CPred) # model matrix
  SA1IDs <- Covs[which(Response$N>0),] %>% mutate(SA1ID = as.integer(factor(SA1))) %>% dplyr::select(SA1ID, SUID)#%>% dplyr::select(SA1, SA1ID) %>% group_by(SA1) %>% summarise(SA1ID = first(SA1ID))
  
  # get posterior samples for each model
  if(verbose) {cat("\nSampling from posterior...\n")}
  PModel_samp <- inla.posterior.sample(N, PModel)
  NModel_samp <- inla.posterior.sample(N, NModel)
  if(verbose) {cat("inla.posterior.sample object size: \n PModel: ", format(object.size(PModel_samp), units = "auto"), 
                   "\n NModel: ",format(object.size(NModel_samp), units = "auto"))}
  
  # get fixed effects samples
  PModel_samp_fixed <- sapply(PModel_samp, function(x) x$latent[(nrow(x$latent)-length(PModel$names.fixed)+1):nrow(x$latent),1]) %>% as.matrix()
  NModel_samp_fixed <- sapply(NModel_samp, function(x) x$latent[(nrow(x$latent)-length(NModel$names.fixed)+1):nrow(x$latent),1]) %>% as.matrix()
  if(verbose) {cat("\n\nFixed effects samples: \n PModel: ", format(object.size(PModel_samp_fixed), units = "auto"), 
                   "\n NModel: ",format(object.size(NModel_samp_fixed), units = "auto"))}
  
  # get random effects samples
  ### Ref: A tutorial in spatial and spatio-temporal models with R-INLA (https://discovery.ucl.ac.uk/id/eprint/1415919/1/Baio_BlaCamBaiRue.pdf) Page 11
  if(verbose) {cat("\n\nExtracting random effects samples...\n")}
  RandEff_PId <- which(PModel$misc$configs$contents$tag == RandEff)
  RandEff_PInd <- PModel$misc$configs$contents$start[RandEff_PId] - 1 + (1:(PModel$misc$configs$contents$length[RandEff_PId]))
  PModel_samp_rand_mt <- sapply(PModel_samp, function(x){x$latent[RandEff_PInd[1:(length(RandEff_PInd)/2)],1]}) %>% as_tibble(rownames = "RName") %>% mutate(SA1ID = as.integer(str_extract(RName, "(?<=:)[0-9]+"))) %>% dplyr::select(-RName) %>% right_join(x= ., y=SA1IDs, by = "SA1ID") %>% arrange(SUID) %>% dplyr::select(-SA1ID, -SUID) %>% as.matrix()
  
  RandEff_NId <- which(NModel$misc$configs$contents$tag == RandEff)
  RandEff_NInd <- NModel$misc$configs$contents$start[RandEff_NId] - 1 + (1:(NModel$misc$configs$contents$length[RandEff_NId]))
  NModel_samp_rand_mt <- sapply(NModel_samp, function(x){x$latent[RandEff_NInd[1:(length(RandEff_NInd)/2)],1]}) %>% as_tibble(rownames = "RName") %>% mutate(SA1ID = as.integer(str_extract(RName, "(?<=:)[0-9]+"))) %>% dplyr::select(-RName) %>% right_join(x= ., y=SA1IDs, by = "SA1ID") %>% arrange(SUID) %>% dplyr::select(-SA1ID, -SUID) %>% as.matrix()
  rm(PModel_samp, NModel_samp)
  if(verbose) {cat("Random effects samples: \n PModel: ", format(object.size(PModel_samp_rand_mt), units = "auto"), 
                   "\n NModel: ",format(object.size(NModel_samp_rand_mt), units = "auto"))}
  
  # Prediction for fixed effects
  if(verbose) {cat("\n\nMatrix multiplication for fixed effects...\n")}
  PredP_samp_fixed_mt <- apply(PModel_samp_fixed , MARGIN = 2, function(Model_samp) {CPred_MM %*% Model_samp }) 
  rm(PModel_samp_fixed)
  PredN_samp_fixed_mt <- apply(NModel_samp_fixed , MARGIN = 2, function(Model_samp) {CPred_MM %*% Model_samp })
  rm(NModel_samp_fixed)
  if(verbose) {cat("Fixed effects prediction matrix size: \n PModel: ", format(object.size(PredP_samp_fixed_mt), units = "auto"), 
                   "\n NModel: ",format(object.size(PredN_samp_fixed_mt), units = "auto"), "\n\n")}
  
  # # Prediction by for both fixed and random effects
  # PredP_lk <- PredP_samp_fixed_mt + PModel_samp_rand_mt
  # PredN_lk <- PredN_samp_fixed_mt + NModel_samp_rand_mt
  gc()
  
  # Calculate the probability
  PredP <- ifelse(is.finite(exp(PredP_samp_fixed_mt + PModel_samp_rand_mt)/(1+exp(PredP_samp_fixed_mt + PModel_samp_rand_mt))),(exp(PredP_samp_fixed_mt + PModel_samp_rand_mt)/(1+exp(PredP_samp_fixed_mt + PModel_samp_rand_mt))), 1 )
  PredN <- ifelse(is.finite(exp(PredN_samp_fixed_mt + NModel_samp_rand_mt)/(1+exp(PredN_samp_fixed_mt + NModel_samp_rand_mt))),(exp(PredN_samp_fixed_mt + NModel_samp_rand_mt)/(1+exp(PredN_samp_fixed_mt + NModel_samp_rand_mt))), 1 )
  PredAll <- PredP * PredN
  rm(PredP_samp_fixed_mt , PModel_samp_rand_mt, PredN_samp_fixed_mt , NModel_samp_rand_mt)
  
  
  if(verbose) {cat("\n\n\nSpatialising predictions...\n")}
  Layer <- SpatUnits[[KMR]] %>% bind_cols(Woody_Clrtype = as.vector(R), Woody = as.vector(Response$N), SUID = Covs$SUID) %>% dplyr::filter(Woody > 0) %>% mutate(ActualProp =  ifelse(Woody > 0, Woody_Clrtype / Woody, NA)) # only for properties with woody cover
  PredictionsP <- rowMeans(PredP) %>% bind_cols(SUID = as.vector(Covs$SUID[which(Response$N>0)])) %>% rename("PredP" = "...1")
  PredictionsN <- rowMeans(PredN) %>% bind_cols(SUID = as.vector(Covs$SUID[which(Response$N>0)])) %>% rename("PredN" = "...1")
  PredictionsAll <- rowMeans(PredAll) %>% bind_cols(SUID = as.vector(Covs$SUID[which(Response$N>0)])) %>% rename("PredAll" = "...1")
  PredictionsCombined <- PredictionsP %>% left_join(PredictionsN, by = join_by(SUID == SUID)) %>% left_join(PredictionsAll, by = join_by(SUID == SUID))
  Layer <- Layer %>% left_join(PredictionsCombined, by = join_by(SUID == SUID)) %>% dplyr::select(-Shape_Length, -Shape_Area)
  if(verbose) {cat("Spatialised predictions: \n\n")
    print(head(Layer))}

  return(list(Layer = Layer, PModel = PModel, NModel = NModel))
}


# function to refit model from selected model ----
# To include control.compute = list(dic = TRUE, config = TRUE)

refit_model <- function(KMR = "CC", ClearType = 1, ModelDir = "output/models/"){

  if (ClearType == 1) {
    CT <- "Ag"
  } else if (ClearType == 2) {
    CT <- "In"
  } else if (ClearType == 3) {
    CT <- "Fo"
  }
  
  # read model selection results for both forward and backward selection
  SelModel_BC <- qread(paste0(ModelDir, "SelModel_", KMR, "_" , CT , "_BC.qs"))
  SelModel_FC <- qread(paste0(ModelDir, "SelModel_", KMR, "_" , CT , "_FC.qs"))
  
  # Find minimum DIC values recorded in the model selection steps (WHILE LOOP)
  MinDIC_BC <- SelModel_BC$DIC_ls[which.min(unlist(SelModel_BC$DIC_ls))]
  MinDIC_FC <- SelModel_FC$DIC_ls[which.min(unlist(SelModel_FC$DIC_ls))]
  MinDIC_BC_DIC <- unlist(MinDIC_BC)
  MinDIC_FC_DIC <- unlist(MinDIC_FC)
  
  # Select the best model (Forward or Backward) based on the minimum DIC values
  Best_Mod <- if_else(MinDIC_BC_DIC < MinDIC_FC_DIC, "BC", "FC")
  cat("Best Model for  ", KMR, ": ", Best_Mod, "\n")
  
  Best_Mod <- if(Best_Mod == "BC"){Best_Mod <- SelModel_BC} else {Best_Mod <- SelModel_FC}
  
  SpatUnits = Best_Mod$SpatUnits
  RespData = Best_Mod$RespData
  CovsCD = Best_Mod$CovsCD
  SA1sPoly = Best_Mod$SA1sPoly
  
  Model <- fit_model2(KMR, ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly, Verbose = FALSE)
  
  # return models
  output_FPath <- file.path(ModelDir, paste0("Model_", KMR, "_" , CT , ".qs"))
  qsave(Model, file = output_FPath)
  
  return(list(PModel = Model$PModel, NModel = Model$NModel, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly))
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
  SpatUnits <- SpatUnits[[KMR]] %>% list() %>% setNames(KMR)
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
  
  # Placeholder for model with error
  ERROR_ls <- list()
  
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
  
  else if (Direction == "Forward Complete" | Direction == "FC"){
    #'######################################
    ## Forward complete model selection ----
    #'######################################
    
    tic("Forward Complete Model Selection")
    # Get a list of Full explanatory variables and explanatory variables yet to be tested
    left_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
    Sel_explV_ls <- list()
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
    Current_DIC <- if(!is.null(ResultP$dic$dic) && !is.null(ResultN$dic$dic)){ResultP$dic$dic + ResultN$dic$dic} else {Inf}
    if(Current_DIC == Inf){ERROR_ls[length(ERROR_ls)+1] <- "H0"}
    Best_DIC <- Current_DIC
    DIC_ls[1] <- Current_DIC 
    dDIC_ls[1] <- Current_DIC ##**##
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- "H0" 
    
    toc(log = TRUE) # null model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", "1", "\n", sep = "")
    
    while(length(left_explV_ls) > 0){
      #### Step forward... ----
      tic(paste0("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      Best_explV <- left_explV_ls[1] 
      Best_DIC <- Inf 
      
      for(x in 1: length(left_explV_ls)){
        
        explV <- paste(c(Sel_explV_ls, left_explV_ls[x]), collapse=" + ")
        
        # fit clearing versus no clearing model 
        ForP_H1 <- as.formula(paste0(paste("P", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjP, scale.model = TRUE)"))
        tic(paste0("ExplV ", x, " of ", length(left_explV_ls)))
        ResultP <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        
        # fit proportion cleared|clearing model
        ForN_H1 <- as.formula(paste0(paste("Prop", explV, sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
        ResultN <- INLA_with_Retry(N_retry=N_retry, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
        toc(log = TRUE) # current model
        
        # Calculate DIC for the current model and compare it with the best model
        Current_DIC <- if(!is.null(ResultP$dic$dic) && !is.null(ResultN$dic$dic)){ResultP$dic$dic + ResultN$dic$dic} else {Inf}
        if(Current_DIC == Inf){ERROR_ls[length(ERROR_ls)+1] <- explV}
        cat("DIC = ", format(round(Current_DIC, 2)), "   ExplV: ", explV, "\n", sep = "")
        DIC_ls[length(DIC_ls)+1] <- Current_DIC 
        dDIC_ls[length(dDIC_ls)+1] <- Current_DIC + 2*length(unlist(strsplit(explV, "\\+"))) ##**##
        names(DIC_ls)[length(DIC_ls)] <- names(dDIC_ls)[length(dDIC_ls)] <- paste(explV) 
        
        # Update the best DIC if the DIC of the current model is better
        if(Current_DIC < Best_DIC){
          Best_DIC <- Current_DIC
          Best_explV <- left_explV_ls[x]
        } else if (Current_DIC == Inf && Best_DIC == Inf){
          ERROR_ls[length(ERROR_ls)+1] <- "H0"
        }
      }
      
      Sel_explV_ls <- c(Sel_explV_ls, Best_explV)
      left_explV_ls <- setdiff(Full_explV_ls, Sel_explV_ls)
      
      toc(log = TRUE) # Step
      cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "") 
      
      cat("Explanatory variable not tested yet: ", paste(left_explV_ls, collapse = ", "), "\n\n", sep = "")
    }
    
    Best_explvs <- names(DIC_ls)[which.min(DIC_ls)]
    if(Best_explvs != "H0"){
      Sel_explV_ls <- unlist(strsplit(Best_explvs, " \\+ "))
    }else{
      Sel_explV_ls <- 1
      warning("Null model selected!")
      }
    toc(log = TRUE) # Forward Complete Model Selection
  }
  
  else if (Direction == "Backward Complete" | Direction == "BC"){
    #'######################################
    ## Backward complete model selection ----
    #'######################################
    
    tic("Backward complete Model Selection")
    
    # Get a list of Full explanatory variables and start by selecting all
    Sel_explV_ls <- Full_explV_ls <- names(CP %>% dplyr::select(-SA1, -SUID, -SA1ID))
    # Placeholder for Selected explanatory variables and Best DIC
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
    Current_DIC <- if(!is.null(ResultP$dic$dic) && !is.null(ResultP$dic$dic)){ResultP$dic$dic + ResultN$dic$dic} else {Inf}
    if(Current_DIC == Inf){ERROR_ls[length(ERROR_ls)+1] <- explV}
    Best_DIC <- Current_DIC
    DIC_ls[1] <- Current_DIC
    dDIC_ls[1] <- Current_DIC + (2*length(Full_explV_ls)) ##**##
    names(DIC_ls)[1] <- names(dDIC_ls)[1] <- paste(unlist(Full_explV_ls), collapse=" + ")
    
    toc(log = TRUE) # Full Model
    cat("DIC=", format(round(Current_DIC, 2)), "   ExplV: ", paste(Full_explV_ls, collapse = ", "), "\n", sep = "")
    
    while(length(Sel_explV_ls) > 0 ){
      #### Step backward... ----
      tic(paste0("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls)))
      cat("Step ", length(Full_explV_ls)-length(Sel_explV_ls)+1 , " of ", length(Full_explV_ls), "\n", sep = "")
      Worst_explV <- Sel_explV_ls[1]
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
        Current_DIC <- if(!is.null(ResultP$dic$dic) && !is.null(ResultN$dic$dic)){ResultP$dic$dic + ResultN$dic$dic} else {Inf}
        if(Current_DIC == Inf){ERROR_ls[length(ERROR_ls)+1] <- explV}
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
      
      toc(log = TRUE) # Step
      cat("DIC=", format(round(Best_DIC, 2)), "   Selected explV: ", paste(Sel_explV_ls, collapse = " , "), "\n", sep = "")
      cat("Explanatory variable removed: ", paste(worst_explV_ls, collapse = ", "), "\n\n", sep = "")
    
    }
    
    Best_explvs <- names(DIC_ls)[which.min(DIC_ls)]
    
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
  ResultP <- INLA_with_Retry(N_retry=N_retry+2, Initial_Tlimit = Initial_Tlimit, ForP_H1, data = DataP, family = "binomial", Ntrials = Ntrials, control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  
  # fit proportion cleared|clearing model
  ForN_H1 <- as.formula(paste0(paste("Prop", paste(Sel_explV_ls, collapse = " + "), sep=" ~ "), " + f(SA1ID, model = 'bym', graph = AdjN, scale.model = TRUE)"))
  cat(deparse(ForN_H1), "\n")
  ResultN <- INLA_with_Retry(N_retry=N_retry+2, Initial_Tlimit = Initial_Tlimit, ForN_H1, data = DataN, family = "beta", control.inla = control.inla(control.vb = INLA::control.vb(enable = FALSE)), control.compute = list(dic = TRUE), control.predictor = list(compute = TRUE, link = 1), verbose = Verbose)
  toc(log = TRUE) # Total Time
  
  # filter CovCD to only include selected explanatory variables (For Predict function)
  CovsCD <- CovsCD[[KMR]] %>% dplyr::select(all_of(Sel_explV_ls)) %>%  list() %>% setNames(KMR)
  
  # return models
  if(!is.null(OutputDir)){
    Model <- list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly, DIC_ls = DIC_ls, ERROR_ls = ERROR_ls)
    output_FPath <- file.path(OutputDir, paste0("SelModel_", KMR, "_", if (ClearType == 1) {"Ag_"} else if (ClearType == 2) {"In_"} else if (ClearType == 3) {"Fo_"} else {"Error"}, Direction, ".qs"))
    qsave(Model, file = output_FPath)
  }
  return(list(PModel = ResultP, NModel = ResultN, KMR = KMR, ClearType = ClearType, SpatUnits = SpatUnits, RespData = RespData, CovsCD = CovsCD, SA1sPoly = SA1sPoly, DIC_ls = DIC_ls, ERROR_ls = ERROR_ls))
}


# Function to combine predictions across KMR and write output ----
Combine_Predictions <- function(ClearType = 1, Prediction_DIR = "output/predictions/", WRITE_SHP = TRUE, WRITE_DATA = TRUE){
  # ClearType: 1 = Agriculture, 2 = Infrastructure, 3 = Forestry
  # Prediction_DIR: Directory containing the prediction files, output from predict_model* function
  # WRITE_SHP: Write shapefile output (Default = TRUE: write output shapefile, FALSE: do not write output shapefile)
  # WRITE_DATA: Write .qs output (Default = TRUE: write output .qs, FALSE: do not write output .qs)

  if (ClearType == 1) {
    CT <- "Ag"
  } else if (ClearType == 2) {
    CT <- "In"
  } else if (ClearType == 3) {
    CT <- "Fo"
  }
  
  Pred_Flist <- list.files(Prediction_DIR, pattern = paste0("^Pred.*", CT ,"\\.qs$"), full.names = TRUE)
  
  Layer_list <- map(Pred_Flist, ~qread(.x)$Layer)
  
  Layer_comb <- do.call(rbind, Layer_list)
  names(Layer_comb) <- c("KMR" , "SA1" , "Wdy_Clr" , "Woody" , "SUID" , "ActlPrp" , "PredP" , "PredN" , "PredAll" , "geometry")
  st_geometry(Layer_comb) <- "geometry"
  
  if (WRITE_SHP) {
    SHP_Filename <- file.path(Prediction_DIR, paste0("Pred_", CT, ".shp"))
    st_write(Layer_comb, SHP_Filename, delete_layer = TRUE)
  }
  if (WRITE_DATA) {
    DATA_Filename <- file.path(Prediction_DIR, paste0("Pred_", CT, ".qs"))
    qsave(Layer_comb, DATA_Filename)
  }
  
  return(Layer_comb)
}

# Function to prepare Koala habitat data ----
Prep_Khab <- function(ZStats_Woody, KMR_name_df){
  # ZStats_Woody: output from section "Select all covariates for each clearing type" *(clearing type data filtering above)
  # KMR_List: A dataframe consist of KMR full name and KMR abbreviation
  # KMR <- names(ZStats_Woody)
  Khab_data <- imap(ZStats_Woody, function(ZStat_Woody, KMR_a) {ZStat_Woody %>% mutate(KMR_a = KMR_a, SUID = 1:n())}) %>% 
    do.call(rbind, .) %>% 
    left_join(KMR_name_df, by = join_by(KMR_a == KMR_a))
  return(Khab_data)
}

# Function to calculate the risk of Koala habitat loss ----
Get_Khab_loss_risk <- function(Pred_data, Khab_data){
  # Pred_data: .qs or .shp output from Combine_Predictions function
  # Khab_data: .qs output from Prep_Khab function. This should consist of KHAB, SUID and KMR
  
  Pred_Khab <- Pred_data %>% left_join(Khab_data, by = join_by(SUID == SUID, KMR == KMR)) %>% 
    mutate(Khab_P = sum.Khab/Woody, KhabRisk = PredAll*Khab_P)
  
  # Warning message if woody vegetation is greater than Koala habitat 
  if(nrow(Pred_Khab[Pred_Khab$Woody < Pred_Khab$Khab,])>0){message("Woody vegetation greater than Koala habitat!!! \n\n")}
  
  Pred_Khab <- Pred_Khab %>% 
    dplyr::select(Woody, Risk = PredAll, Khab_P, KhabRisk, KMR, SUID)
  return(Pred_Khab)
}

# Function to process model fitting and covariate extraction ----
Get_cov_coeff_long <- function(ClearType) {
  # ClearType
  if(ClearType == 1){
    CT <- "Ag"
  } else if(ClearType == 2){
    CT <- "In"
  } else if(ClearType == 3){
    CT <- "Fo"
  }
  
  # Read the required datasets based on the clear type - CT
  SUs <- qread(paste0("output/spatial_units/sus_", CT, ".qs"))
  ZStats_Woody <- qread(paste0("output/data/ZStats_Woody_", CT, ".qs"))
  ZStats_Covs <- qread(paste0("output/data/ZStats_Covs_", CT, ".qs"))
  SA1s <- qread("output/spatial_units/sa1s.qs")
  
  KMRs <- names(ZStats_Covs)
  kmr <- KMRs[1]
  
  # Fit the model using the specified CT
  Model <- fit_model2(KMR = kmr, ClearType = ClearType, SpatUnits = SUs, RespData = ZStats_Woody, CovsCD = ZStats_Covs, 
                      SA1sPoly = SA1s, Explanatory = "All", Verbose = FALSE, N_retry = 3, Initial_Tlimit = 1000, OutputDir = NULL)
  
  # Get summary of the model's fixed effects
  Cov_ls <- summary(Model$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>% 
    dplyr::select(Covariate) %>% arrange(Covariate)
  
  Cov_ls_long <- data.frame()
  
  # Loop through each KMR and gather covariate and coefficient data
  for (kmr in KMRs) {
    MODEL <- qread(paste0("output/models/Model_", kmr, "_", CT, ".qs"))
    
    Cov <- summary(MODEL$PModel)$fixed %>% as.data.frame() %>% rownames_to_column("Covariate") %>%
      dplyr::select(Covariate) %>% unlist()
    
    Cof_PModel <- summary(MODEL$PModel)$fixed[,1]
    Cof_PModel <- if_else(Cof_PModel > 0, 1, -1)
    
    Cof_NModel <- summary(MODEL$NModel)$fixed[,1]
    Cof_NModel <- if_else(Cof_NModel > 0, 1, -1)
    
    Cov_cof <- as.data.frame(cbind(Cov, Cof_PModel, Cof_NModel))
    rownames(Cov_cof) <- NULL
    
    Cov_ls_model <- left_join(Cov_ls, Cov_cof, by = join_by("Covariate" == "Cov")) %>% mutate(kmr = kmr)
    Cov_ls_long <- rbind(Cov_ls_long, Cov_ls_model)
  }
  
  Cof_PModel_CT <- paste0("Cof_PModel_", CT)
  Cof_NModel_CT <- paste0("Cof_NModel_", CT)
  
  # Change all NA to 0
  Cov_ls_long[[Cof_PModel_CT]] <- ifelse(is.na(Cov_ls_long$Cof_PModel), 0, Cov_ls_long$Cof_PModel)
  Cov_ls_long[[Cof_NModel_CT]] <- ifelse(is.na(Cov_ls_long$Cof_NModel), 0, Cov_ls_long$Cof_NModel)
  
  return(Cov_ls_long)
}

# Function to plot Single deforestation risk / Koala habitat loss ----
PLOTMAP_risk_NSW <- function(DATA, FILL, LEGEND_Title = "Deforestation\nrisk", ClearType = 1, FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300){
  # DATA: Spatial polygon data with deforestation risk or Koala habitat loss risk for each polygon
  # FILL: Variable name to be plotted
  # LEGEND_Title: Text for lebeling the 
  # ClearType: 1 = Ag, 2 = In, 3 = Fo
  # FilenamePath_PNG: Directory to save the plot as a PNG file (optional) {set to NULL to not save}
  
  CT <- if(ClearType == 1) {"Ag"} else if(ClearType == 2) {"In"} else if(ClearType == 3) {"Fo"} else {"Error"}
  
  # Plot the deforestation risk or Koala habitat loss risk
  Plot <- ggplot()+
    
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    
    geom_sf(data = DATA, aes(fill = {{FILL}}), color = NA)+
    geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Blues 3" ,rev = TRUE), name = LEGEND_Title)+
    
    # start a new scale
    new_scale_colour() +
    
    geom_sf(data = NSW_urb_sel_pt, colour = "red3", size = 1)+
    geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16), 
                    fontface = "bold", nudge_y = -5, size = 3,
                    color = "black",     # text color
                    bg.color = "grey90", # shadow color
                    bg.r = 0.05)+          # shadow radius
    
    ggspatial::annotation_scale(location = "br", pad_y = unit(1, "cm"))+
    ggspatial::annotation_north_arrow(location = "br", which_north = "true", pad_y = unit(2, "cm"))+
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = c(0.9, 0.3))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = TRUE)
  
  if(!is.null(FilenamePath_PNG)){
    if(tools::file_ext(FilenamePath_PNG) != "png"){warning("Filename extension should be '.png'")}
    ggsave(FilenamePath_PNG, Plot, width = PNG_width, height = PNG_height, dpi = PNG_dpi)
  }
  return(Plot)
}

# Function to plot Maps with Insets deforestation risk / Koala habitat loss ----
PLOTMAP_risk_with_Insets <- function(DATA, FILL, LEGEND_Title = "Koala habitat\nloss risk", ClearType = 1, 
                                     Inset_BL, Inset_dim = 100000,
                                     URB_PT_Main =NULL, URB_PT_SUB1, URB_PT_SUB2,  URB_PT_SUB3, 
                                     FilenamePath_PNG = NULL, PNG_width = 11, PNG_height = 11, PNG_dpi = 300){
  # DATA: Spatial polygon data with deforestation risk or Koala habitat loss risk for each polygon
  # ClearType: 1 = Ag, 2 = In, 3 = Fo
  # Inset_BL: Dataframe with x,y of the Bottom left corner of the inset map
  # Inset_dim: Dimension of the inset map
  # URB_PT_*: list of urban points label for main map, and sub maps (get exact name from NSW_urb_pt$UCL_NAME16)
  # Rename_URB_PT_*: list of names to replace existing name of the urban points label for main map, and sub maps
  # FilenamePath_PNG: Directory to save the plot as a PNG file (optional) {set to NULL to not save}
  
  CT <- if(ClearType == 1) {"Ag"} else if(ClearType == 2) {"In"} else if(ClearType == 3) {"Fo"} else {"Error"}
  
  Map_main <- ggplot() +
    
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    
    geom_sf(data = DATA, aes(fill = {{FILL}}), color = NA)+
    geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Reds 3" ,rev = TRUE), name = LEGEND_Title)+
    
    # start a new scale
    new_scale_colour() +
    
    # selected urban points
    geom_sf(data = NSW_urb_sel_pt, colour = "black", size = 1)+
    geom_text_repel(data = NSW_urb_sel_pt, aes(x = x, y = y , label = UCL_NAME16),
                    nudge_y = -5, size = 4, color = "black", bg.color = "grey90", bg.r = 0.05)+
    
    
    geom_rect(aes(xmin = Inset_BL[1,1], xmax = Inset_BL[1,1]+100000, ymin = Inset_BL[1,2], ymax = Inset_BL[1,2]+100000), fill = NA, colour = "black", linewidth = 0.6) +
    geom_rect(aes(xmin = Inset_BL[2,1], xmax = Inset_BL[2,1]+100000, ymin = Inset_BL[2,2], ymax = Inset_BL[2,2]+100000) , fill = NA, colour = "black", linewidth = 0.6) + 
    geom_rect(aes(xmin = Inset_BL[3,1], xmax = Inset_BL[3,1]+100000, ymin = Inset_BL[3,2], ymax = Inset_BL[3,2] + 100000), fill = NA, colour = "black", linewidth = 0.6) +
    
    geom_text(data = data.frame(x = Inset_BL[1,1] + 100000/2, y = Inset_BL[1,2] + 100000/2), aes(x = x, y = y, label = "A"), size = 4, fontface = "bold", color = "black", bg.color = "grey90", bg.r = 0.05)+
    geom_text(data = data.frame(x = Inset_BL[2,1] + 100000/2, y = Inset_BL[2,2] + 100000/2), aes(x = x, y = y, label = "B"), size = 4, fontface = "bold", color = "black", bg.color = "grey90", bg.r = 0.05)+
    geom_text(data = data.frame(x = Inset_BL[3,1] + 100000/2, y = Inset_BL[3,2] + 100000/2), aes(x = x, y = y, label = "C"), size = 4, fontface = "bold", color = "black", bg.color = "grey90", bg.r = 0.05)+
    
    ggspatial::annotation_scale(location = "bl", pad_y = unit(0.5, "cm"))+
    ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(1, "cm"))+
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = c(0.95, 0.1))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = st_bbox(KMR_shp)[c(1,3)], ylim = st_bbox(KMR_shp)[c(2,4)], expand = FALSE)
  
  Map_sub1 <- ggplot() +
    
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    
    geom_sf(data = DATA, aes(fill = {{FILL}}), color = NA)+
    # geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Reds 3" ,rev = TRUE), name = expression("Koala habitat\nloss risk"))+
    
    # start a new scale
    new_scale_colour() +
    
    geom_sf(data = (NSW_urb_pt %>% filter(UCL_NAME16 %in% URB_PT_SUB1)), colour = "black", size = 2)+
    geom_text_repel(data = (NSW_urb_pt %>% filter(UCL_NAME16 %in% URB_PT_SUB1)), aes(x = x, y = y , label = UCL_NAME16),
                    nudge_y = -5, size = 4, color = "black", bg.color = "grey90", bg.r = 0.05)+
    
    labs(tag = "A") +
    theme(plot.tag.position = "topleft", plot.tag.location = "panel") +
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = c(Inset_BL[1,1],  Inset_BL[1,1]+100000 ), ylim =  c( Inset_BL[1,2], Inset_BL[1,2]+100000), expand = FALSE)
  
  Map_sub2 <- ggplot() +
    
    geom_sf(data = STE, fill = "grey80", color = "white", lwd = 0.2)+
    
    geom_sf(data = DATA, aes(fill = {{FILL}}), color = NA)+
    # geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Reds 3" ,rev = TRUE), name = expression("Koala habitat\nloss risk"))+
    
    # start a new scale
    new_scale_colour() +
    
    geom_sf(data = (NSW_urb_pt %>% filter(UCL_NAME16  %in% URB_PT_SUB2)), colour = "black", size = 2)+
    geom_text_repel(data = (NSW_urb_pt %>% filter(UCL_NAME16 %in% URB_PT_SUB2)), aes(x = x, y = y , label = UCL_NAME16),
                    nudge_y = -5, size = 4, color = "black", bg.color = "grey90", bg.r = 0.05)+
    
    labs(tag = "B") +
    theme(plot.tag.position = "topleft", plot.tag.location = "panel") +
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = c(Inset_BL[2,1],  Inset_BL[2,1]+100000 ), ylim =  c( Inset_BL[2,2], Inset_BL[2,2]+100000), expand = FALSE)
  
  Map_sub3 <- ggplot() +
    
    geom_sf(data = STE, fill = "grey80", color = "black", lwd = 0.2)+
    
    geom_sf(data = DATA, aes(fill = {{FILL}}), color = NA)+
    # geom_sf(data = KMR_shp, fill = NA, color = "grey10", lwd = 0.2)+
    scale_fill_gradientn(colours = hcl.colors(8, palette = "Reds 3" ,rev = TRUE), name = expression("Koala habitat\nloss risk"))+
    
    # start a new scale
    new_scale_colour() +
    
    geom_sf(data = (NSW_urb_pt %>% filter(UCL_NAME16  %in% URB_PT_SUB3)), colour = "black", size = 2)+
    geom_text_repel(data = (NSW_urb_pt %>% filter(UCL_NAME16 %in% URB_PT_SUB3)), aes(x = x, y = y , label = UCL_NAME16),
                    nudge_y = -5, size = 4, color = "black", bg.color = "grey90", bg.r = 0.05)+
    
    labs(tag = "C") +
    theme(plot.tag.position = "topleft", plot.tag.location = "panel") +
    
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(), axis.line.x = element_blank())+
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank())+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_sf(xlim = c(Inset_BL[3,1],  Inset_BL[3,1]+100000 ), ylim =  c( Inset_BL[3,2], Inset_BL[3,2]+100000), expand = FALSE)

  Map_IN <- Map_sub1 + Map_sub2 + Map_sub3
  FINAL_map <- (Map_main / Map_IN) + plot_layout(widths = unit(c(19, 19), c("cm", "cm")), heights = unit(c(15.3, 6), c("cm", "cm")))
  
  ggsave(FilenamePath_PNG, FINAL_map, width = 21, height = 23, units = "cm", dpi = 300)

  return(FINAL_map)
}