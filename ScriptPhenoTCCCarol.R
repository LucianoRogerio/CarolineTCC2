###
## Separar os arquivos de dados em uma pasta chamada data e criar uma pasta output no seu diretório de trabalho
###

suppressWarnings(suppressMessages(library(tidyverse)))
library(reshape2); library(here)
DadosPhen <- read.table(here::here("data", "DadosBAG2022.csv"),
                        header = T, sep = ",", na.strings = "")

colnames(DadosPhen)[c(10:46, 52:61)] <- c("Stand", "PRC", "PRNC", "PTR", "PPA", "DMCsg",
                                          "DRY", "IC", "AP", "NMa", "PA", "Vigor45D",
                                          "Vigor12M", "RF", "PRz", "NRz", "PPD", "Pod",
                                          "PodRCas", "PodRPol", "FrogSkin", "CorPol",
                                          "Friab", "MoscBc", "MoscBr", "MurPA", "MosCm",
                                          "MosNv", "AcSev", "MFT", "MPs", "MBs", "QFs",
                                          "Anth", "BacSev", "FerSev", "SupASev",
                                          paste0("HCN-",1:9), "TCC")

DadosPhen <- DadosPhen %>% mutate(Cook = ifelse(Cook15 == 1,
                                                yes = "15",
                                                no = ifelse(Cook20 == 1,
                                                            yes = "20",
                                                            no = ifelse(Cook25,
                                                                        yes = "25",
                                                                        no = ifelse(Cook30 == 1,
                                                                                    yes = "30",
                                                                                    no = ifelse(Cook40 == 1,
                                                                                                yes = "40",
                                                                                                no = ifelse(Cook40 == 0,
                                                                                                            yes = ">40",
                                                                                                            no = NA))))))) %>% 
  mutate(Cook = factor(Cook, levels = c("15", "20", "25", "30", "40", ">40"), ordered = T)) %>%
  dplyr::select(-Cook15, -Cook20, -Cook25, -Cook30, -Cook40)

colnames(DadosPhen)[c(1:41,51:52)]

DadosParQuant <- DadosPhen[,c(1:30,33:45,56)]
DadosParQual <- DadosPhen[,c(1:9,31:32,57)]

DadosParQuantFin <- DadosParQuant %>% reshape2::melt(data = ., id.vars = c(1:10),
                                                     variable.name = "Trait", value.name = "y") %>%
  filter(!is.na(y)) %>%
  dplyr::mutate(Ano = Ano,
                Campo = Campo,
                Local = Local,
                trial = match(paste(Ano, Campo, Local, sep = "."),
                              unique(paste(Ano, Campo, Local, sep = "."))),
                clone = Genotipos.BGM,
                y = as.numeric(y), .keep = "unused")


DadosParQualFin <- DadosParQual %>% group_by(Genotipos.BGM) %>%
  summarise(CorPol = table(CorPol) %>%
              .[order(., decreasing = T)] %>% .[1] %>% names %>% as.integer(),
            Friab = table(Friab) %>%
              .[order(., decreasing = T)] %>% .[1] %>% names %>% as.character(),
            Cook = table(Cook) %>%
              .[order(., decreasing = T)] %>% .[1] %>% names %>% as.character())

saveRDS(DadosParQualFin, file = here::here("output", "DadosQualCaroline.rds"))

DadosHCN <- DadosPhen[,c(1:9, 47:55)]

DadosHCNFin <- DadosHCN %>% reshape2::melt(data = ., id.vars = c(1:9),
                                           variable.name = "HCN",  value.name = "y") %>%
  filter(!is.na(y)) %>%
  mutate(Ano = Ano,
         Campo = Campo,
         Local = Local,
         trial = match(paste(Ano, Campo, Local, sep = "."),
                       unique(paste(Ano, Campo, Local, sep = "."))),
         clone = Genotipos.BGM, .keep = "unused")

DataCar <- tibble(Traits = unique(DadosParQuantFin$Trait),
                  PhenData = NA)

for(i in DataCar$Traits){
  DataCar$PhenData[DataCar$Traits == i] <- list(DadosParQuantFin %>% filter(Trait == i))
}

DataCar <- DataCar %>% rbind(tibble(Traits = "HCN",
                                    PhenData = list(DadosHCNFin)))

saveRDS(object = DataCar, file = here::here("data", "DadosPhenCar.rds"))


suppressWarnings(suppressMessages(library(tidyverse)))
library(MuMIn)
library(reshape2); library(here)

suppressMessages(source(here::here("code", "MixedModelsFunctions.R")))

DadosCar <- readRDS(file = here::here("data", "DadosPhenCar.rds"))


require(furrr); plan(multisession, workers = 3)
DataCar <- DataCar %>% dplyr::mutate(TrialSel = future_map2(Traits, PhenData, function(Traits, PhenData,...){
  Trials <- unique(PhenData$trial)
  PhenData
  results <- tibble()
  for(i in Trials) {
    try(MixedModels <- analyzeTrial.lme4(PhenData %>% filter(trial %in% i)))
    try(result <- tibble(Trial = i,
                         Trait = Traits,
                         NClones = nrow(unique(PhenData %>%
                                                 filter(trial %in% i) %>% 
                                                 dplyr::select(clone))),
                         VarG = as.data.frame(VarCorr(MixedModels))[,c("grp","vcov")] %>% .[1,2],
                         VarE = as.data.frame(VarCorr(MixedModels))[,c("grp","vcov")] %>% .[2,2],
                         H2 = VarG/(VarG + VarE),
                         Real = suppressWarnings(MuMIn::r.squaredGLMM(MixedModels)[2])))
    try(results <- rbind(results, result))
  }
  return(results = results)
}))

Results <- DataCar %>% dplyr::select(TrialSel) %>% unnest_longer(TrialSel) %>% .[[1]]


library(reactable)
TrialsList <- unique(DadosParQuantFin[,c("Ano", "Campo", "Local", "trial")])

Results2 <- Results %>% right_join(TrialsList, by = c("Trial" = "trial")) %>%
  dplyr::select(Trial, Ano, Campo, Local, everything()) %>% mutate(Selecionado = ifelse(Real > 0.25 & H2 > 0.15, "Sim", "Nao"))

Results2 %>% reactable(groupBy = c("Trait"), columns = list(
  VarG = colDef(format = colFormat(digits = 2, locales = "en-US")),
  VarE = colDef(format = colFormat(digits = 2, locales = "en-US")),
  H2 = colDef(format = colFormat(digits = 3, locales = "en-US")),
  Real = colDef(format = colFormat(digits = 3, locales = "en-US"))))



DataSelPar <- DadosParQuantFin %>% mutate(Trait.Trial = paste(Trait, trial, sep = ".")) %>%
  .[.$Trait.Trial %in% (Results2 %>% mutate(Trait.Trial = paste(Trait, Trial, sep = ".")) %>%
                          filter(Selecionado == "Sim") %>% .$Trait.Trial),]


DataCar <- tibble(Traits = unique(DataSelPar$Trait),
                  PhenData = NA)

for(i in DataCar$Traits){
  DataCar$PhenData[DataCar$Traits == i] <- list(DataSelPar %>% filter(Trait == i))
}





DataSelHCN <- DadosHCNFin %>% mutate(Trait.Trial = paste("HCN", trial, sep = ".")) %>% 
  .[.$Trait.Trial %in% (Results2 %>% mutate(Trait.Trial = paste(Trait, Trial, sep = "."))
                        %>% filter(Selecionado == "Sim", NClones >= 20) %>% .$Trait.Trial),]

DataCar <- DataCar %>% rbind(tibble(Traits = "HCN",
                                    PhenData = list(DataSelHCN)))

saveRDS(object = DataCar, file = here::here("data", "DadosPhenSelCar.rds"))




library(here)
library(furrr)
library(tidyverse)
source(here::here("code", "MixedModelsFunctions.R"))
library(elliptic)

DataCar <- readRDS(here::here("data", "DadosPhenSelCar.rds"))

plan(multisession(workers = 3))

PhenData <- DataCar$PhenData[[27]]

DataCar <- DataCar %>% mutate(Blups = future_map2(Traits, PhenData, function(Traits, PhenData,...){
  RhpcBLASctl::blas_set_num_threads(3)
  print(paste("Trait", Traits, sep = " "))
  PhenData <- PhenData %>%
    mutate(trial = as.character(trial),
           rep = as.character(Bloco),
           repTrial = as.factor(paste(trial, Bloco, sep = ":")),
           LocYear = as.factor(paste(Local, Ano, sep = ":")))
  
  MM <- analyzeTrial.lme4Conj(PhenData)
  blups <- ranef(MM)$clone + fixef(MM)[[1]]
  Blups <- tibble(id = rownames(blups),
                  blups = blups$`(Intercept)`)
  colnames(Blups)[2] <- Traits
  file <- here::here("output", "MixedModels",
                     paste("Blups_", Traits, ".rds", sep = ""))
  saveRDS(object = Blups, file = file)
  return(Blups)
}))


BlupsTraits <- readRDS(here::here("output", "MixedModels", "Blups_HCN.rds")) %>% 
  rename(HCN = `35`)
traits2 <- DataCar$Traits %>% as.character() %>% .[. != "HCN"]

for(i in traits2){
  filename <- paste("Blups_", i, ".rds", sep = "")
  BlupsTraits <- BlupsTraits %>%
    full_join(., readRDS(here::here("output", "MixedModels", filename)), by = "id")
  colnames(BlupsTraits)[ncol(BlupsTraits)] <- i
}

saveRDS(object = BlupsTraits,
        file = here::here("output", "BlupsFenCar.rds"))



library(tidyverse); library(data.table)
Blups <- readRDS(here::here("output", "BlupsFenCar.rds"))

DataPhen <- readRDS(file = here::here("output", "DadosQualCaroline.rds"))


DataTCC <- Blups %>% full_join(DataPhen, by = c("id" = "Genotipos.BGM"))
write.table(DataTCC, file = here::here("output", "BlupsModaTCC.csv"),
            quote = F, row.names = F, sep = ",")

# Conferindo o número de dados perdidos por característica

TraitsSel <- names(colSums(is.na(DataTCC)))[colSums(is.na(DataTCC))/nrow(DataTCC) < 0.80]




TC1Traits <- c("id", "IC", "PTR", "NRz", "DRY", "DMCsg", "PPA", "PA", "PRz")

TC1Data <- DataTCC[, TC1Traits] %>% .[!(rowSums(is.na(.)) >= floor(ncol(.) - 1)/2),]

## Temática 2 - Resistência a Doenças Foliares
TC2Traits <- c("id", "Anth", "MPs", "Vigor45D", "FrogSkin",
               "MosNv", "RF", "AcSev")
TC2Data <- DataTCC[, TC2Traits] %>% .[!(rowSums(is.na(.)) >= floor(ncol(.) - 1)/2),]

## Temática 3 - Qualidade das raízes
# Verificar conjunto de dados de deteorização Fisiológica com Eder e predições de Carotenóides e DMC com Massaine
# Fazer uma coleção para Coloração Branca e uma para coloração amarela
TC3Traits <- c("id", "HCN", "DMCsg", "TCC", "CorPol")
TC3Data <- DataTCC[, TC3Traits] %>% .[!(rowSums(is.na(.)) >= floor(ncol(.) - 1)/2),]



library(reactable)

NaData1 <- tibble(NuclearCol = "RootYield",
                  Trait = colnames(TC1Data)[-1],
                  MissingData = colSums(is.na(TC1Data[,-1]))/nrow(TC1Data[,-1]) %>% as.vector)
NaData2 <- tibble(NuclearCol = "Disease",
                  Trait = colnames(TC2Data)[-1],
                  MissingData = colSums(is.na(TC2Data[,-1]))/nrow(TC2Data[,-1]) %>% as.vector)
NaData3 <- tibble(NuclearCol = "QualityRoot",
                  Trait = colnames(TC3Data)[-1],
                  MissingData = colSums(is.na(TC3Data[,-1]))/nrow(TC3Data[,-1]) %>% as.vector)

reactable::reactable(rbind(NaData1, NaData2, NaData3), columns = list(
  MissingData = colDef(format = colFormat(digits = 2, locales = "en-US", percent = T))))


InpDataFunc <- function(data){
  traits <- colnames(data)[-1]
  data <- as.data.frame(data)
  
  for(i in traits){
    if(class(data[,i]) == "character"){
      # Imputando a moda para dados qualitativos
      data[,i][is.na(data[,i])] <- table(data[,i]) %>%
        .[order(., decreasing = T)] %>% .[1] %>% names %>% as.character
    } else {
      # Imputando a media para dados quantitativos
      data[,i][is.na(data[,i])] <- mean(data[,i], na.rm = TRUE)
    }
  }
  return(tibble(data))
}

## TC1Data

TC1DataImp <- InpDataFunc(TC1Data)

saveRDS(TC1DataImp, file = here::here("output", "TCC1Data.rds"))

## TC2Data

TC2DataImp <- InpDataFunc(TC2Data)

saveRDS(TC2DataImp, file = here::here("output", "TCC2Data.rds"))

## TC3Data - Due to the difficult measure of the TCC, we applied a NIRS prediction
## for TCC to increase the number of clones evaluated for TCC, however the number
## of clones predicted by NIRS is also small compared to the BAG size.
## As TCC is close related to pulp color we will use this trait as a priori to
## impute the remaing clones.

TCCImputMeans <- TC3Data %>% group_by(CorPol) %>%
  summarise(TCCMeans = mean(TCC, na.rm = T)) %>% filter(!is.na(CorPol))
TCCImputMeans[4, 2] <- TCCImputMeans[c(3,5),]$TCCMeans %>% mean

TC3DataImp <- TC3Data %>%
  mutate(TCC = ifelse(test = !is.na(TCC),
                      yes = TCC,
                      no = ifelse(CorPol == 1,
                                  yes = TCCImputMeans[1,]$TCCMeans,
                                  no = ifelse(CorPol == 2,
                                              yes = TCCImputMeans[2,]$TCCMeans,
                                              no = ifelse(CorPol == 3,
                                                          yes = TCCImputMeans[3,]$TCCMeans,
                                                          no = ifelse(CorPol == 4,
                                                                      yes = TCCImputMeans[4,]$TCCMeans,
                                                                      no = TCCImputMeans[5,]$TCCMeans))))),
         CorPol = as.character(CorPol))

TC3DataImp <- InpDataFunc(TC3DataImp)

saveRDS(TC3DataImp, file = here::here("output", "TCC3Data.rds"))

TCCDataImp <- TC1DataImp %>% full_join(TC2DataImp, by = "id") %>% full_join(TC3DataImp, by = "id")
write.table(TCCDataImp, file = here::here("output", "BlupsModaTCCImp.csv"), 
            quote = F, row.names = F, sep = ",")
