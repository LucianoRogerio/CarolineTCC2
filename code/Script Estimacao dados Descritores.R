library(tidyverse)

#DadosPhen <- read.csv2("Dados fenotipicos BAG.CSV", header = T, na.strings = "-")

DadosPhen <- read.table("Teste.txt", header = T, sep = "\t", na.strings = "-")
traits <- colnames(DadosPhen)
traits[c(18, 12, 13, 20, 10, 14, 9, 8, 19, 11, 15, 25, 21)] <- c("PortePl","AltPl", "AltPlsF",
                                                           "Vigor", "PPA", "NMA",
                                                           "PTR", "RRC", "RRNC",
                                                           "IC", "DMCsg", "DMCoven", "RetF")
colnames(DadosPhen)[c(18, 12, 13, 20, 10, 14, 9, 8, 19, 11, 15, 25, 21)] <- c("PortePl","AltPl", "AltPlsF",
                                                             "Vigor", "PPA", "NMA",
                                                             "PTR", "RRC", "RRNC",
                                                             "IC", "DMCsg","DMCoven", "RetF")

TraitSel <- c("PortePl","AltPl", "AltPlsF",
              "Vigor", "PPA", "NMA",
              "PTR", "RRC", "RRNC",
              "IC", "DMCsg", "DMCoven", "RetF",
              "PLOT", "Ano", "Campo",
              "Local", "Genotipos.BGM.", "Genotipos",
              "Bloco")

DadosPhen2 <- DadosPhen[,colnames(DadosPhen)%in%TraitSel]

moda <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

DadosPhen2 %>% group_by(Genotipos) %>% summarize(Genotipos.BGM = unique(Genotipos.BGM.),
                                                 PortePlMean = mean(PortePl, na.rm = TRUE),
                                                 PortePlsd = sd(PortePl, na.rm = TRUE),
                                                 PortePlMod = moda(PortePl),
                                                 AltPlMean = mean(AltPl, na.rm = TRUE),
                                                 AltPlsd = sd(AltPl, na.rm = TRUE),
                                                 AltplsFMean = mean(AltPlsF, na.rm = TRUE),
                                                 AltplsFsd = sd(AltPlsF, na.rm = TRUE),
                                                 VigorMean = mean(Vigor, na.rm = TRUE),
                                                 Vigorsd = sd(Vigor, na.rm = TRUE),
                                                 VigorMod = moda(Vigor),
                                                 PPAMean = mean(PPA, na.rm = TRUE),
                                                 PPAsd = sd(PPA, na.rm = TRUE),
                                                 NMAMean = mean(NMA, na.rm = TRUE),
                                                 NMAsd = sd(NMA, na.rm = TRUE),
                                                 PTRMean = mean(PTR, na.rm = TRUE),
                                                 PTRsd = sd(PTR, na.rm = TRUE),
                                                 RRCMean = mean(RRC, na.rm = TRUE),
                                                 RRCsd = sd(RRC, na.rm = TRUE),
                                                 RRNCMean = mean(RRNC, na.rm = TRUE),
                                                 RRNCsd = sd(RRNC, na.rm = TRUE),
                                                 ICMean = mean(IC, na.rm = TRUE),
                                                 ICsd = sd(IC, na.rm = TRUE),
                                                 DMCsgMean = mean(DMCsg, na.rm = TRUE),
                                                 DMCsgsd = sd(DMCsg, na.rm = TRUE),
                                                 DMCovenMean = mean(DMCoven, na.rm = TRUE),
                                                 DMCovensd = sd(DMCoven, na.rm = TRUE),
                                                 RetFMean = mean(RetF, na.rm = TRUE),
                                                 RetFsd = sd(RetF, na.rm = TRUE))-> Res
head(Res)

write.csv(Res, "DadosProntos.csv", quote = FALSE, row.names = FALSE)

TraitHCNSel <- c("PLOT", "Ano", "Campo",
              "Local", "Genotipos.BGM.", "Genotipos",
              "Bloco")
TraitHCNSel <- c(TraitHCNSel, colnames(DadosPhen)[26:34])
TraitHCNSel[8:16] <- paste("HCN", 1:9, sep ="-")

colnames(DadosPhen)[26:34] <- paste("HCN", 1:9, sep ="-")
DadosHCN <- DadosPhen[,colnames(DadosPhen)%in%TraitHCNSel]


DadosHCN1 <- DadosHCN[,-c(9:16)]
DadosHCN2 <- DadosHCN[,-c(8, 10:16)]
DadosHCN3 <- DadosHCN[,-c(8:9, 11:16)]
DadosHCN4 <- DadosHCN[,-c(8:10, 12:16)]
DadosHCN5 <- DadosHCN[,-c(8:11, 13:16)]
DadosHCN6 <- DadosHCN[,-c(8:12, 14:16)]
DadosHCN7 <- DadosHCN[,-c(8:13, 15:16)]
DadosHCN8 <- DadosHCN[,-c(8:14, 16)]
DadosHCN9 <- DadosHCN[,-c(8:15)]

DadosHCN1$Rep <- rep(1, times = nrow(DadosHCN))
DadosHCN2$Rep <- rep(2, times = nrow(DadosHCN))
DadosHCN3$Rep <- rep(3, times = nrow(DadosHCN))
DadosHCN4$Rep <- rep(4, times = nrow(DadosHCN))
DadosHCN5$Rep <- rep(5, times = nrow(DadosHCN))
DadosHCN6$Rep <- rep(6, times = nrow(DadosHCN))
DadosHCN7$Rep <- rep(7, times = nrow(DadosHCN))
DadosHCN8$Rep <- rep(8, times = nrow(DadosHCN))
DadosHCN9$Rep <- rep(9, times = nrow(DadosHCN))

colnames(DadosHCN1)[8] <- "HCN"
colnames(DadosHCN2)[8] <- "HCN"
colnames(DadosHCN3)[8] <- "HCN"
colnames(DadosHCN4)[8] <- "HCN"
colnames(DadosHCN5)[8] <- "HCN"
colnames(DadosHCN6)[8] <- "HCN"
colnames(DadosHCN7)[8] <- "HCN"
colnames(DadosHCN8)[8] <- "HCN"
colnames(DadosHCN9)[8] <- "HCN"


DadosHCNCom <- rbind(DadosHCN1, DadosHCN2, DadosHCN3,
                     DadosHCN4, DadosHCN5, DadosHCN6,
                     DadosHCN7, DadosHCN8, DadosHCN9)

DadosHCNCom %>% group_by(Genotipos) %>% summarise(GenotiposBGM = unique(Genotipos.BGM.),
                                                  HCNMean = mean(HCN, na.rm = TRUE),
                                                  HCNsd = sd(HCN, na.rm = TRUE),
                                                  HCNmoda = moda(HCN)) -> Res2

write.csv(Res2, "DadosHCNPronts.csv", quote = FALSE, row.names = F)

Res %>% left_join(Res2, by = c("Genotipos" = "Genotipos")) %>% select(-GenotiposBGM)-> Res3

write.csv(Res3, "DadosProntos.CSV", quote = F, row.names = F)

#### Colecao Nuclear com dados fenotipicos
library(corehunter); library(StatMatch)

DataCar <- read.table(here::here("data", "DadosSelCaroline.CSV"), header = T, sep = ";", dec = ".",
                      na.strings = "NA", colClasses =  c("factor",
                                                        rep("character", times = 24),
                                                        rep("numeric", times = 27)))

DistCar <- gower.dist(data.x = DataCar[,-1])
row.names(DistCar) <- colnames(DistCar) <- DataCar$Acessos
DistCar2 <- distances(DistCar)

CCC <- sampleCore(DistCar2, objective(type = "EN", measure = "PD"), size = 338)

DataCarCCPheno1 <- DataCar[DataCar$Acessos%in%CCC$sel,]

write.table(DataCarCCPheno1, "DadosPheno1CCPheno1.CSV", quote = F, sep = ";", dec = ".")

colnames(DataCar)[1] <- "NAME"
DataCar$NAME <- as.character(DataCar$NAME)
DataCarP <- phenotypes(data = as.data.frame(DataCar[,-1]))
core <- sampleCore(DataCarP, objective(type = "EN", measure = "GD"), size = 338)

DataCarCCPheno2 <- DataCar[core$sel,]

write.table(DataCarCCPheno2, file = "DadosPheno2CCPheno1.CSV", quote = F, sep = ";", dec = ".")


library(phangorn)
DistCar <- as.dist(DistCar)
writeDist(x = DistCar, file = "MatrizGowerCaroline.phy", format = "phylip", upper = T, diag = T)

NSelDARwin <- read.csv(file = "ClonesNSelCCDARwin.csv", header = T)
DataCarCCPheno3 <- DataCar[-NSelDARwin$Clone,]

write.table(DataCarCCPheno3, file= "DadosPheno3CCPheno1.CSV", quote = F, sep = ";", dec = ".")


#### Boxplot de colecoes nucleares

Method1 <- read.table("DadosPheno1CCPheno1.CSV", header = T, sep = ";", dec = ".",
                     colClasses = c(rep("factor", times = 36),
                                    rep("numeric", times = 16)))
Method1$HCNMod <- as.factor(Method1$HCNMod)
Method1$DMCsg <- as.numeric(as.character(Method1$DMCsg))


Method2 <- read.table("DadosPheno2CCPheno1.CSV", header = T, sep = ";", dec = ".",
                      colClasses = c(rep("factor", times = 36),
                                     rep("numeric", times = 16)))
Method2$HCNMod <- as.factor(Method2$HCNMod)
Method2$DMCsg <- as.numeric(as.character(Method2$DMCsg))


Method3 <- read.table("DadosPheno3CCPheno1.CSV", header = T, sep = ";", dec = ".",
                      colClasses = c(rep("factor", times = 36),
                                     rep("numeric", times = 16)))
Method3$HCNMod <- as.factor(Method3$HCNMod)
Method3$DMCsg <- as.numeric(as.character(Method3$DMCsg))


AllBAG <- read.table("DadosSelCaroline.CSV", header = T, sep = ";", dec = ".",
                     colClasses = c(rep("factor", times = 36),
                                    rep("numeric", times = 16)))
AllBAG$HCNMod <- as.factor(AllBAG$HCNMod)
AllBAG$DMCsg <- as.numeric(as.character(AllBAG$DMCsg))


colnames(Method1) == colnames(Method2)

AllBAG$Data <- rep("BAG", times = nrow(AllBAG))
AllBAG$Method <- rep("BAG", times = nrow(AllBAG))

Method1$Data <- rep("Pheno", times = nrow(Method1))
Method1$Method <- rep("Gower + CH", times = nrow(Method1))

Method2$Data <- rep("Pheno", times = nrow(Method2))
Method2$Method <- rep("CoreHunter", times = nrow(Method2))

Method3$Data <- rep("Pheno", times = nrow(Method3))
Method3$Method <- rep("Gower + DARwin", times = nrow(Method3))

Alldataset <- rbind(AllBAG, Method1, Method2, Method3)

library(ggplot2); library(reshape2)
QualityTrait <- colnames(Alldataset)[2:36]


AlldataSetQuant <- Alldataset %>% select(-QualityTrait)
AlldataSetQuali <- Alldataset %>% select(Acessos, QualityTrait, Data, Method)

AlldataSetQuant2 <- melt(AlldataSetQuant, id.vars = c("Acessos", "Data", "Method"),
       variable.name = "Trait", value.name = "Value", )
head(AlldataSet2)

filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]

  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

AlldataSetQuant3 <- AlldataSetQuant2 %>% group_by(Trait, Method) %>%
  mutate(Value2 = filter_lims(Value))

AlldataSetQuant3 %>% ggplot(mapping = aes(y = Value2, color = Data, x = Method)) +
  geom_boxplot() + facet_wrap(~Trait, scale = "free")


write.table(as.character(AllBAG$Acessos), file = "ClonesFenotipadosCNCaroline.csv",
            sep = "\t", dec = ".", quote = F, row.names = F)

FenNames <- read.table(file = "ClonesFenotipadosCNCaroline.csv", header = T, sep = "\t")

GBSNames <- read.table("GBSNames.CSV", header = T, sep = ";")

GBSNamesSel <- GBSNames[GBSNames$CloneUpd%in%FenNames$Clones,]
rownames(GBSNamesSel) <- 1:nrow(GBSNamesSel)
GBSDataCar <- GBSdata[row.names(GBSdata)%in%GBSNamesSel$FullSampleName,]

GBSRow <- data.frame(Clones = row.names(GBSDataCar), rank = 1:nrow(GBSDataCar))
GBSNamesSel %>% left_join(GBSRow, by = c("FullSampleName" = "Clones")) %>%
  arrange(rank) -> GBSNamesSel2


row.names(GBSDataCar) <- GBSNamesSel2$CloneUpd
GBSDataCar2 <- GBSDataCar[!duplicated(row.names(GBSDataCar)),]

saveRDS(GBSDataCar2, file = "data/DadosGBSCaroline.RDS")

GBSDataCar <- readRDS("data/DadosGBSCaroline.RDS")

Freq <- colMeans(GBSDataCar, na.rm = T)

MAF <- Freq
for(i in 1:length(MAF)){
  if(MAF[i] >0.5){MAF[i] <- 1 - MAF[i]} else {MAF[i] <- MAF[i]}
}

GBSDataCar2 <- as.data.frame(GBSDataCar[,(MAF >= 0.05)])
GBSDataCar3 <- GBSDataCar2

GBSDataCar3[GBSDataCar3==0] <- "0/0"
GBSDataCar3[GBSDataCar3==1] <- "0/1"
GBSDataCar3[GBSDataCar3==2] <- "1/1"

GBSDataCar3[1:10,1:20]

M_Darwin_Car <- apply(GBSDataCar3, 1, function(x) do.call(cbind, strsplit(x, split = "/")))

dim(M_Darwin_Car)

rownames(M_Darwin_Car) <- rep(colnames(GBSDataCar3), each=2)

M_Darwin_Car[1:10,1:10]

write.table(t(M_Darwin_Car), file = "output/CarDarwinMatrix.txt", quote = F,
            sep = "\t")

rm(M_Darwin_Car); rm(GBSDataCar); rm(GBSDataCar2); rm(GBSDataCar3)
##### Core Hunter
library(corehunter); library(tidyverse); library(data.table)

CHCar <- genotypes(GBSDataCar2, format = "biparental")

## Shannon optimization - Modified Rogers Distance
core3 <- sampleCore(CHCar, objective(type = "SH", measure = "MR"), size = 338)
CloneSel <- core3$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "í", replacement = "i") %>%
  gsub(pattern = "CPAFRO.35", replacement = "CPAFRO-35=BRS-COLONIAL")

core3$sel[core3$sel%like%"Aip"]
DataCar$Acessos[DataCar$Acessos%like%"Gau"]

CloneSel[!CloneSel%in%DataCar$Acessos]


DataCarCCGeno4 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno4, file= "DadosGeno1CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Shannon optimization - Cavalli-Sforza and Edwards Distance
core4 <- sampleCore(CHCar, objective(type = "SH", measure = "CE"), size = 338)
CloneSel <- core4$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "í", replacement = "i") %>%
  gsub(pattern = "CPAFRO.35", replacement = "CPAFRO-35=BRS-COLONIAL")


CloneSel[!CloneSel%in%DataCar$Acessos]

DataCarCCGeno5 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno5, file= "DadosGeno2CCGeno1.CSV", quote = F, sep = ";", dec = ".")



## Expected Proportion of Heterozigous Loci optimization - Modified Rogers Distance
core5 <- sampleCore(CHCar, objective(type = "HE", measure = "MR"), size = 338)
CloneSel <- core5$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "í", replacement = "i") %>%
  gsub(pattern = "CPAFRO.35", replacement = "CPAFRO-35=BRS-COLONIAL")

CloneSel[!CloneSel%in%DataCar$Acessos]

DataCar$Acessos[DataCar$Acessos%like%"Vass"]

DataCarCCGeno6 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno6, file= "DadosGeno3CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Expected Proportion of Heterozigous Loci optimization - Cavalli-Sforza and Edwards Distance
core6 <- sampleCore(CHCar, objective(type = "HE", measure = "CE"), size = 338)
CloneSel <- core6$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "í", replacement = "i") %>%
  gsub(pattern = "CPAFRO.35", replacement = "CPAFRO-35=BRS-COLONIAL")

CloneSel[!CloneSel%in%DataCar$Acessos]

DataCarCCGeno7 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno7, file= "DadosGeno4CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Allele coverage optimization - Modified Rogers Distance
core7 <- sampleCore(CHCar, objective(type = "CV", measure = "MR"), size = 338)
CloneSel <- core7$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "CPAFRO-05", replacement = "CPAFRO-05=ENTALA-GATO") %>%
  gsub(pattern = "CPAFRO-09", replacement = "CPAFRO-09=PAO-DO-ACRE")

CloneSel[!CloneSel%in%DataCar$Acessos]

DataCarCCGeno8 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno8, file= "DadosGeno5CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Allele coverage optimization - Cavalli-Sforza and Edwards Distance
core8 <- sampleCore(CHCar, objective(type = "CV", measure = "CE"), size = 338)
CloneSel <- core8$sel %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "CPAFRO.35", replacement = "CPAFRO-35=BRS-COLONIAL") %>%
  gsub(pattern = "CPAFRO.", replacement = "CPAFRO-") %>%
  gsub(pattern = "CPAFRO-03", replacement = "CPAFRO-03=PIRARUCU") %>%
  gsub(pattern = "CPAFRO-09", replacement = "CPAFRO-09=P?O-DO-ACRE") %>%
  gsub(pattern = "CPAFRO-12", replacement = "CPAFRO-12=CIP?") %>%
  gsub(pattern = "CPAFRO-18", replacement = "CPAFRO-18=IM-946") %>%
  gsub(pattern = "CPAFRO-28", replacement = "CPAFRO-28=PARATI") %>%
  gsub(pattern = "CPAFRO-40", replacement = "CPAFRO-40=EAB-451")

CloneSel[!CloneSel%in%DataCar$Acessos]
DataCar$Acessos[DataCar$Acessos%like%"CPAFRO"]

DataCarCCGeno9 <- DataCar[DataCar$Acessos%in%CloneSel,]

write.table(DataCarCCGeno9, file= "DadosGeno6CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Allele Match Distance for MLST - DARwin
saveRDS(GBSDataCar2, file = "data/DadosGBSCaroline.RDS")

GBSDataCar <- readRDS("data/DadosGBSCaroline.RDS")

Freq <- colMeans(GBSDataCar, na.rm = T)

MAF <- Freq
for(i in 1:length(MAF)){
  if(MAF[i] >0.5){MAF[i] <- 1 - MAF[i]} else {MAF[i] <- MAF[i]}
}

GBSDataCar2 <- as.data.frame(GBSDataCar[,(MAF >= 0.05)])
GBSDataCar3 <- GBSDataCar2

GBSDataCar3[GBSDataCar3==0] <- "0/0"
GBSDataCar3[GBSDataCar3==1] <- "0/1"
GBSDataCar3[GBSDataCar3==2] <- "1/1"

GBSDataCar3[1:10,1:20]

M_Darwin_Car <- apply(GBSDataCar3, 1, function(x) do.call(cbind, strsplit(x, split = "/")))

dim(M_Darwin_Car)

rownames(M_Darwin_Car) <- rep(colnames(GBSDataCar3), each=2)

M_Darwin_Car[1:10,1:10]

write.table(t(M_Darwin_Car), file = "output/CarDarwinMatrix.txt", quote = F,
            sep = "\t")

rm(M_Darwin_Car); rm(GBSDataCar); rm(GBSDataCar2); rm(GBSDataCar3)

CCSelAM.DARwin <- read.table(here::here("output", "CarDarwinMatrix.don"),
                             sep = "\t", header = T,skip = 2,
                             nrows = 1363)

ClonesSel <- CCSelAM.DARwin %>% filter(CC338 == "Kept") %>%
  select(Nº)
ClonesSel <- ClonesSel$Nº %>%
  gsub(pattern = "BGM.", replacement = "BGM-") %>%
  gsub(pattern = "BRS.", replacement = "BRS-") %>%
  gsub(pattern = "Aipim.", replacement = "Aipim-") %>%
  gsub(pattern = "SE.", replacement = "SE-") %>%
  gsub(pattern = "X2011.24.", replacement = "2011-24-") %>%
  gsub(pattern = "X2011.34.", replacement = "2011-34-") %>%
  gsub(pattern = "BRS-Poti.Branca", replacement = "BRS-Poti-Branca") %>%
  gsub(pattern = "BGM-1277.1", replacement = "BGM-1277") %>%
  gsub(pattern = "BGM-1550.1", replacement = "BGM-1550") %>%
  gsub(pattern = "BGM-1616.1", replacement = "BGM-1616") %>%
  gsub(pattern = "BGM-1680.1", replacement = "BGM-1680") %>%
  gsub(pattern = "CPAFRO.03", replacement = "CPAFRO-03=PIRARUCU") %>%
  gsub(pattern = "CPAFRO.18", replacement = "CPAFRO-18=IM-946") %>%
  gsub(pattern = "CPAFRO.29", replacement = "CPAFRO-29") %>%
  gsub(pattern = "CPAFRO.33", replacement = "CPAFRO-33") %>%
  gsub(pattern = "CPAFRO.42", replacement = "CPAFRO-42")

ClonesSel[ClonesSel%in%DataCar$Acessos]
DataCarCCGeno10 <- DataCar[DataCar$Acessos%in%ClonesSel,]

write.table(DataCarCCGeno10, file= "DadosGeno7CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## Modified Rodgers distances for MLST - DARwin
library(poppr); library(adegenet); library(here); library(tidyverse)
GBSDataCar <- readRDS("data/DadosGBSCaroline.RDS")

Freq <- colMeans(GBSDataCar, na.rm = T)

MAF <- Freq
for(i in 1:length(MAF)){
  if(MAF[i] >0.5){MAF[i] <- 1 - MAF[i]} else {MAF[i] <- MAF[i]}
}

GBSDataCar2 <- as.data.frame(GBSDataCar[,(MAF >= 0.05)])
rm(GBSDataCar)
GenIdCar <- df2genind(GBSDataCar2, ploidy = 2, type = "PA", ncode = 1)


MRDist <- rogers.dist(GenIdCar)
MRDist <- MRDistance
writeDist(x = MRDist, file = here::here("output","MatrizMDCaroline.phy"),
          format = "phylip", upper = F, diag = T)

CCSelMR.DARwin <- read.table(here::here("output", "MatrizMDCaroline.DON"),
                             sep = "\t", header = T,skip = 2,
                             nrows = 1336)

ClonesSel <- CCSelMR.DARwin %>% filter(Tree.Sample_338 == "Kept") %>%
  select(Name)
ClonesSel <- ClonesSel$Name %>%
             gsub(pattern = "CPAFRO-03", replacement = "CPAFRO-03=PIRARUCU") %>%
             gsub(pattern = "CPAFRO-18", replacement = "CPAFRO-18=IM-946") %>%
             gsub(pattern = "CPAFRO-28", replacement = "CPAFRO-28=PARATI")
DataCarCCGeno11 <- DataCar[DataCar$Acessos%in%ClonesSel,]

write.table(DataCarCCGeno11, file= "DadosGeno8CCGeno1.CSV", quote = F, sep = ";", dec = ".")

## for MLST - DARwin

library(poppr); library(adegenet); library(here); library(phangorn)
library(data.table)
CEDist <- edwards.dist(GenIdCar)

writeDist(x = CEDist, file = here::here("output","MatrizCECaroline.phy"),
          format = "phylip", upper = F, diag = T)

CCSelCE.DARwin <- read.table(here::here("output","MatrizCECaroline.DON"),
                             sep = "\t", header = T, skip = 2,
                             nrows = 1336)
ClonesSel <- CCSelCE.DARwin %>% filter(Tree.Sample_338 == "Kept") %>%
  select(Name)
ClonesSel <- ClonesSel$Name %>%
  gsub(pattern = "CPAFRO-03", replacement = "CPAFRO-03=PIRARUCU") %>%
  gsub(pattern = "CPAFRO-05", replacement = "CPAFRO-05=ENTALA-GATO") %>%
  gsub(pattern = "CPAFRO-18", replacement = "CPAFRO-18=IM-946")

DataCarCCGeno12 <- DataCar[DataCar$Acessos%in%ClonesSel,]

write.table(DataCarCCGeno12, file= "DadosGeno9CCGeno1.CSV", quote = F, sep = ";", dec = ".")
