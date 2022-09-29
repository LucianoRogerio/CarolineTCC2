
library(tidyverse); library(reshape2)

DataQual  <- read.csv2("DadosCarolineQualitativos.csv", header = T, dec = ".")
DataQuan1 <- read.csv2("DadosCarolineQuantitativos1.csv", header = T, dec = ".")
DataQuan2 <- read.csv2("DadosCarolineQuantitativos2.csv", header = T, dec = ".")

source("Shannon-Weaver.R")


Caract <- c("CorFolhaApical" = 4, "Pubescencia" = 2, "FormaLobuloCentral" = 9,
            "CorPeciolo" = 6, "CorFolhaDesenv" = 5, "NLobulos" = 5,
            "PosicaoPeciolo" = 4, "CorNerv" = 4, "Sinuosidade" = 2,
            "CorCortexCaule" = 4, "CorExternaCaule" = 7, "ComprFilotaxia" = 3,
            "CorEpidCaule" = 4, "HabCresCaule" = 2, "CorRamosTerm" = 3,
            "TipoPlan" = 4, "PorteMod" = 5, "ComprEstip" = 2,
            "MargEstip" = 2, "HabRamif" = 4, "AnguloRamif" = 22,
            "NiveisRam" = 4, "ProemCicatrizesFolhas" = 2, "FormaRzs" = 4,
            "CorExternaRzs" = 4, "CorCortexRzs" = 4, "CorPolpaRzs" = 5,
            "TxtEpidermeRzs" = 2, "PresPedunculoRzs" = 3, "ConstrRzs" = 3,
            "DestqPelicRz" = 2, "DestqCortexRz" = 3, "PosicaoRzs" = 3,
            "Flowering" = 2, "VigorMod" = 3, "HCNMod" = 9)


DataQual %>% select(-Acessos) %>% apply(FUN = Shannon.Weaver.QL, nclass = Caract, 2) %>%
  diag() %>% as.data.frame -> Res
row.names(Res) <- DataQual %>% select(-Acessos) %>% colnames

write.csv2(Res, "IndicesShannonWeaverQualCaroline.csv", quote = F, row.names = T)

head(DataQuan1)

DataQuan11 <- DataQuan1 %>% dplyr::select(Acessos, ComprLobuloPl1, ComprLobuloPl2, ComprLobuloPl3,
                                   ComprLobuloPl4, ComprLobuloPl5)

DataQuan12 <- DataQuan1 %>% dplyr::select(Acessos, LargLobulo1, LargLobulo2, LargLobulo3,
                                   LargLobulo4, LargLobulo5)
DataQuan13 <- DataQuan1 %>% dplyr::select(Acessos, RelComprLarPl1, RelComprLarPl2, RelComprLarPl3,
                                   RelComprLarPl4, RelComprLarPl5)
DataQuan14 <- DataQuan1 %>% dplyr::select(Acessos, ComprPecioloPl1, ComprPecioloPl2, ComprPecioloPl3,
                                   ComprPecioloPl4, ComprPecioloPl5)
DataQuan15 <- DataQuan1 %>% dplyr::select(Acessos, EspessuraEntreCascaPl1, EspessuraEntreCascaPl2,
                                   EspessuraEntreCascaPl3)
DataQuan16 <- DataQuan1 %>% dplyr::select(Acessos, ComprMedRzsPl1, ComprMedRzsPl2, ComprMedRzsPl3)
DataQuan17 <- DataQuan1 %>% dplyr::select(Acessos, DiamMedRzsPl1, DiamMedRzsPl2, DiamMedRzsPl3)
DataQuan18 <- DataQuan1 %>% dplyr::select(Acessos, NRPl1, NRPl2, NRPl3)
DataQuan19 <- DataQuan1 %>% dplyr::select(Acessos, FRYPl1, FRYPl2, FRYPl3)

colnames(DataQuan11)[2:6] <- 1:5
colnames(DataQuan12)[2:6] <- 1:5
colnames(DataQuan13)[2:6] <- 1:5
colnames(DataQuan14)[2:6] <- 1:5
colnames(DataQuan15)[2:4] <- 1:3
colnames(DataQuan16)[2:4] <- 1:3
colnames(DataQuan17)[2:4] <- 1:3
colnames(DataQuan18)[2:4] <- 1:3
colnames(DataQuan19)[2:4] <- 1:3


DataQuan11F <- DataQuan11 %>% melt(id.vars = "Acessos", value.name = "ComprLobulo", variable.name = "Plant")
DataQuan12F <- DataQuan12 %>% melt(id.vars = "Acessos", value.name = "LargLobulo", variable.name = "Plant")
DataQuan13F <- DataQuan13 %>% melt(id.vars = "Acessos", value.name = "RelComprLar", variable.name = "Plant")
DataQuan14F <- DataQuan14 %>% melt(id.vars = "Acessos", value.name = "ComprPeciolo", variable.name = "Plant")
DataQuan15F <- DataQuan15 %>% melt(id.vars = "Acessos", value.name = "EspessuraEntreCasca", variable.name = "Plant")
DataQuan16F <- DataQuan16 %>% melt(id.vars = "Acessos", value.name = "ComprMedRzs", variable.name = "Plant")
DataQuan17F <- DataQuan17 %>% melt(id.vars = "Acessos", value.name = "DiamMedRzs", variable.name = "Plant")
DataQuan18F <- DataQuan18 %>% melt(id.vars = "Acessos", value.name = "NR", variable.name = "Plant")
DataQuan19F <- DataQuan19 %>% melt(id.vars = "Acessos", value.name = "FRY", variable.name = "Plant")

DataQuan11F1 <- DataQuan11F %>% group_by(Acessos) %>% summarise(ComprLobulo = mean(ComprLobulo, na.rm = T))
DataQuan12F1 <- DataQuan12F %>% group_by(Acessos) %>% summarise(LargLobulo = mean(LargLobulo, na.rm = T))
DataQuan13F1 <- DataQuan13F %>% group_by(Acessos) %>% summarise(RelComprLar = mean(RelComprLar, na.rm = T))
DataQuan14F1 <- DataQuan14F %>% group_by(Acessos) %>% summarise(ComprPeciolo = mean(ComprPeciolo, na.rm = T))
DataQuan15F1 <- DataQuan15F %>% group_by(Acessos) %>% summarise(EspEntreCasca = mean(EspessuraEntreCasca, na.rm = T))
DataQuan16F1 <- DataQuan16F %>% group_by(Acessos) %>% summarise(ComprMedRzs = mean(ComprMedRzs, na.rm = T))
DataQuan17F1 <- DataQuan17F %>% group_by(Acessos) %>% summarise(DiamMedRzs = mean(DiamMedRzs, na.rm = T))
DataQuan18F1 <- DataQuan18F %>% group_by(Acessos) %>% summarise(NRPl = mean(NR, na.rm = T))
DataQuan19F1 <- DataQuan19F %>% group_by(Acessos) %>% summarise(FRYPl = mean(FRY, na.rm = T))

DataQuan1F <- DataQuan11F1 %>% inner_join(DataQuan12F1, by = "Acessos") %>%
  inner_join(DataQuan13F1, by = "Acessos") %>% inner_join(DataQuan14F1, by = "Acessos") %>%
  inner_join(DataQuan15F1, by = "Acessos") %>% inner_join(DataQuan16F1, by = "Acessos") %>%
  inner_join(DataQuan17F1, by = "Acessos") %>% inner_join(DataQuan18F1, by = "Acessos") %>%
  inner_join(DataQuan19F1, by = "Acessos")

DataQuan <- DataQuan1F %>% inner_join(DataQuan2, by = "Acessos")

DataQuan[DataQuan=="NaN"] <- NA

Shannon.Weaver.QT(DataQuanTeste$ComprLobulo)

Res <- DataQuan %>% dplyr::select(-Acessos) %>% apply(FUN = Shannon.Weaver.QT, 2)

write.csv2(Res, file = "IndicesShannonWeaverQuanCaroline.txt", quote = F)

DataTot <- DataQual %>% left_join(DataQuan, by = "Acessos")

write.table(DataTot, "DadosTodosCaroline.CSV", row.names = F, quote = F, sep = ";", dec = ".")

