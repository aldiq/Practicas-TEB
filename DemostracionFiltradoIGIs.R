#Filtrado de Genoma Ldonovani
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LdonovaniBPK282A1/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_LdonovaniBPK282A1.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_LdonovaniBPK282A1.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_LdonovaniBPK282A1.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)


#################################################################################

#Filtrado de Genoma Lmajor
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LmajorFriedlin/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_LmajorFriedlin.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_LmajorFriedlin.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85#defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_LmajorFriedlin.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)


##########################################################################

#Filtrado de Genoma Linfantum
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LinfantumJPCM5/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_LinfantumJPCM5.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_LinfantumJPCM5.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85 #defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_LinfantumJPCM5.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)

###################################################################

#Filtrado de Genoma T. Evansi
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/T. evansi STIB805/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_T. evansi STIB805.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_T. evansi STIB805.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85 #defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_T. evansi STIB805.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)


###################################################################################

#Filtrado de Genoma Tbrucei Gambiense
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TbruceigambienseDAL972/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_TbruceigambienseDAL972.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_TbruceigambienseDAL972.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85 #defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_TbruceigambienseDAL972.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)

#################################################################################

#Filtrado de Genoma Tbrener Esmeraldo
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TcruziCLBrenerEsmeraldo/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_TcruziCLBrenerEsmeraldo.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_TcruziCLBrenerEsmeraldo.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85 #defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_TcruziCLBrenerEsmeraldo.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)


#######################################################################################

#Filtrado de Genoma T.Brener Non Esmeraldo Like
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TcruziCLBrenerNon-Esmeraldo-like/PQSs")
pqsByIGIs = read.csv("pqsByIGIs_TcruziCLBrenerNon-Esmeraldo-like.csv", stringsAsFactors = FALSE)
pqsPositions = read.csv("pqsPositions_TcruziCLBrenerNon-Esmeraldo-like.csv", stringsAsFactors = FALSE)
#Reducción básica según ocurrencia o no de pqs
pqsByIGIsReducida = pqsByIGIs[pqsByIGIs$pqsNumber>0,]
#reducción fina según el score de PQSs
tresholdScore = 85 #defino punto de corte
tablaTemp = pqsPositions[pqsPositions$score>=tresholdScore,]
secuenciasUnicas = unique(tablaTemp$seq_id)
write.csv(secuenciasUnicas, "IGIs_filtrados_TcruziCLBrenerNon-Esmeraldo-like.csv", row.names = F)
# #reducción fina más recuento de pqs
# x = aggregate(tablaTemp, list("seq_id"), FUN=sum)
