#CARGO MIS FUNCIONES ARTESANALES
source("F:/Aldana/PRACTICANATO/PROYECTO/funciones/igiCalculator.R")
source("F:/Aldana/PRACTICANATO/PROYECTO/funciones/igiExtractor.R")
##################################################################################
# Demostración con L.major
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/Trypanosoma cruzi CL Brener/Non Esmeraldo-like")
#Calculo y guardo las IGIs del primer organismos
igisDF = igiCalculator(gffFile="TriTrypDB-48_TcruziCLBrenerNon-Esmeraldo-like.gff")
#Puedo guardar la tabla en formato csv (o el que considere más adecuado)
write.csv(igisDF, "IGIs_TcruziCLBrenerNon-Esmeraldo-like.csv", row.names = FALSE)
#xtraigo y guardo las secuencias de las IGIS calculadas previamente
igisSequences = igiExtractor(fastaFile="TriTrypDB-48_TcruziCLBrenerNon-Esmeraldo-like_Genome.fasta", positions=igisDF)
#Puedo guardar las secuencias en formato fasta
Biostrings::writeXStringSet(igisSequences, filepath="IGIs_TcruziCLBrenerNon-Esmeraldo-like_Genome.fasta")
#####################################################################################
