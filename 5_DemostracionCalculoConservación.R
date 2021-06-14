#El analisis se realiza en 4 pasos
#1) Carga de datos de IGIs
#2) ALINEAMIENTOS
#3) BUSQUEDA DE PQSs en las secuencias consenso
#4) Generacion de LOGOS de los PQSs de interes
################################################

#RETOQUES PRELIMINARES:####
#Creo un vector con los nombres resumidos de los genomas, que me ayudara a trabajar
#ATENCION: Estos nombres deben coincidir entre los archivos a utilizar (de la carpeta outputs_1)
# y las columnas de la tabla de ortologos
genomesNames = c("Linfantum","Ldonovani","Lmajor","TcruziNE" , "Tbrucei","Tevansi", "TcruziE")
#Abro la tabla de salida de Roary con los nombres de las proteinas ortologas
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff")
ortologos = read.csv("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/Resultado script 5/Galaxy75-[Roary_on_data_51,_data_50,_and_others_Gene_Presence_Absence].csv",stringsAsFactors = F)
#Reviso la distribucion de secuencias por grupo de ortologos
table(ortologos$No..sequences)
#Grafico la distribucion de secuencias por grupo de ortologos
plot(table(ortologos$No..sequences), main="ortologos, blast70, MCL0.5")
####
#Renombro las columnas con el vector de nombres creado al principio
colnames(ortologos)[15:21]
#Me aseguro que el vector de nombres esta en el orden correcto (mismo que las columnas)
#(Esto debe ser verificado manualmente)
genomesNamesReordered=genomesNames
#Renombro las columnas
colnames(ortologos)[15:21]=genomesNamesReordered
#Obtengo los ortologos con representantes en los x genomas y reduzco la tabla
x=3
ortologos = ortologos[ortologos$No..sequences==x,c(15:21)]

#1)Carga de datos de IGIs ############################################################################
#Cargo la lista de archivos de la carpeta 1 (Primero los csv y luego los fasta)
IGIsCsvs=list.files("./Resultado script 1", pattern = "*.csv")
IGIsFasta=list.files("./Resultado script 1", pattern = "*.fasta")
#Leo las secuencias IGIs completas y sus anotaciones (me interesa conocer las direcciones de las hebras)list.files("./Outputs_1")
#Cargo todo en una lista
IGIsData=list()
IGIsData$csv=list()
IGIsData$fasta=list()
setwd("./Resultado script 1")
for (i in 1:7){
  #Tomo los nombres individuales de los genomas
  genomeName=genomesNamesReordered[i]
  #Obtengo los paths de cada archivo (de a un genoma por vez)
  pathCsv=IGIsCsvs[grep(genomeName,IGIsCsvs)]
  pathFasta=IGIsFasta[grep(genomeName,IGIsFasta)]
  #Abro y guardo en la lista los archivos csv y fasta de cada genoma
  print(paste0("Cargando datos de ", genomeName, "..."))
  IGIsData$csv[[genomeName]]=read.csv(pathCsv, stringsAsFactors = FALSE)
  IGIsData$fasta[[genomeName]]=Biostrings::readDNAStringSet(pathFasta)
  
}

#2)ALINEAMIENTOS###########################################################################
# Para cada grupo de ortologos, tomo las secuencias de cada genoma, 
# las agrupo en un mismo objeto DNAStringSet y aplico los diferentes
# algoritmos de alineamiento
library(msa)
alignments=list()
#alignments$clustalOmega=list()
#alignments$clustalW=list()
alignments$muscle=list()
ortologos=ortologos[,!ortologos[1,]==""]
for(i in 1:nrow(ortologos)){#bucle 1 para ir analizando cada grupos e ortologos
  multiSeq=Biostrings::DNAStringSet()
  k=0;#para controlar cuantas secuencias se invierten
  for (j in 1:ncol(ortologos)){#bucle 2, toma los IGIs IDs de cada columna de la tabla de ortologos
    #Tomo los nombres individuales de los genomas
    # genomeName=genomesNamesReordered[j]
    genomeName=colnames(ortologos)[j]
    #Tomo el ID del IGI para el grupo de ortologos (j) y genoma actual (i)
    IGIname=ortologos[i,colnames(ortologos)==genomeName]
    #Extraigo la secuencia
    seq=IGIsData$fasta[[genomeName]][IGIname]
    #IMPORTANTE: hago correccion de hebra (reverseComplement para hebras negativas)
    strand=as.character(IGIsData$csv[[genomeName]][IGIsData$csv[[genomeName]]$ID==IGIname,]$strand)
    if(strand=="-"){
      seq=Biostrings::reverseComplement(seq)
      k=k+1
    }
    multiSeq=c(multiSeq, seq)
  }
  print(paste0("Alineando ortologos numero ",i,". (Se han invertido ",k," secuencias antes del alineamiento)..."))
  print(multiSeq)
  #alignments$clustalOmega[[i]] = msa(multiSeq,"ClustalOmega")
  #alignments$clustalW[[i]] = msa(multiSeq,"ClustalW")
  ti = Sys.time()
  alignments$muscle[[i]] = msa(multiSeq,"Muscle")
  te=Sys.time()
  print(paste0("Alineamiento numero ",i,"exitoso. Se ha demorado: ", te-ti, " min"))
}
#Guardo los alineamientos como objeto R###
# save(alignments,file="../Resultado script 6/alignments.Rdata")
load("../Resultado script 6/alignments.Rdata")
###
list.files()
setwd("Resultado script 6")
#3) BUSQUEDA DE PQSs en las secuencias consenso##########################################################
#Obtengo las secuencias consenso de los alineamientos y hago
#una busqueda de PQSs en dichas secuencias.
#Adicionalmente obtengo los valores de conservacion de cada nt
#(Variarcion segun como se construya la matriz de sustitucion "mat")
mat <- nucleotideSubstitutionMatrix(4, -1)
mat <- cbind(rbind(mat, "-"=-8), "-"=-8)
library(pqsfinder)
#Creo df donde almacenaron toda la info de los pqs
PQSfullData=data.frame(ortologo="",pqsID="",start="",width="",score="",
                       strand="", pqs="", conservationMeanIn="",
                       conservationMeanOut="", conservationDiff="",
                       stringsAsFactors = F)
#Bucleo para cada alineamiento y calculo los PQSs
for(i in 1:length(alignments$muscle)){
  print(paste0("Calculando pqs de alineamiento ",i))
  alineamiento = alignments$muscle[[i]]
  consenso = msaConsensusSequence(alineamiento)
  #Obtengo los valores de conservacion de cada nucleotido
  conservation = msaConservationScore(alineamiento, mat, gapVsGap=0)
  #Calculo pqs a partir de consenso
  consenso=gsub("\\?","N",consenso)
  consenso=DNAString(consenso)
  consensusPQSs=pqsfinder(consenso)
  print(paste0("Se han encontrado ", length(consensusPQSs), " cuadruples en la secuencia consenso ",i))
  if(length(consensusPQSs)>0){
    pqsData=data.frame()
    for(j in 1:length(consensusPQSs)){
      start=start(consensusPQSs[j])
      end=start(consensusPQSs[j])+width(consensusPQSs[j])-1
      #Obtengo las conservaciones medias para el pqs y su entorno inmediato (+-25nt)
      conservationMeanIn=mean(conservation[start:end])
      conservationMeanOut=mean(c(conservation[(start-25):start],conservation[end:(end+25)]))
      conservationDiff=conservationMeanIn-conservationMeanOut
      #Cargo todo en un vector y lo agrego a la tabla general
      x=c(ortologo=i, pqsID=j, start=start, width=width(consensusPQSs[j]),
          score=score(consensusPQSs[j]), strand=strand(consensusPQSs[j]),
          pqs=as.character(subseq(consensusPQSs[j]@subject, start, end)),
          conservationMeanIn=conservationMeanIn, conservationMeanOut=conservationMeanOut,
          conservationDiff=conservationDiff)
      PQSfullData=rbind(PQSfullData,x)
    }
  }
}
PQSfullData=PQSfullData[-(PQSfullData$ortologo==""),]
#Guardo en formato CSV
write.csv(PQSfullData,"../Resultado script 6/PQSs_fullData2.csv", row.names = FALSE)

#4) Generacion de LOGOS de los PQSs de interes############################################################
library(ggseqlogo)
#Grafico el entorno los PQSs consenso
for(i in 1:nrow(PQSfullData)){
  alignNumber=as.integer(PQSfullData[i,]$ortologo)
  pqsID=as.integer(PQSfullData[i,]$pqsID)
  start=as.integer(PQSfullData[i,]$start)
  end=start+as.integer(PQSfullData[i,]$width)-1
  alineamiento = alignments$muscle[[alignNumber]]
  subseqs=c()
  for(j in 1:5){
    temp=as.character(subseq(alineamiento@unmasked[j], start-10,end+10))
    subseqs=c(subseqs,temp)
  }
  png(filename=paste0("../Resultado script 6/PQSlogo_ortoGroup_",alignNumber,"_pqsID_",pqsID,".png"),width = 1000, height = 400)
  print(ggseqlogo(subseqs, method = 'bits', seq_type='DNA'))
  dev.off()
}

#OPCIONAL
# Grafico el entorno los PQSs consenso (solo en genomas de Leishmanias)
for(i in 1:nrow(PQSfullData)){
  alignNumber=as.integer(PQSfullData[i,]$ortologo)
  start=as.integer(PQSfullData[i,]$start)
  end=start+as.integer(PQSfullData[i,]$width)-1
  alineamiento = alignments$muscle[[alignNumber]]
  subseqs=c()
  for(j in 4:5){#Las 3 ultimas secuencias del alineamiento corresponden a los Leishmania
    temp=as.character(subseq(alineamiento@unmasked[j], start-10,end+10))
    subseqs=c(subseqs,temp)
  }
  png(filename=paste0("../Resultado script 6/PQSlogo_ortoGroup_",alignNumber,"_pqs_",i,".png"),width = 1000, height = 400)
  print(ggseqlogo(subseqs, method = 'bits', seq_type='DNA'))
  dev.off()
}
