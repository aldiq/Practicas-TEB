#Adaptacion de GFF para que funcionen como input de Roary
#Incluye: 
#1)reduccion de genes totales
#2)reduccion de secuencias totales a fin de disminuir el tama�o de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/L. infantum JPCM5")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_LinfantumJPCM5.gff")
igisFiltrados = read.csv("IGIs_filtrados_LinfantumJPCM5.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_LinfantumJPCM5_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interes)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guion bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posicion
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guion bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Funcion accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guion bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}








setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LinfantumJPCM5/PQSs")
list.files()

#######################################################################
#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/L.donovani BPK282A1")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_LdonovaniBPK282A1.gff")
igisFiltrados = read.csv("IGIs_filtrados_LdonovaniBPK282A1.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_LdonovaniBPK282A1_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LdonovaniBPK282A1/PQSs")
list.files()

############################################################################

#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/Leishmania major Friedlin")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_LmajorFriedlin.gff")
igisFiltrados = read.csv("IGIs_filtrados_LmajorFriedlin.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_LmajorFriedlin_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/LmajorFriedlin/PQSs")
list.files()

################################################################################
#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/T. brucei gambiense DAL972")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_TbruceigambienseDAL972.gff")
igisFiltrados = read.csv("IGIs_filtrados_TbruceigambienseDAL972.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_TbruceigambienseDAL972_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TbruceigambienseDAL972/PQSs")
list.files()

###############################################################################

#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/T. evansi STIB805")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_TevansiSTIB805.gff")
igisFiltrados = read.csv("IGIs_filtrados_T. evansi STIB805.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_TevansiSTIB805_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/T. evansi STIB805/PQSs")
list.files()
############################################################################

#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/Trypanosoma cruzi CL Brener/Esmeraldo-like")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_TcruziCLBrenerEsmeraldo-like.gff")
igisFiltrados = read.csv("IGIs_filtrados_TcruziCLBrenerEsmeraldo.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_TcruziCLBrenerEsmeraldo-like_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TcruziCLBrenerEsmeraldo/PQSs")
list.files()
#################################################################################

#Adaptaciónd e GFF para que funcionen como input de Roary
#Incluye: 
#1)reducción de genes totales
#2)reducción de secuencias totales a fin de disminuir el tamaño de los archivos
setwd("F:/Aldana/PRACTICANATO/PROYECTO/archivos fasta gff/Trypanosoma cruzi CL Brener/Non Esmeraldo-like")
list.files()
annotations = rtracklayer::readGFF("TriTrypDB-48_TcruziCLBrenerNon-Esmeraldo-like.gff")
igisFiltrados = read.csv("IGIs_filtrados_TcruziCLBrenerNon-Esmeraldo-like.csv")
genomeLeishmania = Biostrings::readDNAStringSet("TriTrypDB-48_TcruziCLBrenerNon-Esmeraldo-like_Genome.fasta")

#Obtengo gff reducido (filtrado con los IGIs de interés)
filteredAnnotations=annotations[annotations$ID%in%igisFiltrados$x,]
filteredAnnotations$type="CDS"
filteredAnnotations$phase=0
# system("rm Prueba2.gff")
# system("touch Prueba2.gff")
write("##gff-version 3",file="Prueba2.gff")
# write("##gff-version 3",file="Prueba2.gff",append=TRUE)
for(i in 1:length(genomeLeishmania)){
  #for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo una longitud aproximada
  #length = max(filteredAnnotations[filteredAnnotations$seqid==name,]$end)
  length = length(genomeLeishmania[[grep(name, names(genomeLeishmania))]])
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name = gsub("\\.","_",name)
  line = paste0("##sequence-region ",name," 1 ",length)
  write(line,file="Prueba2.gff",append=TRUE)
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){#Para prueba con los dos primeros cromosomas
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #obtengo las anotaciones de a un cromosoma y las ordeno por posición
  filteredAnnotationsSubset = filteredAnnotations[filteredAnnotations$seqid==name,]
  filteredAnnotationsSubset = filteredAnnotationsSubset[order(filteredAnnotationsSubset$start),]
  #Agrego las lineas al archivo prueba
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  filteredAnnotationsSubset$seqid = gsub("\\.","_",filteredAnnotationsSubset$seqid)
  rtracklayer::export.gff3(filteredAnnotationsSubset, "Prueba2.gff",  append=TRUE)
  
}
write("##FASTA",file="Prueba2.gff",append=TRUE)
#Agrego las secuencias fasta
#Función accesoria para recortar las secuencias en lineas de 60 caracteres (por las dudas)
tensiSplit <- function(string,size) {
  stringr::str_extract_all(string, paste0('.{1,',size,'}'))
}
for(i in 1:length(genomeLeishmania)){
  # for(i in 1:2){
  seqs = levels(filteredAnnotations$seqid)
  seqs = seqs[order(seqs)]
  name = seqs[i]
  #Cualquier punto en el nombre de la secuencia debe ser reemplazado con guin bajo
  name_modified = gsub("\\.","_",name)
  write(paste0(">",name_modified),file="Prueba2.gff",append=TRUE)
  #Recorto la secuencia y bucleo
  seqLines = tensiSplit(as.character(genomeLeishmania[[grep(name, names(genomeLeishmania))]]),60)[[1]]
  for (j in 1:length(seqLines)){
    line = seqLines[j]
    write(line,file="Prueba2.gff",append=TRUE)
  }
  print(i)
}
setwd("C:/Users/Daniel/Dropbox/Practicas_Aldana/tablas y secuencias obtenidas/TcruziCLBrenerNon-Esmeraldo-like/PQSs")
list.files()
