# Tenemos 6 muestras del grupo 1, 16 del grup 2 y 27 del grupo 3
# Leemos el archivo "Asignacion_de_Muestras_A_Grupos.csv" que hemos creado
# a partir del archivo excel. Yo lo he creado separado por ";"

setwd("/home/alex/Dropbox/Classes/UOC/AnalisisDeMicroArrays/Curs-2012/Exercicis i PECS 2012/PEC-1")
muestras <- read.csv2("Asignacion_de_muestras_a_grupos.csv", head=T)
misMuestras <- as.character (muestras$Sample.ID)
paraAnalisis <- c(sample(misMuestras[1:6], 5), 
                  sample(misMuestras[7:20], 5),
                  sample(misMuestras[23:48], 5))