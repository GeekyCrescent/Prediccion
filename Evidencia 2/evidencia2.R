# CORNONAVIRUS DE ANIMALES

library('seqinr')

# -------------------------------------------------- CONSEGUIMOS LOS DATOS

# Bat - Desmodus Rotundus
bat <- read.fasta("Fastas/bat.fasta")
# Pacific Salmon - Oncorhynchus tshawytscha
pacific_salmon <- read.fasta("Fastas/pacific_salmon.fasta")
# Duck - Anatidae
duck <- read.fasta("Fastas/duck.fasta")
# Canada Goose - Branta canadensis
goose <- read.fasta("Fastas/goose.fasta")
# Shrew - Sorex araneus
shrew <- read.fasta("Fastas/shrew.fasta")
# Rodent - Myodes rufocanus
rodent <- read.fasta("Fastas/rodent.fasta")
# Porcine - Sus scrofa
porcine <- read.fasta("Fastas/porcine.fasta")
# Bulbul - Pycnonotus jocosus
bulbul <- read.fasta("Fastas/bulbul.fasta")
# Rat - Rattus norvegicus
rat <- read.fasta("Fastas/rat.fasta")
# Ferret - Ferret 
ferret <- read.fasta("Fastas/ferret.fasta")
# Swine - Sus scrofa
swine <- read.fasta("Fastas/swine.fasta")
# Camel - Camelus
camel <- read.fasta("Fastas/camel.fasta")
# Mink - Neovison vison
mink <- read.fasta("Fastas/mink.fasta")
# Rabbit- Oryctolagus cuniculus
rabbit <- read.fasta("Fastas/rabbit.fasta")
# Night heron - Ardeidae
night_heron <- read.fasta("Fastas/night_heron.fasta")
# Wigeon - Mareca
wigeon <- read.fasta("Fastas/wigeon.fasta")
# Thrush - Turdus hortulorum
thrush <- read.fasta("Fastas/thrush.fasta")
# Turkey - Meleagris gallopavo
turkey <- read.fasta("Fastas/turkey.fasta")
# Beluga Whale - Delphinapterus leucas
beluga_whale <- read.fasta("Fastas/beluga_whale.fasta")
# Sparrow - Passeridae
sparrow <- read.fasta("Fastas/sparrow.fasta")
# Human - Homosapiens
human <- read.fasta("Fastas/human.fasta")


# -------------------------------------------------- VECTORES PARA LOS FORS

# Vector que incluye todas las secuencias
todas_secuencias <- c(
  bat, pacific_salmon, duck, goose, shrew, rodent, porcine, 
  bulbul, rat, ferret, swine, camel, mink, rabbit, night_heron, 
  wigeon, thrush, turkey, beluga_whale, sparrow, human
)

# Lista de nombres de las secuencias
nombres_secuencias <- c(
  "bat", "pacific_salmon", "duck", "goose", "shrew", "rodent", "porcine", 
  "bulbul", "rat", "ferret", "swine", "camel", "mink", "rabbit", "night_heron", 
  "wigeon", "thrush", "turkey", "beluga_whale", "sparrow", "human"
)


# --------------------------------------------------  LONGITUDES DE SECUENCIAS

# Loop para imprimir los nombres de las secuencias
for (i in 1 : length(todas_secuencias)) {
  mensaje <- sprintf("%s : %d", nombres_secuencias[i], length(todas_secuencias[i][[1]]))
  print(mensaje)
}

# -------------------------------------------------- Conseguimos las bases


# Crear las variables de conteo para cada secuencia
bat_c <- count(bat[[1]], 1)
pacific_salmon_c <- count(pacific_salmon[[1]], 1)
duck_c <- count(duck[[1]], 1)
goose_c <- count(goose[[1]], 1)
shrew_c <- count(shrew[[1]], 1)
rodent_c <- count(rodent[[1]], 1)
porcine_c <- count(porcine[[1]], 1)
bulbul_c <- count(bulbul[[1]], 1)
rat_c <- count(rat[[1]], 1)
ferret_c <- count(ferret[[1]], 1)
swine_c <- count(swine[[1]], 1)
camel_c <- count(camel[[1]], 1)
mink_c <- count(mink[[1]], 1)
rabbit_c <- count(rabbit[[1]], 1)
night_heron_c <- count(night_heron[[1]], 1)
wigeon_c <- count(wigeon[[1]], 1)
thrush_c <- count(thrush[[1]], 1)
turkey_c <- count(turkey[[1]], 1)
beluga_whale_c <- count(beluga_whale[[1]], 1)
sparrow_c <- count(sparrow[[1]], 1)
human_c <- count(human[[1]], 1)

# Crear el vector de conteos combinando cada variable individualmente
conteos_vector <- c(
  bat_c, pacific_salmon_c, duck_c, goose_c, shrew_c, rodent_c, porcine_c, 
  bulbul_c, rat_c, ferret_c, swine_c, camel_c, mink_c, rabbit_c, night_heron_c, 
  wigeon_c, thrush_c, turkey_c, beluga_whale_c, sparrow_c, human_c
)

# -------------------------------------------------- GRÁFICA QUE MUESTRAS LAS BASES


library(ggplot2)


# Crear el data frame combined_counts
combined_counts <- data.frame(
  Variante_Virus = rep(nombres_secuencias, each = 21/5),  
  Nucleotidos = rep(c('a', 'c', 'g', 't'), times = 21), 
  nucleotide = rep(LETTERS[1:4], each = 21),            
  Cantidad = conteos_vector                            
)


# Gráfico de barras con ggplot2 y labels del eje x orientados horizontalmente
ggplot(combined_counts, aes(fill = Nucleotidos, y = Cantidad, x = Variante_Virus)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# ---------------------------------------------------- Hacer el árbol filogenético


# Importamos librerías necesarias y creamos un vector con todos los accessions
library(Biostrings)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

accession_virus <- c('NC_055953', 'NC_076697', 'NC_048214', 'NC_046965', 'NC_046955', 
                     'NC_046954', 'NC_039208', 'NC_011547', 'NC_032730', 'NC_030292',
                     'NC_028806', 'NC_028752', 'NC_023760', 'NC_017083', 'NC_016994',
                     'NC_016995', 'NC_011549', 'NC_010800', 'NC_010646', 'NC_016992',
                     'NC_045512')

# Creamos un objeto que lee las secuencias y lo anotamos en un archivo

virus_sequences <- read.GenBank(accession_virus)

write.dna(virus_sequences, file='virus_seqs.fasta', format='fasta', append=
            FALSE, nbcol=6, colsep=' ', colw=10)

# Cargamos las secuencias

virus_seq_not_align <- readDNAStringSet('virus_seqs.fasta', format='fasta')

virus_seq_not_align <- OrientNucleotides(virus_seq_not_align)

virus_seq_algin <- AlignSeqs(virus_seq_not_align)

# Anotamos la secuencia en otro archivo

writeXStringSet(virus_seq_algin, file='virus_seq_align.fasta')

# Obtenemos el nuevo archivo y creamos una matriz de distancia

virus_aligned <- read.alignment('virus_seq_align.fasta', format='fasta')

matriz_distancia <- dist.alignment(virus_aligned, matrix='similarity')

# Creación del árbol filogenético

virus_tree <- nj(matriz_distancia)
class(virus_tree)

virus_tree <- ladderize(virus_tree)

# Plot del árbol

plot(virus_tree, cex=0.6)
title('COVID')

# Recreamos el árbol con los nombres de los animales

text.string<-
"((((((Mink, Ferret) ,Swine) ,Camel, Bat), Rat), Shrew),((Rabbit, Rodent) ,Human), ((((Turkey, Duck), Rodent), Beluga Whale), (((Bulbul, Thrush), ((Porcine, Sparrow), Night Heron)), Wigeon), (Pacific Salmon)));"
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE,edge.width=2)


