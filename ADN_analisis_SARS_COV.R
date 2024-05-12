## Libraries used

# BiocManager::install("ggtree")
# BiocManager::install("DECIPHER")
# install.packages("viridis")

suppressMessages(library(seqinr))
suppressMessages(library(adegenet))
suppressMessages(library(ape))
suppressMessages(library(ggtree))
suppressMessages(library(DECIPHER))
suppressMessages(library(viridis))
suppressMessages(library(ggplot2))

# Funciones: Regresa la cantidad de cada nucleótido en la secuencai de ADN
cantidad_nucleotidos <- function(secuencia)
{
  library(stringr)
  cant_A <- sum(str_count(secuencia, "A"))
  cant_T <- sum(str_count(secuencia, "T"))
  cant_C <- sum(str_count(secuencia, "C"))
  cant_G <- sum(str_count(secuencia, "G"))
  total <- nchar(secuencia)
  
  # porcentaje_A <- cant_A / total * 100
  # porcentaje_T <- cant_T / total * 100
  # porcentaje_C <- cant_C / total * 100
  # porcentaje_G <- cant_G / total * 100
  # 
  return(c(cant_A,cant_T,cant_C,cant_G)) 
}

# Leer el archivo fasta
covid <- readDNAStringSet("covid_file.fasta", format = "fasta")

covid

## Tabla 1: longitud de secuencias
suppressMessages(covid <- OrientNucleotides(covid))
suppressMessages(alineada <- AlignSeqs(covid))
writeXStringSet(alineada, 
                # file="C:/Users/Santiago/Documents/R_Codes/Etapa2/covid.fasta")
                file="covid.fasta")

secuencias <- as.character(covid)
long <- sapply(secuencias, nchar)
long <- long-21000

secuencias_df <- data.frame(
  Secuencias = c("China", "Alemania", "Brasil","EUA", "Francia", "India", "Italia", "Korea", "Nigeria",
                 "SA"),
  # Secuencias <- names(long),
  Longitud = long,
  Continentes = c("Asia", "Europa", "America", "America", "Europa", "Asia", "Europa", "Asia", "Africa", "Africa")
)

secuencia_df <- secuencias_df[order(secuencias_df$Continentes), ]

# Primera tabla, comparación en la longitud de las secuencias. 
tabla <- ggplot(secuencia_df, aes(x = reorder(Secuencias, order(Continentes)), y = Longitud, fill = Continentes, label = Secuencias)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "País", y = "´Número de caracteres (+21000)", title = "Longitud de genoma")
tabla

## Tabla 2: número de baes en cada secuencia de ADN
covid_seq <- as.character(covid)
nucleotidos_df <- data.frame()
for(i in 1:length(covid_seq))
{
  nuc <- cantidad_nucleotidos(covid_seq[[i]])
  nucleotidos <- data.frame("Location" = names(covid_seq[i]),
                            "Nucleotidos" = c("A", "T","C","G"),
                            "Bases" = c(nuc[1],nuc[2],nuc[3],nuc[4]))
  nucleotidos_df <- rbind(nucleotidos_df, nucleotidos)
}
num_bases <- ggplot(nucleotidos_df, aes(x = Bases, y = Location, fill = Nucleotidos)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis_d()

num_bases

## Heat-Map: Similitudes entre las diferentes variantes. 
dna <- read.dna("covid.fasta", format = "fasta")

D <- dist.dna(dna, model="TN93") #TN93 es un modelo en base a la probabilidad de cada nucleótido de mutar en uno diferente, teniendo diferentes valores conseguidos en base a ecuaciones que toman en cuenta la transición y transposición de los mismos.  

similitud <- as.matrix(D)

colores <- colorRampPalette(c("cadetblue1", "midnightblue"))
heatmap(similitud, Rowv = NA, Colv = NA, col =colores(50), margins =c(10, 10))

# Árbol filogenético
tree <- nj(D)

tree <- ladderize(tree)

myBoots <- boot.phylo(tree, dna, function(e) root(nj(dist.dna(e, model = "TN93")),1))
myBoots <- c(rep(0, times=10), myBoots)


ggtree(tree) + 
  geom_tiplab(hjust = -0.3, size=4, align = TRUE)+
  xlim(0,.005)

tip_color_var = viridis(18)
ggtree(tree) +
  geom_nodepoint() +
  geom_nodelab(aes(label = myBoots), repel=TRUE, align=TRUE) +
  geom_tiplab(aes(color = tip_color_var), hjust = -0.3, size=4, align = TRUE) +
  scale_color_manual(values = c("#FF6663", "#FF6663", "#F1C40F", "#6a7fdb", "#F1C40F", "#6a7fdb", "#57e2e5", "#F1C40F", "#57e2e5", "#6a7fdb"), guide = NULL) +
  xlim(0, 0.005)
