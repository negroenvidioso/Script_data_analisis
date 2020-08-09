x <- c(mirna_humanos$ID, mirna_humanos$Mature1_ID, mirna_humanos$Mature2_ID)  # vector con los ID de los miRNA
y <- c(mirna_humanos$Accession, mirna_humanos$Mature1_Acc, mirna_humanos$Mature2_Acc)  # vector con los numeros de accsion

df <- data.frame(x,y)  # data frame con los id y los numeros de accesion

sapply(df, function(df) sum(is.na(df)))  # contar los valores vacios entre las columnas

delete.na <- function(df, n=0) {
  df[rowSums(is.na(df)) <= n,]
}
delete.na(df)  # borra todos las filas con los valores nulos

datos <- na.omit(df)

sapply(datos, function(datos) sum(is.na(datos)))

write.table(datos, file = "miRNA_ID_acc.txt")
