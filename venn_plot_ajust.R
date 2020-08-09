
setwd('~/Documentos/')
length(arteria_0.05$gene)
length(arterial_pulmonar_0.05$gene)
length(pulmonar_0.05$gene)

library(grid)
library(futile.logger)
library(VennDiagram)
library(scales)

venn.diagram(
  list(arteria_0.05$gene, arterial_pulmonar_0.05$gene, pulmonar_0.05$gene),
  category.names = c("Hipertensión arterial (1215)", "Hipertensión arterial pulmonar (138)", "Hipertensión pulmonar (127)"),
  filename = "interaction_metavolcano.png",
  output = TRUE,
  imagetype = "png",
  height = 480,  # dimensiones de la imagen
  width = 480,   # dimensiones de la imagen
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#4b54cd", '#c49e2d', '#39c895'),
  fill = c(alpha("#4b54cd",0.3), alpha('#c49e2d',0.3), alpha('#39c895',0.3)),
  cex = 0.5,
  frontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-20, 18, 180),
  cat.dist = c(0.055, 0.055, 0.03),
  cat.fontfamily = "sans",
  cat.col = c("#4b54cd", '#c49e2d', '#39c895'),
  rotation = 1,
  print.mode = "raw",
  main = "Genes perturbados consistentemente en hipertensión esencial",
  main.cex = 0.3)

diferencias_1 = intersect(arteria_0.05$gene, arterial_pulmonar_0.05$gene)    
diferencias_2 = intersect(diferencias_1, pulmonar_0.05$gene)
# "hsa-miR-188-5p"  "hsa-miR-548b-3p" "hsa-miR-486-5p"