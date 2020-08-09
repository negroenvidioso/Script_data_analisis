library(ggplot2)
# Create a basic bar
pie = ggplot(hp, aes(x="", y=hp$V2, fill=hp$V1))+
  geom_bar(stat="identity", width=1)
# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value), "%")), position = position_stack(vjust = 0.5),size = 4)
# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "Porcentaje de insidencia de \n hipertension")
# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "#666666", size = 30))
pie