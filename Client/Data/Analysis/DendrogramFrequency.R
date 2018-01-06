setwd("C:\\Projects\\Analytics\\Anomaly\\Client\\Data\\Output\\Analysis")
library(dendextend)

locations = list.dirs(path = ".")

for(location in locations)
{
  if(location == ".")
  {
    next()
  }
  location = gsub("/", "\\\\", location)
  
  distanceMatrices = list.files(path = paste(".\\", location, sep = ""), pattern = "^DistanceFrequencyMatrix")
  
  for(distanceMatrix in distanceMatrices)
  {
    day = substr(distanceMatrix, nchar("DistanceFrequencyMatrix") + 1, nchar(distanceMatrix) - 4)
    fileName = paste(location, "\\DistanceFrequencyMatrix", day, ".csv", sep = "")
    tb <- read.csv(fileName, header=TRUE);
    
    dend <- tb %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k=1) %>% set("branches_lwd", 3.5) %>% set("labels_cex", .9) %>% set("leaves_pch", 19) %>% set("leaves_col", c("blue", "red"))
    jpeg(paste(location, "\\frequency_", substr(location, 3, nchar(location) - 4), "_", day, ".jpg", sep = ""), width = 2000, height = 1200)
    plot(dend, main = paste(substr(location, 3, nchar(location) - 4), "   Day", day), cex.main = 2.5)
    dev.off()
    plot(dend, main = paste(substr(location, 3, nchar(location) - 4), "   Day", day))
  }
}

