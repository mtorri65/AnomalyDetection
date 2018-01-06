setwd("C:\\Projects\\Analytics\\AnomalyDetection\\Client\\Data\\Output\\Analysis")

locations = list.dirs(path = ".")

for(location in locations)
{
  if(location == ".")
  {
    next()
  }
  location = gsub("/", "\\\\", location)
  timeEvolutions = list.files(path = paste(".\\", location, sep = ""), pattern = "^TimeEvolution")  
  for(timeEvolution in timeEvolutions)
  {
    fileName = paste(location, "\\", timeEvolution, sep = "")
    numFields = max(count.fields(fileName, sep = ","))
    dat <- read.table(file=fileName, sep=",",col.names = paste0("V",seq_len(numFields)), fill = TRUE)
    plot(dat$V1, dat$V2, pch=19, type="o", col="red",xlab="Day", ylab="Deviation")
    lines(dat$V1, dat$V3, pch=19, type="o", col="black")
    title(timeEvolution)
    linearFit = lm(dat$V3 ~ dat$V1, data = dat)
    abline(linearFit)
    dev.print(jpeg, paste(location, "\\timeEvolution_", substr(location, 3, nchar(location) - 4), ".jpg", sep = ""), width = 2000, height = 1200)
#        jpeg(paste(location, "\\timeEvolution_", substr(location, 3, nchar(location) - 4), "_", day, ".jpg", sep = ""), width = 2000, height = 1200)
    dev.off   

    readline(prompt="Press [enter] to continue")
  }
}  
