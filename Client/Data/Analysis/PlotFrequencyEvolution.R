setwd("C:\\Projects\\Analytics\\AnomalyDetection\\Client\\Data\\Output\\Analysis")

locations = list.dirs(path = ".")

for(location in locations)
{
  if(location == ".")
  {
    next()
  }
  location = gsub("/", "\\\\", location)
  timeEvolutions = list.files(path = paste(".\\", location, sep = ""), pattern = "^FrequencyEvolution")  
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
    
    readline(prompt="Press [enter] to continue")
  }
}  
