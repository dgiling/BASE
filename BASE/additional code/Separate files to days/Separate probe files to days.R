

# Script for separating multiple days into separate diel periods and saving in individual csv files


#----------------------------------

# UPDATE DIRECTORIES ONLY:

input.directory<-"[your directory]/additional code/Separate files to days/input"
output.directory <- "[your directory]/additional code/Separate files to days/output"


#----------------------------------



setwd(input.directory)
fname.list<-list.files(input.directory)


for (fname in fname.list)
{

setwd(input.directory)
data<-read.csv(fname, sep=",", header=T)

unique.dates<-unique(data$Date)
unique.dates<-unique.dates[1:(length(unique.dates))] 

    for (day in unique.dates) 
    {

    min.row<-min(as.numeric(rownames(data[data$Date==day,])))
    max.row<-max(as.numeric(rownames(data[data$Date==day,])))

    subset.data<-data[min.row:(max.row+1),]
    
              
                  new.fname<-paste(substr(fname,1,nchar(fname)-4),"_", day, ".csv", sep="")
                  
                  setwd(output.directory)
                  
                  write.csv(subset.data, new.fname, row.names=F)
                  
       }

}

