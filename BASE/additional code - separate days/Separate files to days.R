

# Script for separating multiple days into separate diel periods and saving in individual csv files


#----------------------------------

# UPDATE DIRECTORY ONLY:
# set location of "additional code - separate days" folder

directory <- "C:/Desktop/Analysis/BASE/additional code - separate days"  # example directory

#----------------------------------


fname.list<-list.files(file.path(directory,"/input"))

for (fname in fname.list)
    {
    
    data<-read.csv(file.path(directory,"/input",fname), sep=",", header=T)
    
    unique.dates<-unique(data$Date)
    unique.dates<-unique.dates[1:(length(unique.dates))] 
    
        for (day in unique.dates) 
        {
    
        min.row<-min(as.numeric(rownames(data[data$Date==day,])))
        max.row<-max(as.numeric(rownames(data[data$Date==day,])))
    
        subset.data<-data[min.row:(max.row),]
        
        new.fname<-paste(substr(fname,1,nchar(fname)-4),"_", day, ".csv", sep="")
        write.csv(subset.data, file.path(directory,"/output",new.fname), row.names=F)
                      
        }
    
    }

