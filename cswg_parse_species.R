download.file("https://docs.google.com/spreadsheets/d/1g5mb9zL5LH32msA8Igsyap9AwR9UwS4NuOiItUdQgZ8/pub?gid=0&single=true&output=csv",
               destfile="output.csv",
               quiet=T)
               
t <- read.csv("output.csv")
file.remove("output.csv")

t[,4:7][is.na(t[,4:7])] <- 0 # SWAP State NA codes with zeros

spp <- unique(as.character(t[,1]))
out <- data.frame()
for(i in 1:length(spp)){    

    focal <- cbind(t[which(grepl(t$Species,pattern=spp[i]))[1],1:3],
                   matrix(colSums(t[grepl(t$Species,pattern=spp[i]),4:7]),nrow=1)
                   )

    names(focal) <- names(t[which(grepl(t$Species,pattern=spp[i]))[1],])
    
    if(length(out) == 0){
      out <- focal
    } else {
      out <- rbind(out,focal)
    }
    
}

write.csv(out,"out.csv",row.names=F)
