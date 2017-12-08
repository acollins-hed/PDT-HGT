#!/usr/bin/Rscript

library("RColorBrewer")

data <- read.table("datahgt.txt",sep=",")
datan <- read.table("datanohgt.txt",sep=",")

if(max(data[,82])>=max(datan[,82]))
    MAX = max(data[,82]) else
    MAX = max(datan[,82])

par(mfrow=c(1,2))

plot(datan[,1],datan[,82],col="red",main="Mean Code Distance Over Time",xlab="Time Steps",ylab="Mean Code Distance",ylim=c(0,MAX))
points(data[,1],data[,82],col="blue")

colors = c("blue","red")
names = c("HGT","No HGT")
legend("topright",names,text.col=colors,bty="n")

for(i in 2:81){
    if(i == 2){
        smallest = min(data[,i])
        biggest = max(data[,i])
        si = i
        } else {
            if(smallest > min(data[,i])){
                smallest = min(data[,i])
                si = i
            }
            if(biggest < max(data[,i]))
                biggest = max(data[,i])
        }
}
smallest

for(i in 2:81){
    if(i == 2){
        smallestn = min(datan[,i])
        biggestn = max(datan[,i])
        sin = i
        } else {
        if(smallestn > min(datan[,i])){
            smallestn = min(data[,i])
            sin = i
        }
        if(biggestn < max(datan[,i]))
            biggestn = max(datan[,i])
        }
}
smallestn

if(biggest < biggestn)
    biggest = biggestn

if(smallest <= smallestn)
    plot(data[,1],data[,si],type = "l",col="blue",main="Mean Amino Acid Distance Over Time",xlab="Time Steps",ylab="Mean Amino Acid Distance Between Codon Neighbors",ylim=c(smallest,biggest)) else
    plot(datan[,1],datan[,sin],type = "l",col="red",main="Mean Amino Acid Distance Over Time",xlab="Time Steps",ylab="Mean Amino Acid Distance Between Codon Neighbors",ylim=c(smallestn,biggest))

for(i in 2:(si-1)){
    lines(data[,1],data[,i],col="blue")
}
for(i in (si+1):81){
    lines(data[,1],data[,i],col="blue")
}
for(i in 2:(sin-1)){
    lines(datan[,1],datan[,i],col="red")
}
for(i in (sin+1):81){
    lines(datan[,1],datan[,i],col="red")
    }
legend("topright",names,text.col=colors,bty="n")
