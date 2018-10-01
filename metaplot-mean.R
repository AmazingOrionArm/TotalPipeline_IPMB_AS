args = commandArgs(trailingOnly=TRUE)
samplelist=read.table("samples_list.txt",sep='\t')

#n_rain=nrow(samplelist)
#cl <- rainbow(n_rain)
metafile1=read.table(paste(samplelist[1,1],"_CG.matrix.gz",sep=''),skip=1)
metafile2=read.table(paste(samplelist[1,1],"_CHG.matrix.gz",sep=''),skip=1)
metafile3=read.table(paste(samplelist[1,1],"_CHH.matrix.gz",sep=''),skip=1)
y_cg=max(apply(as.matrix(metafile1[,7:86]),2,mean,na.rm=T))*1.2*100
y_chg=max(apply(as.matrix(metafile2[,7:86]),2,mean,na.rm=T))*1.2*100
y_chh=max(apply(as.matrix(metafile3[,7:86]),2,mean,na.rm=T))*1.2*100

outfile=paste("./metaplot_mean_CG.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CG methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(0,y_cg),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CG.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si] && args[2]!=args[3]) {
  av=av2-av1
}
hy_cg=ifelse(max(av)>=0, max(av)*1.2*100, max(av)*0.8*100)
ly_cg=ifelse(min(av)>=0, min(av)*0.8*100, min(av)*1.2*100)
plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_cg,hy_cg),xaxt="n",cex=1.5)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(vector("numeric", length(av)),xlab=NA,type='h',cex.main = 1,xaxt="n",cex=1.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

#legend("topright", legend = t(samplelist[1]), col = t(palette()), lwd = 1,cex = 0.5)
dev.off()

#---

outfile=paste("./metaplot_mean_CHG.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CHG methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile2)),ylim=c(0,y_chg),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CHG.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si] && args[2]!=args[3]) {
  av=av2-av1
}
hy_chg=ifelse(max(av)>=0, max(av)*1.2*100, max(av)*0.8*100)
ly_chg=ifelse(min(av)>=0, min(av)*0.8*100, min(av)*1.2*100)
plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_chg,hy_chg),xaxt="n",cex=1.5)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(vector("numeric", length(av)),xlab=NA,type='h',cex.main = 1,xaxt="n",cex=1.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

#legend("topright", legend = t(samplelist[1]), col = t(palette()), lwd = 1,cex = 0.5)
dev.off()

#---

#CHH
outfile=paste("./metaplot_mean_CHH.png",sep="")
bitmap(outfile,type="png16m",height=3.5,width=3.5,res=450)
#plot(NA,ylab='CHH methylation levels(%)',xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile3)),ylim=c(0,y_chh),xaxt="n",cex=1.5)

#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
#palette()[1] "black"   "red"     "green3"  "blue"    "cyan"    "magenta" "yellow" "gray"
si=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,3]
  if(i>1)
    if(name!=samplelist[i-1,3])
      si=i
}
av=0
av1=0
av2=0
for (i in 1:nrow(samplelist)){
  name=samplelist[i,1]
  file=paste(name,"_CHH.matrix.gz",sep='')
  metafile=read.table(file,sep='\t',skip=1)
  if(i<si){
    av1=av1+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }else{
    av2=av2+apply(as.matrix(metafile[,7:86]),2,mean,na.rm=T)
  }
}
av1=av1/(si-1)
av2=av2/((length(samplelist)+1)-(si-1))
if(args[2]==samplelist[si-1,3] && args[2]!=args[3]){
  av=av1-av2
}else if (args[2]==samplelist[si] && args[2]!=args[3]) {
  av=av2-av1
}
hy_chh=ifelse(max(av)>=0, max(av)*1.2*100, max(av)*0.8*100)
ly_chh=ifelse(min(av)>=0, min(av)*0.8*100, min(av)*1.2*100)
plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab=NA,type='l',cex.main = 1,xlim=c(0,ncol(metafile1)),ylim=c(ly_chh,hy_chh),xaxt="n",cex=1.5)
lines(av*100,xlab=NA,type='l',cex.main = 1,xaxt="n",cex=1.5,col=palette()[i])
lines(vector("numeric", length(av)),xlab=NA,type='h',cex.main = 1,xaxt="n",cex=1.5,col='black') #draw zero line

xx<-c(20,40,60)
ind<-c("TSS",args[1],"TES")

abline(v=20,lty=2,col='grey70')
abline(v=60,lty=2,col='grey70')
axis(1,at=xx,labels=ind)

#legend("topright", legend = t(samplelist[1]), col = t(palette()), lwd = 1,cex = 0.5)
dev.off()


