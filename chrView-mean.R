#chrView-mean plot

chrvwlist=read.table("chrView_mean_list.txt",sep='\t', stringsAsFactors=FALSE)

n_rain=nrow(chrvwlist)
cl <- rainbow(n_rain)
viewfile1=read.table(chrvwlist[1,2],sep='\t',header=T)
ym_cg=max(viewfile1$meanCG)*1.5*100
ym_chg=max(viewfile1$meanCHG)*1.5*100
ym_chh=max(viewfile1$meanCHH)*1.5*100

ys_cg=min(viewfile1$meanCG)*1.5*100
ys_chg=min(viewfile1$meanCHG)*1.5*100
ys_chh=min(viewfile1$meanCHH)*1.5*100


#CG
outfile=paste("./chrView_mean_CG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CG methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_cg,ym_cg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5)
#axis(2,at=seq(ys_cg,ym_cg))
for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCG*100,lwd=0.5,col= cl[i])
		lines(sub$Position,vector("numeric", length(sub$Position)),type='h',lwd=0.5,col='black') #draw zero line
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()


#CHG
outfile=paste("./chrView_mean_CHG.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CHG methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_chg,ym_chg),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCHG*100,lwd=0.5,col= cl[i])
		lines(sub$Position,vector("numeric", length(sub$Position)),type='h',lwd=0.5,col='black') #draw zero line
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()


#CHH
outfile=paste("./chrView_mean_CHH.png",sep="")
bitmap(outfile,type="png16m",height=2,width=8,res=400)

plot(NA,ylab=expression(paste(Delta,' CHH methylation levels (%)')),xlab="Chromosome",col='white',pch=19,type='p',cex.main = 1,frame = FALSE,ylim=c(ys_chh,ym_chh),xlim=c(0,nrow(viewfile1)),xaxt="n",cex=1.5)

for (i in 1:nrow(chrvwlist)){
	name=chrvwlist[i,1]
	file=chrvwlist[i,2]
	viewfile=read.table(file,sep='\t',header=T)
	
	chr=unique(viewfile$chromosome)
	for (j in chr){
		sub=viewfile[viewfile$chromosome==j,]
		lines(sub$Position,sub$meanCHH*100,lwd=0.5,col= cl[i])
		lines(sub$Position,vector("numeric", length(sub$Position)),type='h',lwd=0.5,col='black') #draw zero line
		#last row
		abline(v=sub[nrow(sub),1]+0.5,col="grey",lwd=0.5)
		#label
		labsite=sub[as.integer(nrow(sub)/2),1]
		lab=sub$chromosome[1]
		sta_p=sub[1,1]
		end_p=sub[(nrow(sub)),1]

axis(1,at=labsite ,label=lab,tick=F)
axis(1,at=sta_p ,label=FALSE)
axis(1,at=end_p ,label=FALSE)

xx<-c(sta_p,end_p)
axis(1,at=xx,labels=F)
}
}
legend("topright", legend = t(chrvwlist[1]), col = t(cl), lwd = 1,cex = 0.5)
dev.off()


