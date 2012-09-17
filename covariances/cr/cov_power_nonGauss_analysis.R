library(graphics)
library(lattice)
library(plotrix)
library(fields)

arcmin=2.90888208665721580e-4
Npower <- 15 
Ncl <-20

ReadMatrix <- function(filecolumn,npar=TRUE,print=TRUE)
{
  mat<-array(0, dim=c(Npower*Ncl,Npower*Ncl))
  
  m<-1
  for (i in c(1:Npower)){
    for (j in c(i:Npower)){
      for (k in c(1:Ncl)){
	for (l in c(k:Ncl)){
	  mat[(i-1)*(Ncl)+k,(j-1)*(Ncl)+l]=filecolumn[m]
	  mat[(i-1)*(Ncl)+l,(j-1)*(Ncl)+k]=filecolumn[m]
	
	  mat[(j-1)*(Ncl)+l,(i-1)*(Ncl)+k]=filecolumn[m]
	  mat[(j-1)*(Ncl)+k,(i-1)*(Ncl)+l]=filecolumn[m]
	m=m+1
	}
      }
    }
  }
  
  return(mat)
}


ALL <- read.table(paste("../../../data/baryons/OWLS_cov_power_lmin=3.409e+01_lmax=4.400e+03_Ncl=20_tomo_NonGauss"))

gauss <-ReadMatrix(c(ALL$V7))
trilin  <- ReadMatrix(c(ALL$V8))          
cov1h  <- ReadMatrix(c(ALL$V9))          
cov2h  <- ReadMatrix(c(ALL$V10))          
hsv  <- ReadMatrix(c(ALL$V11))          
covtotNG <- gauss+trilin+cov1h+cov2h+hsv

normcov <-covtotNG

trilinplot <-trilin/normcov
cov1hplot  <-cov1h/normcov
cov2hplot  <-cov2h/normcov
hsvplot <-hsv/normcov


inv <- solve(covtot)
eigeninv <- eigen(inv, only.values=TRUE)
eigencov <- eigen(covtot, only.values=TRUE)

ell <- c(ALL$V6[1:Ncl])



covtot <- cov1hplot 
zr<- range(min(log(covtot)),max(log(covtot)))

postscript(paste("../../../data/baryons/plots/cov1hplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
#par(mfrow=c(Npower,Npower),omi=c(.8,.8,0.2,0.2),mai=c(.0,.0,0.0,0.0),mgp=c(1.5,.2,0),tck=-0.01,cex=1.0,cex.main=1.0)
par(omi=c(.8,.8,0.0,0.0),mai=c(.0,.0,0.0,0.0),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)

split.screen( rbind(c(0.05,.85,0,.95), c(0.85,1,0,.95)))
 split.screen(c(15,15), screen=1)-> ind
k <- 1
 for(i in c(1:Npower))
{
  for(j in c(1:Npower))
  { 
  screen( ind[k])
    flip <-as.matrix(covtot[(((i-1)*Ncl)+1):(i*Ncl),(j*Ncl):(((j-1)*Ncl)+1)])
    image(log(ell),log(ell),log(flip),col=tim.colors(), zlim=zr,xlab="",ylab="",xaxt="n",yaxt="n")
  k=k+1
  }
}
screen( 2)
  logd <-(max(log(covtot))-min(log(covtot)))/5 
  tickvals <- formatC(c(min(covtot),exp(min(log(covtot))+logd),exp(min(log(covtot))+2*logd),exp(min(log(covtot))+3*logd),exp(min(log(covtot))+4*logd),max(covtot)), digits = 1, format = "e")
  ticks <-seq(min(log(covtot)),max(log(covtot)),(max(log(covtot))-min(log(covtot)))/5)
#  tickvals <-c(round(min(covtot),digits=4),round(max(covtot),digits=4))
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.3,.43,0.05,0.95),axis.args=list( at=ticks, labels=tickvals), col=tim.colors())

   close.screen( all=TRUE)
  test <- c("z11","z12","z13","z14","z15","z22","z23","z24","z25","z33","z34","z35","z44","z45","z55")
  spot <- seq(0.0, 0.81, length.out= 15)
  spot2 <- seq(-0.06, 0.91, length.out= 15)
  mtext(test,at=spot,side=3,line=-2.0,cex=1,outer=TRUE)
  mtext(rev(test),at=spot2,side=2,line=1.7,cex=1,outer=TRUE)
dev.off()   



covtot <- cov2hplot 
zr<- range(min(log(covtot)),max(log(covtot)))

postscript(paste("../../../data/baryons/plots/cov2hplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
#par(mfrow=c(Npower,Npower),omi=c(.8,.8,0.2,0.2),mai=c(.0,.0,0.0,0.0),mgp=c(1.5,.2,0),tck=-0.01,cex=1.0,cex.main=1.0)
par(omi=c(.8,.8,0.0,0.0),mai=c(.0,.0,0.0,0.0),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)

split.screen( rbind(c(0.05,.85,0,.95), c(0.85,1,0,.95)))
 split.screen(c(15,15), screen=1)-> ind
k <- 1
 for(i in c(1:Npower))
{
  for(j in c(1:Npower))
  { 
  screen( ind[k])
    flip <-as.matrix(covtot[(((i-1)*Ncl)+1):(i*Ncl),(j*Ncl):(((j-1)*Ncl)+1)])
    image(log(ell),log(ell),log(flip),col=tim.colors(), zlim=zr,xlab="",ylab="",xaxt="n",yaxt="n")
  k=k+1
  }
}
screen( 2)
  logd <-(max(log(covtot))-min(log(covtot)))/5 
  tickvals <- formatC(c(min(covtot),exp(min(log(covtot))+logd),exp(min(log(covtot))+2*logd),exp(min(log(covtot))+3*logd),exp(min(log(covtot))+4*logd),max(covtot)), digits = 1, format = "e")
  ticks <-seq(min(log(covtot)),max(log(covtot)),(max(log(covtot))-min(log(covtot)))/5)
#  tickvals <-c(round(min(covtot),digits=4),round(max(covtot),digits=4))
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.3,.43,0.05,0.95),axis.args=list( at=ticks, labels=tickvals), col=tim.colors())

   close.screen( all=TRUE)
  test <- c("z11","z12","z13","z14","z15","z22","z23","z24","z25","z33","z34","z35","z44","z45","z55")
  spot <- seq(0.0, 0.81, length.out= 15)
  spot2 <- seq(-0.06, 0.91, length.out= 15)
  mtext(test,at=spot,side=3,line=-2.0,cex=1,outer=TRUE)
  mtext(rev(test),at=spot2,side=2,line=1.7,cex=1,outer=TRUE)
dev.off()   





covtot <- trilinplot
zr<- range(min(covtot),max(covtot))

postscript(paste("../../../data/baryons/plots/trilinplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
#par(mfrow=c(Npower,Npower),omi=c(.8,.8,0.2,0.2),mai=c(.0,.0,0.0,0.0),mgp=c(1.5,.2,0),tck=-0.01,cex=1.0,cex.main=1.0)
par(omi=c(.8,.8,0.0,0.0),mai=c(.0,.0,0.0,0.0),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)

split.screen( rbind(c(0.05,.85,0,.95), c(0.85,1,0,.95)))
 split.screen(c(15,15), screen=1)-> ind
k <- 1
 for(i in c(1:Npower))
{
  for(j in c(1:Npower))
  { 
  screen( ind[k])
    flip <-as.matrix(covtot[(((i-1)*Ncl)+1):(i*Ncl),(j*Ncl):(((j-1)*Ncl)+1)])
    image(log(ell),log(ell),flip,col=tim.colors(), zlim=zr,xlab="",ylab="",xaxt="n",yaxt="n")
  k=k+1
  }
}
screen( 2)
  del <-(max(covtot)-min(covtot))/5 
  tickvals <- formatC(c(min(covtot),min(covtot)+del,min(covtot)+2*del,min(covtot)+3*del,min(covtot)+4*del,max(covtot)), digits = 1, format = "e")
  ticks <-seq(min(covtot),max(covtot),(max(covtot)-min(covtot))/5)
#  tickvals <-c(round(min(covtot),digits=4),round(max(covtot),digits=4))
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.3,.43,0.05,0.95),axis.args=list( at=ticks, labels=tickvals), col=tim.colors())

   close.screen( all=TRUE)
  test <- c("z11","z12","z13","z14","z15","z22","z23","z24","z25","z33","z34","z35","z44","z45","z55")
  spot <- seq(0.0, 0.81, length.out= 15)
  spot2 <- seq(-0.06, 0.91, length.out= 15)
  mtext(test,at=spot,side=3,line=-2.0,cex=1,outer=TRUE)
  mtext(rev(test),at=spot2,side=2,line=1.7,cex=1,outer=TRUE)
dev.off()   



covtot <- hsvplot
zr<- range(min(covtot),max(covtot))

postscript(paste("../../../data/baryons/plots/hsvplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
#par(mfrow=c(Npower,Npower),omi=c(.8,.8,0.2,0.2),mai=c(.0,.0,0.0,0.0),mgp=c(1.5,.2,0),tck=-0.01,cex=1.0,cex.main=1.0)
par(omi=c(.8,.8,0.0,0.0),mai=c(.0,.0,0.0,0.0),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)

split.screen( rbind(c(0.05,.85,0,.95), c(0.85,1,0,.95)))
 split.screen(c(15,15), screen=1)-> ind
k <- 1
 for(i in c(1:Npower))
{
  for(j in c(1:Npower))
  { 
  screen( ind[k])
    flip <-as.matrix(covtot[(((i-1)*Ncl)+1):(i*Ncl),(j*Ncl):(((j-1)*Ncl)+1)])
    image(log(ell),log(ell),flip,col=tim.colors(), zlim=zr,xlab="",ylab="",xaxt="n",yaxt="n")
  k=k+1
  }
}
screen( 2)
  del <-(max(covtot)-min(covtot))/5 
  tickvals <- formatC(c(min(covtot),min(covtot)+del,min(covtot)+2*del,min(covtot)+3*del,min(covtot)+4*del,max(covtot)), digits = 1, format = "e")
  ticks <-seq(min(covtot),max(covtot),(max(covtot)-min(covtot))/5)
#  tickvals <-c(round(min(covtot),digits=4),round(max(covtot),digits=4))
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.3,.43,0.05,0.95),axis.args=list( at=ticks, labels=tickvals), col=tim.colors())

   close.screen( all=TRUE)
  test <- c("z11","z12","z13","z14","z15","z22","z23","z24","z25","z33","z34","z35","z44","z45","z55")
  spot <- seq(0.0, 0.81, length.out= 15)
  spot2 <- seq(-0.06, 0.91, length.out= 15)
  mtext(test,at=spot,side=3,line=-2.0,cex=1,outer=TRUE)
  mtext(rev(test),at=spot2,side=2,line=1.7,cex=1,outer=TRUE)
dev.off()   


covtot <- cov2cor(covtotNG)
zr<- range(min(covtot),max(covtot))

postscript(paste("../../../data/baryons/plots/corrplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
#par(mfrow=c(Npower,Npower),omi=c(.8,.8,0.2,0.2),mai=c(.0,.0,0.0,0.0),mgp=c(1.5,.2,0),tck=-0.01,cex=1.0,cex.main=1.0)
par(omi=c(.8,.8,0.0,0.0),mai=c(.0,.0,0.0,0.0),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)

split.screen( rbind(c(0.05,.85,0,.95), c(0.85,1,0,.95)))
 split.screen(c(15,15), screen=1)-> ind
k <- 1
 for(i in c(1:Npower))
{
  for(j in c(1:Npower))
  { 
  screen( ind[k])
    flip <-as.matrix(covtot[(((i-1)*Ncl)+1):(i*Ncl),(j*Ncl):(((j-1)*Ncl)+1)])
    image(log(ell),log(ell),flip,col=tim.colors(), zlim=zr,xlab="",ylab="",xaxt="n",yaxt="n")
  k=k+1
  }
}
screen( 2)
  del <-(max(covtot)-min(covtot))/5 
  tickvals <- formatC(c(min(covtot),min(covtot)+del,min(covtot)+2*del,min(covtot)+3*del,min(covtot)+4*del,max(covtot)), digits = 1, format = "e")
  ticks <-seq(min(covtot),max(covtot),(max(covtot)-min(covtot))/5)
#  tickvals <-c(round(min(covtot),digits=4),round(max(covtot),digits=4))
   image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.3,.43,0.05,0.95),axis.args=list( at=ticks, labels=tickvals), col=tim.colors())

   close.screen( all=TRUE)
  test <- c("z11","z12","z13","z14","z15","z22","z23","z24","z25","z33","z34","z35","z44","z45","z55")
  spot <- seq(0.0, 0.81, length.out= 15)
  spot2 <- seq(-0.06, 0.91, length.out= 15)
  mtext(test,at=spot,side=3,line=-2.0,cex=1,outer=TRUE)
  mtext(rev(test),at=spot2,side=2,line=1.7,cex=1,outer=TRUE)
dev.off()   





