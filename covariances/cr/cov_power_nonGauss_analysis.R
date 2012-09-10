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
    for (j in c(1:Npower)){
      for (k in c(1:Ncl)){
	for (l in c(k:Ncl)){
	  mat[(i-1)*(Ncl)+k,(j-1)*(Ncl)+l]=filecolumn[m]
	  mat[(j-1)*(Ncl)+l,(i-1)*(Ncl)+k]=filecolumn[m]
	  m=m+1
	}
      }
    }
  }
  
  return(mat)
}


ALL <- read.table(paste("../../../predictions/test20cov_power_no_tomo_Start0_End1"))

Cl1 <- ReadMatrix(c(ALL$V1))
Cl2 <- ReadMatrix(c(ALL$V2))
noise <-ReadMatrix(c(ALL$V3))
clin  <-ReadMatrix(c(ALL$V4))          
cquad <-ReadMatrix(c(ALL$V5))          
gauss_tot <-noise+clin+cquad

trilin  <- ReadMatrix(c(ALL$V6))          
cov1h  <- ReadMatrix(c(ALL$V7))          
cov2h  <- ReadMatrix(c(ALL$V8))          
hsvterm  <- ReadMatrix(c(ALL$V9))          
covtot <- noise+clin+cquad+trilin+cov1h+cov2h+hsvterm


test <- cbind(matrix(Cl1,Npower*Ncl*Npower*Ncl,1),matrix(Cl2,Npower*Ncl*Npower*Ncl,1),matrix(noise,Npower*Ncl*Npower*Ncl,1),matrix(clin,Npower*Ncl*Npower*Ncl,1),matrix(cquad,Npower*Ncl*Npower*Ncl,1),matrix(trilin,Npower*Ncl*Npower*Ncl,1),matrix(cov1h,Npower*Ncl*Npower*Ncl,1),matrix(cov2h,Npower*Ncl*Npower*Ncl,1),matrix(hsvterm,Npower*Ncl*Npower*Ncl,1),matrix(covtot,Npower*Ncl*Npower*Ncl,1))

write(t(test),file=paste("../covtot_easyread"),ncolumns=10)
NEW <- read.table(paste("../covtot_easyread"))


cg1 <- diag(noise+clin+cquad)
cg2 <- diag(noise+clin+cquad)

normcovgauss <-sqrt(cg1%*%t(cg2))

trilinplot <-trilin/normcovgauss
cov1hplot  <-cov1h/normcovgauss
cov2hplot  <-cov2h/normcovgauss
hsvplot <-hsv/normcovgauss


inv <- solve(covtot)
eigeninv <- eigen(inv, only.values=TRUE)
eigencov <- eigen(covtot, only.values=TRUE)

ell <- c(ALL$V2[1:Ncl])

covtot <- cov2hplot 
zr<- range(min(log(covtot)),max(log(covtot)))

postscript(paste("../plots/cov2hplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
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

postscript(paste("../plots/trilinplot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
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



# postscript(paste("../plots/covhsv2.eps",sep=""),height=8,width=8,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
# # par(omi=c(0.9,0.9,0.2,0.7),mai=c(0.5,0.5,0.0,0.0),font=2,tck=-0.02,cex.axis=1.2)
# covtotflip <-as.matrix(hsvplot[,300:1])
# ver <- c("z11", "z12", "z13", "z14", "z15","z22", "z23", "z24", "z25", "z33","z34", "z35", "z44", "z45", "z55")
# hor <- paste(rev(ver))
# #dimnames(covtotflip)<-list(hor,ver)
# 
# rgb.palette <- colorRampPalette(c("blue","green","yellow", "red"), space = "rgb")
# #leg <-expression(paste(gamma , gamma),paste(mu , mu),paste(p , p),paste(gamma , mu),paste(gamma ,p),paste(mu ,p))
# levelplot(covtotflip,xlab="",ylab="",panel = function(...) {
  #         panel.levelplot(...)
#         panel.abline(h=c(0,20,40,60,80,5,100,120,140,160,180,200,220,240,260,280,300))
#         panel.abline(v=c(0,20,40,60,80,5,100,120,140,160,180,200,220,240,260,280,300))
# 	panel.axis(side=c("top"), at=c(10, 30, 50, 70,90,110,130,150,170,190,210,230,250,270,290),labels=leg,rot=0, half=FALSE, ticks=FALSE)
# 	panel.axis(side=c("left"), at=c(10, 30, 50, 70,90,110,130,150,170,190,210,230,250,270,290),labels=rev(leg),rot=0, half=FALSE, ticks=FALSE)
# },col.regions=rgb.palette(500),cuts=500)
# 
# dev.off()
# 
# 
# 




postscript(paste("../plots/covhsv2.eps",sep=""),height=8,width=8,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")
# par(omi=c(0.9,0.9,0.2,0.7),mai=c(0.5,0.5,0.0,0.0),font=2,tck=-0.02,cex.axis=1.2)
covtotflip <-as.matrix(hsvplot[,300:1])
ver <- c("z11", "z12", "z13", "z14", "z15","z22", "z23", "z24", "z25", "z33","z34", "z35", "z44", "z45", "z55")
hor <- paste(rev(ver))
#dimnames(covtotflip)<-list(hor,ver)

rgb.palette <- colorRampPalette(c("blue","green","yellow", "red"), space = "rgb")
#leg <-expression(paste(gamma , gamma),paste(mu , mu),paste(p , p),paste(gamma , mu),paste(gamma ,p),paste(mu ,p))
levelplot(covtotflip,xlab="",ylab="",xaxt="n",yaxt="n",col.regions=rgb.palette(500),cuts=500)

dev.off()







postscript(paste("../plots/covtot.eps",sep=""),height=9.2,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")

par(mai=c(.4,.4,0.2,0.1),mgp=c(0.92,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)
image.plot(ell,ell,log(covtot),col=tim.colors(),xlab="",ylab="",xaxt="n",yaxt="n")
dev.off()   




postscript(paste("../plots/cov1h.eps",sep=""),height=6.6,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")

par(mai=c(.4,.4,0.2,0.1),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)
image.plot(ell,ell,log(cov1hplot),col=tim.colors(),xlab="",ylab="",xaxt="n",yaxt="n")

dev.off()   

postscript(paste("../plots/cov2h.eps",sep=""),height=6.6,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")

par(mai=c(.4,.4,0.2,0.1),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)
image.plot(ell,ell,log(cov2hplot),col=tim.colors(),xlab="",ylab="",xaxt="n",yaxt="n")

dev.off()   

postscript(paste("../plots/covtrilin.eps",sep=""),height=6.6,width=10,encoding = "TeXtext.enc",onefile=FALSE, horizontal=FALSE, paper="special")

par(mai=c(.4,.4,0.2,0.1),mgp=c(0.9,.1,0),tck=-0.01,cex=1.0,cex.main=1.0)
image.plot(ell,ell,log(trilinplot),col=tim.colors(),xlab="",ylab="",xaxt="n",yaxt="n")

dev.off()   


