#######################################
## PLOT MODEL CURVES AND DERIVATIVES ##
#######################################

library(tikzDevice)
library(ggplot2)
library(epigrowthfitPNAS)

## plot.setup(pdfname,width=7,height=7)

## texname = paste0(rtargetname, ".tex")
texname <- "plotmodels.tex"  ## avoid make magic
## plot.setup(texname,width=7,height=7)
tikz(texname,width=7,height=7, standAlone=TRUE) ## avoid make magic

##mylwd <- 4 ## pdf
mylwd <- 6 ## tikz
col.richards <- "black"
col.logistic <- "red"
col.exp <- "grey"
xmax <- 10
seg.len <- 4.5

c0 <- 1e-2

## high and low values
s.hi <- 1.3
s.lo <- 0.7
K.hi <- 1.3
K.lo <- 0.7

x0fun <- function(c0,K,s) {(c0/K)^s}

plot.configure <- function(ymax=1,xlab="Time",ylab="",...) {
    plot(NA,NA,xlim=c(0,xmax),ylim=c(0,ymax),xaxs="i",yaxs="i",bty="L",las=1,
         xlab=xlab,ylab=ylab,...)
}

plot.richards <- function(lty="solid",...) {
    mr <- get_model("richards",s_lower=0)
    if (names(dev.cur()) != "tikz output")
        plot(mr,add=TRUE,col="grey",lwd=1,lty="solid",...)
    plot(mr,add=TRUE,col=col.richards,lwd=mylwd,lty=lty,...)
}
plot.logistic <- function(...) {
    plot(get_model("logistic"),add=TRUE,col=col.logistic,lwd=mylwd,...)
}
plot.exp <- function(...) {
    plot(get_model("exp"),add=TRUE,col=col.exp,lwd=mylwd,...)
}

par(mfrow=c(2,2))
par(cex=1.15)
## margins:  [default is c(b,l,t,r) = c(5, 4, 4, 2) + 0.1]
mar.orig <- par("mar")
names(mar.orig) <- c("bottom", "left", "top", "right")
mymar <- mar.orig
mymar["right"] <- mymar["right"]*0.3
mymar["top"] <- mymar["top"]*0.15
mymar["bottom"] <- mymar["bottom"]*0.4
par(mar=mymar)

#################################
## CUMULATIVE INCIDENCE MODELS ##
#################################

plot.configure(ymax=1.05,ylab="Cumulative Incidence")
plot.exp()
K <- 1
s <- s.hi; plot.richards(par=c(s=s,x0=x0fun(c0,K,s)),lty="dashed")
s <- 1;    plot.richards(par=c(s=s,x0=x0fun(c0,K,s)),lty="solid")
s <- s.lo; plot.richards(par=c(s=s,x0=x0fun(c0,K,s)),lty="dotted")
legend("topleft",bty="n",legend="$K=1$")
if (names(dev.cur()) != "tikz output")
  legend("bottomright",bty="n",lwd=1,lty="solid",col="grey",
         legend=c(s.hi,1.0,s.lo),title="",seg.len=seg.len)
legend("bottomright",bty="n",lwd=mylwd,lty=c("dashed","solid","dotted"),
       legend=c(s.hi,1.0,s.lo),title="$s$",seg.len=seg.len)

plot.configure(ymax=1.05,ylab="Cumulative Incidence")
plot.exp()
s <- 1
K <- K.hi; plot.richards(par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="dashed")
K <- 1;    plot.richards(par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="solid")
K <- K.lo; plot.richards(par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="dotted")
legend("topleft",bty="n",legend="$s=1$")
if (names(dev.cur()) != "tikz output")
  legend("bottomright",bty="n",lwd=1,lty="solid",col="grey",
         legend=c(K.hi,1.0,K.lo),title="",seg.len=seg.len)
legend("bottomright",bty="n",lwd=mylwd,lty=c("dashed","solid","dotted"),
       legend=c(K.hi,1.0,K.lo),title="$K$",seg.len=seg.len)

######################
## INCIDENCE MODELS ##
######################

plot.configure(ymax=0.33,ylab="Incidence")
plot.exp(cumulative=FALSE)
K <- 1
s <- s.hi; plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s)),lty="dashed")
s <- 1;    plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s)),lty="solid")
s <- s.lo; plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s)),lty="dotted")

plot.configure(ymax=0.33,ylab="Incidence")
plot.exp(cumulative=FALSE)
s <- 1
K <- K.hi; plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="dashed")
K <- 1;    plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="solid")
K <- K.lo; plot.richards(cumulative=FALSE,par=c(s=s,x0=x0fun(c0,K,s),K=K),lty="dotted")

##plot.close(pdfname)
dev.off() # avoid make magic
