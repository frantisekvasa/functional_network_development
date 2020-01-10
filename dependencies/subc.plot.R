# Function to plot subcortical regions (assumes bilateral regions, ordered as LH -> RH)
#
# Inputs:
# vec       vector of subcortical values to plot
# subc.ord  order of subcortical regions along y-axis
# nm.subc   names of subcortical regions
# colBar    colorbar
# col.l     lower limit of colorbar
# col.u     upper limit of colorbar                                                                   
#
# Additional inputs for "subc.plot.sig"
# p.vec     vector of P-values
# thr       significance threshold (only markers for regions "i" with "vec[i][p.vec[i]] < thr" are plotted)
#
# Frantisek Vasa, fdv247@gmail.com, August 2018

subc.plot = function(vec,subc.ord,nm.subc,colBar,col.l,col.u) {

  nsubc = length(vec)/2 # number of (bilateral) subcortical regions
  
  col.sc = matrix(colBar[ceiling(100*((vec-col.l)/(col.u-col.l)))])                           # create color-scale
  
  plot(0, type='n', xaxt='n', yaxt='n', xlim=c(col.l,col.u), ylim=c(0.5,nsubc+0.5),xlab='',ylab='') # initialise plot
  
  # tick marks on axes
  if ((col.l+col.u)==0) {
    axis(1, at=seq(col.l,col.u,by=col.u)); axis(1, at=c(col.l/2,col.u/2),labels=c('',''))
  } else {
    axis(1, at=c(col.l,(col.l+col.u)/2,col.u)); axis(1, at=c( col.l+(col.u-col.l)/4 , col.u-(col.u-col.l)/4 ) ,labels=c('',''))
  }
  axis(2, at=1:8, labels=nm.subc[rev(subc.ord)], las=1)
  
  abline(v=0,lty=2,col='grey') # vertical line at 0
  
  # LH / RH region ID's
  lh.id = 1:nsubc               # LH regions
  rh.id = (nsubc+1):length(vec) # RH regions
  
  # 1:8 = left hemsiphere (symbols are arrows pointing left)
  points(vec[lh.id][rev(subc.ord)],1:nsubc,col='grey',pch=60,cex=2);                        # first plot markers in grey, slightly larger (=="outline")
  points(vec[lh.id][rev(subc.ord)],1:nsubc,col=col.sc[lh.id][rev(subc.ord)],pch=60,cex=1.75); # then plot markers in color
  # 9:16 = right hemisphere (symbols are arrows pointing right)
  points(vec[rh.id][rev(subc.ord)],1:nsubc,col='grey',pch=62,cex=2);                         # first plot markers in grey, slightly larger (=="outline")
  points(vec[rh.id][rev(subc.ord)],1:nsubc,col=col.sc[rh.id][rev(subc.ord)],pch=62,cex=1.75); # then plot markers in color
   
}

subc.plot.sig = function(vec,p.vec,thr,subc.ord,nm.subc,colBar,col.l,col.u) {
  
  nsubc = length(vec)/2 # number of (bilateral) subcortical regions
  
  col.sc = matrix(colBar[ceiling(100*((vec-col.l)/(col.u-col.l)))])                           # create color-scale
  
  plot(0, type='n', xaxt='n', yaxt='n', xlim=c(col.l,col.u), ylim=c(0.5,nsubc+0.5),xlab='',ylab='') # initialise plot
  
  # tick marks on axes
  if ((col.l+col.u)==0) {
    axis(1, at=seq(col.l,col.u,by=col.u)); axis(1, at=c(col.l/2,col.u/2),labels=c('',''))
  } else {
    axis(1, at=c(col.l,(col.l+col.u)/2,col.u)); axis(1, at=c( col.l+(col.u-col.l)/4 , col.u-(col.u-col.l)/4 ) ,labels=c('',''))
  }
  axis(2, at=1:8, labels=nm.subc[rev(subc.ord)], las=1)
  
  abline(v=0,lty=2,col='grey') # vertical line at 0
  
  # LH / RH region ID's
  lh.id = 1:nsubc               # LH regions
  rh.id = (nsubc+1):length(vec) # RH regions
  
  for (i in 1:length(vec)) {
    if (p.vec[i]<=thr) {
      # 1:8 = left hemsiphere (symbols are arrows pointing left)
      points(vec[lh.id][rev(subc.ord)][i],c(1:nsubc)[i],col='grey',pch=60,cex=2);                        # first plot markers in grey, slightly larger (=="outline")
      points(vec[lh.id][rev(subc.ord)][i],c(1:nsubc)[i],col=col.sc[lh.id][rev(subc.ord)][i],pch=60,cex=1.75); # then plot markers in color
      # 9:16 = right hemisphere (symbols are arrows pointing right)
      points(vec[rh.id][rev(subc.ord)][i],c(1:nsubc)[i],col='grey',pch=62,cex=2);                         # first plot markers in grey, slightly larger (=="outline")
      points(vec[rh.id][rev(subc.ord)][i],c(1:nsubc)[i],col=col.sc[rh.id][rev(subc.ord)][i],pch=62,cex=1.75); # then plot markers in color
    }
  }
  
}