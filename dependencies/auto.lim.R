# Automatic axis limits when plotting single input variable x,
# Frantisek Vasa (fdv247@gmail.com)

auto.lim=function(x) {
  lim.d = max(x)-min(x)   # difference between max and min
  lim.incr = lim.d/20;    # increment from max and min (hardcoded as 20th of lim.d)
  lim.l = min(x)-lim.incr # lower limit
  lim.u = max(x)+lim.incr # upper limit
  return(c(lim.l,lim.u))
}
