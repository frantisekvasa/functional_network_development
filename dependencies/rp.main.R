# Formats r-squared and p-values within titles (mainly due to bugs in interaction of rstudio and the b-quote function)

# written by Frantisek Vasa (fdv247@gmail.com)

# Pearson's r + P
rp.main = function(rsq,p,nsig) {
  main = bquote(r^2 ~ '=' ~ .(toString(signif(rsq,nsig))) * ', p =' ~ .(toString(signif(p,nsig)))) 
  return(main)
}

# Spearman's rho + P
rp.main.sp = function(sp.rho,p,nsig) {
  if (p < 1e-10) {
    main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p < 1e-10') 
  } else {
    main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p =' ~ .(toString(signif(p,nsig))))
  }
  return(main)
}

# Spearman's rho + P + Pspin
rp.main.sp.spin = function(sp.rho,p,p.spin,nsig) {
  if (p < 1e-10) {
    #main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p < 1e-10')
    if (p.spin < 1e-4) {
      main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p < '*10^{-10} * paste(', p'['spin'],' < '*10^{-4})) 
    } else { # if p.spin > 1e-4
      main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p < '*10^{-10} * paste(', p'['spin'],' = '*.(toString(signif(p.spin,nsig))) ))
    }
  } else { # if p > 1e-10
    if (p.spin < 1e-4) {
      main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p =' ~ .(toString(signif(p,nsig))) * paste(', p'['spin'],' < '*10^{-4})) 
    } else { # if p.spin > 1e-4
      main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p =' ~ .(toString(signif(p,nsig))) * paste(', p'['spin'],' = '*.(toString(signif(p.spin,nsig))) ))
    }
    return(main)
  }
}

