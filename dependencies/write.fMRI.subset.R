# Write text files with regional values, for plotting of cortical surface plots using code by Kirstie J Whitaker; 
# see https://github.com/KirstieJane/DESCRIBING_DATA/wiki/Making-Surface-Plots-with-Pysurfer

# Frantisek Vasa (fdv247@gmail.com)

write.fMRI.subset=function(v,ids,nroi,path) {
  # v = vector of subset values
  # ids of subset values in 1:nroi
  # assumes 16 subcortical regions
  # path to write file out to
  v.subset = array(NA,dim=c(nroi,1))
  v.subset[ids] = v
  write(v.subset, file = path, ncolumns = 1)
}