#!/usr/bin/env Rscript

args = commandArgs(TRUE)


sharead = function(chidata,oardata,sh_chifile,sh_oarfile) {
  chi_btb=read.table(chidata, sep = "\t")
  oar_btb=read.table(oardata,sep = "\t")
  oar_btb$id= 1: nrow(oar_btb)
  shared2=merge(chi_btb,oar_btb, by="V4")
  pshared=shared2[order(shared2$id),]
  sh_chi=data.frame("v1"= pshared$V1.x, "v2"=pshared$V2.x, "v3"=pshared$V3.x,     "v4"=pshared$V4, "v5"=pshared$V5.x, "v6"=pshared$V6.x)
  sh_oar=data.frame("v1"= pshared$V1.y, "v2"=pshared$V2.y, "v3"=pshared$V3.y, "v4"=pshared$V4, "v5"=pshared$V5.y, "v6"=pshared$V6.y)
  write.table(sh_chi,file = sh_chifile, sep = "\t", col.names = F,row.names = F, quote = F)
  write.table(sh_oar,file = sh_oarfile, sep = "\t", col.names = F,row.names = F, quote = F)
  
}


sharead(args[1],args[2],args[3],args[4])
