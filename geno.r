#!/usr/bin/env Rscript

#total up the lofs per population, total up lofs per individual 
#usage ./mini_geno.r pop.indivs.raw

#get arguments
require(data.table)
arguments <- commandArgs(trailingOnly=TRUE)
for (i in 1:length(arguments)) {
  print(paste("arg",as.character(i),"=",arguments[i]))
  indivstowrite=paste('out.indivs.',strsplit(arguments[i],'.indivs.raw')[[1]],'.tsv',sep='')
  lofstowrite=paste('out.lofs.', strsplit(arguments[i],'.indivs.raw')[[1]],'.tsv',sep='') 
  dat=fread(arguments[i],header=T,check.names=F)
  print (arguments[i])
  row.names(dat)=dat$IID

  # Number homozygous LoFs per individual
  count.2.across=function(index) length(grep(2,dat[index,7:ncol(dat),with=F]))
  num_lof_variants=as.integer(sapply(1:nrow(dat),count.2.across))
  lofs=data.frame(num_lof_variants, row.names=row.names(dat))
  write.table(lofs,lofstowrite,sep='\t',quote=F,row.names=T)
  
  # Number individuals per homozygous LoF 
  count.2.down=function(index) length(grep(2,as.matrix(dat[,index,with=F])))  
  variant=colnames(dat)[7:ncol(dat)]
  total_indivs_w_lof=sapply(7:ncol(dat),count.2.down)
  indivs=data.frame(variant, total_indivs_w_lof)
  write.table(indivs,indivstowrite,sep='\t',quote=F,row.names=F)
  print (indivstowrite)
  print (lofstowrite)

}

# example of duplicate row to filter out: 1:1246013-G-A_A