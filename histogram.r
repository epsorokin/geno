# Visualize sampling distributions of Fst and PGx Fst values 
require(ggplot2)
rm(list=ls())
curdir='~/path/to/todays/experiment'

myfiles=list.files(curdir,pattern='Try1.tsv')

dat=data.frame(Group=as.character(),Fst=as.character())
for (f in myfiles) {
   tmp=read.table(paste('results',f,sep='/'),header=T,sep='\t',stringsAsFactors=F)
   dat=rbind(dat,tmp)
   rm(tmp)
}

pgxFile='out.RatioOfAvgs.ADME-core.tsv' #FIXME: is in another directory now
pgx=read.table(pgxFile,header=T)


df=data.frame(Group=as.character(),Mean_backgd.Fst=as.numeric(),Mean_backgd.sd=as.numeric(), Pgx.Fst=as.numeric(),Zscore=as.numeric(),pvalue=as.numeric())

for (g in pgx$Group) {
  m=mean(dat[dat$Group==g,]$Fst)
  s=sd(dat[dat$Group==g,]$Fst)
  tmp=data.frame(Group=g, Mean_backgd.Fst=m, Mean_backgd.sd=s, Pgx.Fst=pgx[pgx$Group==g,]$Fst,
    Zscore=(pgx[pgx$Group==g,]$Fst-m)/s, pvalue=1-pnorm ((pgx[pgx$Group==g,]$Fst-m)/s ) )
  df=rbind(df,tmp)
  rm(tmp)
}

write.table(df,'out.ADME-core.Zscore.pval.tsv',sep='\t',quote=F,row.names=F)

# Visualizations 
pdf('temp.PGx.v.SamplingDist.082117.pdf',width=10,height=7)

p1=ggplot(dat,aes(x=Fst))+geom_histogram(fill='lavender')+facet_wrap(~Group,scales='free')+geom_vline(data=pgx,aes(xintercept=Fst),color='red',linetype='dashed' ) +theme_classic()
print(p1)
dev.off()

pdf('temp.SamplingDensity.per-Pop.pdf',width=6,height=6)
p2=ggplot(dat,aes(x=Fst))+geom_density(aes(color=Group))+theme_classic()+scale_color_brewer(palette='Spectral')
print(p2)
dev.off()
