########

# Visualize site frequency spectrum for each of eleven largest populations within a large set of populations

#setwd()
require(data.table)
require(ggplot2)
require(grid)

rm(list=ls())

# Super populations
superpop=c('CentralAmerican','CentralSouthAsia','Cuba','DominicanRepublic','Mexico','AA','Asian','NatAm','NatHw','PuertoRico','SouthAmerican')


# Read in PAGE data; FIXME: Adjust population labels 
page=fread("/path/to/freq/ALL.POPS.TXT",header=T,sep="\t")
page=data.frame(page)
names(page)[271:272]=c('MAF.total','NCHROBS.total')
names(page)[1:4]=c("CHR.AA","SNP.AA", "A1.AA","A2.AA")
toMatch= c("CHR.AA", "SNP.AA", "A1.AA","A2.AA", "total",superpop)
matches <- grep(paste(toMatch,collapse="|"), colnames(page),value=T)  # Keep AF data get rid of redundant cols
page <- subset(page, select = matches)

Pops=c('total',superpop)
# Update 0627: Get rid of total 
Pops=superpop

# For each population, grep out a new column 
dat=data.frame(x="NA",y=1:dim(page)[1])
for (p in Pops) {
  # Subset out columns of interest 
  print(p)
  toMatch=grep(p, colnames(page),value=T)
  #print(length(toMatch)) # Make sure each match only picking up 2 columns
  temp=as.data.frame(subset(page,select =toMatch))
  temp.name=paste("AlC",p,sep=".")
  temp$V3=round(temp[,1] * temp[,2] ,digits=1)
  colnames(temp)[3]=temp.name
  head(temp)
  dat=cbind(dat, temp)
}

dat=cbind(page$SNP.AA, dat)



# What SNPs have more than 2 counts in a single population
popCountCols = grep( 'AlC.', colnames(dat), value=T )
popCounts = subset(dat,select=c("page$SNP.AA",popCountCols))
t=melt(popCounts) ## Need to split up data frame and put MAF as new column. 
names(t)=c("SNP","Population","Counts")
t$Population=gsub('AlC.','',t$Population)
MAFcols=grep('MAF.',colnames(dat),value=T)
popMAF = subset (dat,select=c("page$SNP.AA",MAFcols))
t2=melt(popMAF)
names(t2)=c("SNP","Population","MAF")
dfm=cbind(t,t2)
dfm=dfm[,c(1,2,3,6)]
dfm$Population=gsub('AlC.','',dfm$Population)
dfm$Population=mapvalues(dfm$Population,from=c('total','PAGE.AA','PAGE.Asian','CentralAmerican','CentralSouthAsia','Cuba','DominicanRepublic','PAGE.NatAm','PAGE.NatHw','Mexico',
                      'PuertoRico','SouthAmerican'),to=c('Total','AA','Asian','CentralAm','SouthAsn','Cuba','D.R.','NatAm','NatHw','Mexico',
                                                            'P.R.','SouthAm' ) )

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# All CRVs
m=ddply(dfm, 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'

s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/29680 * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',limits=c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%'))+
  theme(axis.text.y=element_text(vjust=0.6))+ylab('CRVs in PAGE/MEGA (%)') +
  theme(legend.title=element_blank(),legend.position='top',axis.title.y=element_blank()) + coord_flip()+
  scale_x_discrete(limits=rev(levels(s$Population)))
ggsave(filename='SFS_allCRV_legend.png',plot = last_plot(), width=8,height=6,units='in',dpi=300)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Facet by pathogenicity
clinvar=fread("~/path/to/may2017.dl/clinvar_allele_trait_pairs.single.b37.tsv",sep = '\t',header=T)
clinvar$partialid <- paste(clinvar$chrom, clinvar$pos, sep = ":")
clinvar$bimid <- paste(clinvar$partialid, clinvar$ref, clinvar$alt, sep="-")
benign=unique(clinvar[clinvar$clinical_significance %in% c("Benign","Likely benign"),]$bimid)
pathogenic=unique(clinvar[clinvar$clinical_significance %in% c("Pathogenic"),]$bimid)
cv=unique(clinvar$bimid)


# All CLinVar
m=ddply(dfm[dfm$SNP %in% cv,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% cv,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)

a=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('ClinVar variants (%)') +
  theme(axis.title.y=element_blank()) +
  coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_Clinvar.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

# Benign only
m=ddply(dfm[dfm$SNP %in% benign,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/8423 * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
b=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('Benign variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_benignLB.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

# Pathogenic
m=ddply(dfm[dfm$SNP %in% pathogenic,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/7699 * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
c=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('Pathogenic variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_pathogenic.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Facet by fxnal assertion
cadd=fread("~/path/to/cadd/COMBINED/out.lof+crv.cadd.tsv",header=T,sep="\t")
cadd$partial=paste(cadd$`#Chrom`, cadd$Pos,sep=":")
cadd$bimid=paste(cadd$partial,cadd$Ref,cadd$Alt,sep="-")
missense=unique(cadd[cadd$Consequence=="NON_SYNONYMOUS",]$bimid)
synon=unique(cadd[cadd$Consequence=="SYNONYMOUS",]$bimid)
ptv=unique(cadd[cadd$Consequence=="STOP_GAINED",]$bimid)

# Missense
m=ddply(dfm[dfm$SNP %in% missense,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% missense,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population,ordered=T)
h=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('Missense variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
ggsave(filename='SFS_missense.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

# Synonymous
m=ddply(dfm[dfm$SNP %in% synon,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% synon,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
g=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('Synonymous variants (%)') +
  theme(axis.title.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
ggsave(filename='SFS_synon.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

# PTV
m=ddply(dfm[dfm$SNP %in% ptv,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% synon,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
i=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('PTV variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
ggsave(filename='SFS_synon.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# # Facet by LoF
# lofs=fread('~/path/to/LoF-analyses.Apr01/LOF.ANNOTATIONS/awked.Freeze1.HC.lofs.txt',sep=' ',header=F)
# lof=unique(lofs$V1)
# m=ddply(dfm[dfm$SNP %in% lof,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
# m$c=as.character(m$c)
# m[m$Counts==1,]$c='singleton'
# m[m$MAF==0,]$c='absent'
# s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% lof,]$SNP)) * 100, digits=1) )
# s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
#            labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
# s$Population=factor(s$Population)
# i=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
#   scale_fill_brewer(palette='BrBG',direction=-1)+
#   theme(axis.text.x=element_text(vjust=0.6))+ylab('predicted LoF variants (%)') +
#   theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
#   guides(fill=F)+
#   scale_x_discrete(limits=rev(levels(s$Population)))
# #ggsave(filename='SFS_lof.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Facet by PGx 
rsid=fread("~/path/to/crv.freeze1/id.maps/id1.idE.txt",header=T)
pharm=fread("~/path/to/page/CONTENT_0.CRV/PHARMGKB/final.pharmgkb.annotated.rsid.txt",header=F)
head(pharm)
pgx=page[page$SNP.AA %in% pharm$V2,]
pgx=merge(pgx, rsid, by.x="SNP.AA",by.y="mine",all.x=TRUE)
pharm=read.table('~/Dropbox/page/Fst/pgx/1_annotate/out.merged-page-pharmgkb-table.tsv',header=T,sep='\t',fill = T)
pharm=pharm[,c('rsID','SNP.AA','Gene','Drug','Evidence','Disease','Type')]
pgx=merge(pgx,pharm,by.x='rsID',by.y='rsID')
pgx=unique(pharm$SNP.AA)


# All PGx
m=ddply(dfm[dfm$SNP %in% pgx,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% pgx,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
d=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('PGx variants (%)') +
  theme(axis.title.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_PGx.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)


# Highevidence
strong=unique(pharm[pharm$Evidence %in% c("1A","1B"),]$SNP.AA)
 
m=ddply(dfm[dfm$SNP %in% strong,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% strong,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
e=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('PGx Level1 variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_PGx1A1B.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)

adme=read.csv('~/Dropbox/page/ADME_ref_panel/adme_panel.csv',header=T)
pgx=merge(page,rsid,by.x='SNP.AA',by.y='mine')
adme.var=pgx[pgx$rsID %in% adme$ASSAY,]$SNP.AA

#ADME
m=ddply(dfm[dfm$SNP %in% adme.var,], 'Population', transform, c=cut (MAF, breaks =c (0.000001, 0.01,0.05,1)))
m$c=as.character(m$c)
m[m$Counts==1,]$c='singleton'
m[m$MAF==0,]$c='absent'
s=ddply(m, c('Population','c'), plyr::summarize, total_content=length(SNP), percent_content=round(length(SNP)/length(unique(dfm[dfm$SNP %in% adme.var,]$SNP)) * 100, digits=1) )
s$c=factor(s$c,levels=rev(c("absent","singleton","(1e-06,0.01]", "(0.01,0.05]","(0.05,1]")),
           labels = rev(c('absent','singleton','rare: <1%','frequent: 1-5%','common: >5%')),ordered=T)
s$Population=factor(s$Population)
f=ggplot(s, aes(x=Population,y=percent_content,fill=c))+geom_bar(stat='identity',position='stack',color='black',lwd=0.2)+theme_classic() +
  scale_fill_brewer(palette='BrBG',direction=-1)+
  theme(axis.text.x=element_text(vjust=0.6))+ylab('ADME PGx variants (%)') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +coord_flip()+
  guides(fill=F)+
  scale_x_discrete(limits=rev(levels(s$Population)))
#ggsave(filename='SFS_adme.png',plot = last_plot(), width=2,height=2,units='in',dpi=300)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

pdf('Summary.SFS.pdf',width=6,height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 10)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(a, vp = vplayout(1, 1:4))  # key is to define vplayout
print(b, vp = vplayout(1, 5:7))  
print(c, vp = vplayout(1, 8:10))  
print(d, vp = vplayout(2, 1:4))  
print(e, vp = vplayout(2, 5:7))  
print(f, vp = vplayout(2, 8:10))  
print(g, vp = vplayout(3, 1:4))  
print(h, vp = vplayout(3, 5:7)) 
print(i, vp = vplayout(3, 8:10)) 
dev.off()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Summary --ClinVar
above=unique(clinvar[,c('bimid','clinical_significance')])
names(above)=c('SNP','clinical_significance')
above=unique(above[above$SNP %in% page$SNP.AA,])
above$clinical_significance=as.character(above$clinical_significance)
# Conflicted 
conflicted=unique(above[above$clinical_significance=="Conflicted",]$SNP)
above[above$SNP %in% conflicted,]$clinical_significance="Conflicted"  # Get rid of other assertions for Conflicted
above=above[!duplicated(above),]
above[above$clinical_significance %in% c("Benign","Likely benign"),]$clinical_significance= "Benign/LB"

above=above[!duplicated(above),]
# Variants with multiple assertions should be conflicted, except if with B/LB not provided
above=ddply(above,'SNP',transform,reps=length(SNP))
dups=above[above$reps>1,]
d=ddply(dups,'SNP',plyr::summarize,cat=paste(clinical_significance[1],clinical_significance[2]))
conflicted=d[!d$cat %in% c("Benign/LB not provided","not provided Benign/LB"),]$SNP
above[above$SNP %in% conflicted,]$clinical_significance="Conflicted"  # Get rid of other assertions for Conflicted

benign=d[d$cat %in% c("Benign/LB not provided","not provided Benign/LB"),]$SNP
above[above$SNP %in% benign,]$clinical_significance="Benign/LB"  # Get rid of other assertions for Conflicted
above=above[,c(1:3)]
above=above[!duplicated(above),]

above[!above$clinical_significance %in% c ("Pathogenic","Likely pathogenic","Conflicted",
                                           "VUS","not provided","Benign/LB"),]$clinical_significance='Other'

above$clinical_significance=mapvalues(above$clinical_significance,from=c("Pathogenic","Likely pathogenic","Conflicted","VUS","not provided","Benign/LB","Other"), 
                                      to=c('P','LP','C','VUS','np','B/LB','oth'))
above$clinical_significance=factor(above$clinical_significance,levels=rev(c('P','LP','C','VUS','np','B/LB','oth'),ordered=T)

above=na.omit(above)
dim(above) # Should be about 23K variants 

a=ddply(above,'clinical_significance',plyr::summarize,total= length(SNP), pct= length(SNP)/dim(above)[1]* 100 )
aa=ggplot(a,aes(x=clinical_significance,y=pct))+geom_bar(fill='grey60',stat='identity')+theme_classic()+ylab('% Clinical Assertions')+xlab('')
aa
ggsave(plot=last_plot(),filename='ClinVar-assertions.png',dpi=330,height=2,width=2,units='in')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#Summary -- Fxnal 
m=unique(cadd[,c('bimid','Consequence')])
m=cadd[cadd$bimid %in% page$SNP.AA,]
# Make some adjustments
m$Consequence=as.character(m$Consequence)
regulatory=c('3PRIME_UTR','5PRIME_UTR','DOWNSTREAM','INTERGENIC','NONCODING_CHANGE','REGULATORY','UPSTREAM','INTRONIC')
splice=c('CANONICAL_SPLICE','SPLICE_SITE')
m[m$Consequence %in% regulatory,]$Consequence='Non-coding'
m[m$Consequence %in% splice,]$Consequence='Splicing'
m[m$Consequence=="NON_SYNONYMOUS",]$Consequence='Missense'
m[m$Consequence=="SYNONYMOUS",]$Consequence='Synonymous'
m[m$Consequence=="STOP_GAINED",]$Consequence='Stop gain'
m[m$Consequence=="FRAME_SHIFT",]$Consequence='Frameshift'

# Take only the largest classes with at least 100 variants
small=c('INFRAME','STOP_LOST')
m=m[!m$Consequence %in% small,]
s=ddply(m, 'Consequence', plyr::summarize, total_content=length(bimid), percent_content=length(bimid)/dim(m)[1] * 100 )
s=s[order(-s$total_content),]
s$Consequence=factor(s$Consequence,levels=s$Consequence, ordered=T)
ggplot(s, aes(x=factor(Consequence),y=percent_content))+geom_bar(stat='identity',position='identity',fill='grey60')+
  geom_text(aes(label=total_content),vjust=-0.2,size=2)+
  theme_classic() + theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,size=6))+ylab('CRVs in PAGE (%)')
ggsave(plot=last_plot(),filename='Consequences.png',dpi=330,height=2,width=2,units='in')
