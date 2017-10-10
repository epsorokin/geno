require(data.table)
require(stringr)
require(plyr)

## Read in SNP, HGDP Ref Freq, Consequence type

ann=fread('/path/to/Data.Csq.MAF.MEGA.tsv',header=T)
ann=data.frame(ann)
ann=ann[,c("SNP","MAF.HGDP_Eur","Consequence")]

## Read in 
# Add in PGx category
pgx=fread('/path/to/mega.pgx-ids',sep='\t',header=F)

ann$pgx='No'
ann[ann$SNP %in% pgx$V1,]$pgx='Yes'

# Add MAF bin
ann$bucket=cut(ann$MAF.HGDP_Eur,breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

# Modify Consequence for complex variants
temp=ann[grep('&',ann$Consequence),]
temp=ddply(temp,'Consequence',transform, count=str_count(Consequence,'&')+1)
temp=ddply(temp,'Consequence',transform,Consequence2=sample( str_split_fixed (Consequence,'&',count) ))
temp=temp[,c('SNP','Consequence2','bucket','pgx')]
names(temp)[2]='Consequence'
not=ann[!ann$SNP %in% temp$SNP,]
not=not[,c('SNP','Consequence','bucket','pgx')]
d=rbind(temp,not)

# We can remove ones where MAF is zero in HGDP Eur
d=d[!is.na(d$bucket),]
# Also remove duplicates
d=d[!duplicated(d),]
# Separate into p and notp
p=d[d$pgx=='Yes',]
notp=d[d$pgx=='No',]

# Read in Fst data
fst=fread('Fst.HGDP_Eur_v_PAGE.V2.tsv',header=T)
fst=data.frame(fst)
fst=fst[,c('SNP','Fstnum','Fstdenom','Group')]

# Now lets calc Fst for PGx vars

pharma=fst[fst$SNP %in% pgx$V1,]
result= ddply(pharma, 'Group', summarize, Fst=sum(Fstnum,na.rm=T)/sum(Fstdenom,na.rm=T))
result
write.table(result,'out.RatioOfAvgs.PGx.tsv',sep='\t',quote=F,row.names=F)
savehistory('merge_annotations.Fst-Pgx.r')
