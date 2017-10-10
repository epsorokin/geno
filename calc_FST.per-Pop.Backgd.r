require(data.table)
require(stringr)
require(plyr)

## Read command line argument

args=commandArgs(trailingOnly=T)
e=args[1] # Ethnicity - first and only argument passed to script

## Read in SNP, HGDP Ref Freq, Consequence type, bucket, ADME/PGx categories
ann=data.frame(fread('out.Try1.Merge-ADME-annotations.tsv',header=T))

## Remove MAF=0 

ann=ann[!is.na(ann$bucket),]

## Simplify this data frame

ann=ann[,c('SNP','pgx','ADME','bucket','Consequence2')]
ann=ann[!duplicated(ann),]

## Separate ADME and non-ADME variants

p=ann[ann$ADME=='Core',] # Define SNPs within ADME - core bin
notp=ann[!ann$ADME %in% c('Core','Extended') & ann$pgx=="No",] # And the converse

## Now read in the Fst data
fst=data.frame(fread('/path/to/Fst.HGDP_Eur_v_PAGE.V2.tsv',header=T))
fst=fst[,c('SNP','Fstnum','Fstdenom','Group')]

# For each group; for each cons and bucket; sample equal grp of notp snps; append; list length of final thing 

result=data.frame(Group=as.character(),Fst=as.numeric())

for (i in 1:1000) {
print(i)
mysample=data.frame(SNP=as.character())

for (c in unique(na.omit(p$Consequence2))){
   for (f in unique(p$bucket)){
      vars=unique(p[p$Consequence2==c & p$bucket==f,]$SNP)
      backgd=data.frame(SNP=sample(notp[notp$Consequence2==c & notp$bucket==f,]$SNP, size=length(vars) )    )
      #print(length(vars));print(length(backgd))
      mysample=rbind(mysample,backgd)
      rm(backgd)
      }    
   }

background=fst[fst$SNP %in% mysample$SNP & fst$Group==e,]
result=rbind(result,ddply(background,'Group',summarize,Fst=sum(Fstnum,na.rm=T)/sum(Fstdenom,na.rm=T)))
rm(background) 
#print(length(mysample$SNP))

}

write.table(result,paste('results/out.BackgdFst.',e,'.Try1.tsv',sep=''),sep='\t',quote=F,row.names=F)



