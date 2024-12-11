require(poolfstat)
##recup info: det.all=nom des pools ou individus
## det.pools=true pools (with name extension _tp) + fake pools (with name extension _fp)
true.pools=c(26, 290, 291, 292, 293, 294, 295, 296, 329, 330, 331, 332, 333, 334, 335)
fake.pools.CNXJ=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
tmp.nom=read.table("/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles_2.txt",stringsAsFactors = F)
tmp.nom=matrix(unlist(strsplit(as.character(tmp.nom[,1]),split="\\.")),ncol=2,byrow=T)[,1]
tmp.nomtp=data.frame(nom=paste0(tmp.nom[true.pools],"_tp"),pop=tmp.nom[true.pools], fullid=tmp.nom[true.pools], hapsize=100,type="TruePool",stringsAsFactors = FALSE)
tmp.ind.CNXJ=cbind(tmp.nom[fake.pools.CNXJ],matrix(unlist(strsplit(tmp.nom[fake.pools.CNXJ],split=",")), ncol=2,byrow=T))
colnames(tmp.ind.CNXJ)=c("FullID","POP","ID")


dum.CNXJ=as.matrix(table(tmp.ind.CNXJ[,2]))
tmp.nomfp.CNXJ=data.frame(nom=paste0(rownames(dum.CNXJ),"_fp"),pop=rownames(dum.CNXJ),hapsize=2*dum.CNXJ, type="FakePool",stringsAsFactors = F)
det.pools.CNXJ=rbind(tmp.nomtp,tmp.nomfp.CNXJ)
det.all.CNXJ=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CNXJ[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CNXJ[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CNXJ))), stringsAsFactors = F)

rm(dum.CNXJ,tmp.nom,tmp.nomfp.CNXJ,tmp.ind.CNXJ,tmp.nomtp)

#############
fake.pools.idx.CNXJ=which(det.pools.CNXJ$type=="FakePool")
nfake.pools.CNXJ=length(fake.pools.idx)
fakepool.ind.idx.CNXJ=list()
for(i in 1:nfake.pools){
tmp.nom.CNXJ=det.pools.CNXJ$pop[fake.pools.idx.CNXJ[i]]
fakepool.ind.idx.CNXJ[[i]]=which(det.all.CNXJ$pop==tmp.nom.CNXJ & det.all.CNXJ$seqtype=="IndSeq")
}
#######
vcf.files=list.files(path="./results/",pattern="bcftools_all.vcf.gz")
nctg.gg=0 ; nctg.pp=0
for(f in vcf.files){#loop over each contig vcf
tmp<-vcf2pooldata(f,poolsizes = rep(10,287),nlines.per.readblock = 1e4,min.maf=0.05)

#creation de fake pool
nctg.pp=nctg.pp+1
tmp.refcount=tmp@randomallele.pca[,true.pools]
tmp.readcov=tmp@readcoverage[,true.pools]
for(i in 1:nfake.pools){
tmp.refcount=cbind(tmp.refcount,rowSums(tmp@randomallele.pca[,fakepool.ind.idx.CNXJ[[i]]]))
tmp.readcov=cbind(tmp.readcov,rowSums(tmp@readcoverage[,fakepool.ind.idx.CNXJ[[i]]]))
}
if(nctg.pp==1){
all.pools=tmp
all.pools@randomallele.pca=tmp.refcount
all.pools@readcoverage=tmp.readcov
all.pools@npools=nrow(det.pools)
all.pools@poolsizes=det.pools$hapsize
all.pools@poolnames=det.pools$nom
}else{
all.pools@randomallele.pca=rbind(all.pools@randomallele.pca,tmp.refcount)
all.pools@readcoverage=rbind(all.pools@readcoverage,tmp.readcov)
all.pools@nsnp=all.pools@nsnp+tmp@nsnp
all.pools@snp.info=rbind(all.pools@snp.info,tmp@snp.info)
}
}
rm(tmp)
