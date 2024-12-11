require(poolfstat)
##recup info: det.all=nom des pools ou individus
## det.pools=true pools (with name extension _tp) + fake pools (with name extension _fp)
true.pools=c(26, 290, 291, 292, 293, 294, 295, 296, 329, 330, 331, 332, 333, 334, 335)
fake.pools.ISR=c(297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328)
tmp.nom=read.table("/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles_2.txt",stringsAsFactors = F)
tmp.nom=matrix(unlist(strsplit(as.character(tmp.nom[,1]),split="\\.")),ncol=2,byrow=T)[,1]
tmp.nomtp=data.frame(nom=paste0(tmp.nom[true.pools],"_tp"),pop=tmp.nom[true.pools], fullid=tmp.nom[true.pools], hapsize=100,type="TruePool",stringsAsFactors = FALSE)
tmp.ind.ISR=cbind(tmp.nom[fake.pools.ISR],matrix(unlist(strsplit(tmp.nom[fake.pools.ISR],split=",")), ncol=2,byrow=T))
colnames(tmp.ind.ISR)=c("FullID","POP","ID")


dum.ISR=as.matrix(table(tmp.ind.ISR[,2]))
tmp.nomfp.ISR=data.frame(nom=paste0(rownames(dum.ISR),"_fp"),pop=rownames(dum.ISR),hapsize=2*dum.ISR, type="FakePool",stringsAsFactors = F)
det.pools.ISR=rbind(tmp.nomtp,tmp.nomfp.ISR)
det.all.ISR=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.ISR[,3]), pop=c(tmp.nom[true.pools],tmp.ind.ISR[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.ISR))), stringsAsFactors = F)

rm(dum.ISR,tmp.nom,tmp.nomfp.ISR,tmp.ind.ISR,tmp.nomtp)

#############
fake.pools.idx.ISR=which(det.pools.ISR$type=="FakePool")
nfake.pools.ISR=length(fake.pools.idx)
fakepool.ind.idx.ISR=list()
for(i in 1:nfake.pools){
tmp.nom.ISR=det.pools.ISR$pop[fake.pools.idx.ISR[i]]
fakepool.ind.idx.ISR[[i]]=which(det.all.ISR$pop==tmp.nom.ISR & det.all.ISR$seqtype=="IndSeq")
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
tmp.refcount=cbind(tmp.refcount,rowSums(tmp@randomallele.pca[,fakepool.ind.idx.ISR[[i]]]))
tmp.readcov=cbind(tmp.readcov,rowSums(tmp@readcoverage[,fakepool.ind.idx.ISR[[i]]]))
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
