require(poolfstat)
##recup info: det.all=nom des pools ou individus
## det.pools=true pools (with name extension _tp) + fake pools (with name extension _fp)
true.pools=c(26, 290, 291, 292, 293, 294, 295, 296, 329, 330, 331, 332, 333, 334, 335)
fake.pools.CnQTP=c(77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)
tmp.nom=read.table("/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles_2.txt",stringsAsFactors = F)
tmp.nom=matrix(unlist(strsplit(as.character(tmp.nom[,1]),split="\\.")),ncol=2,byrow=T)[,1]
tmp.nomtp=data.frame(nom=paste0(tmp.nom[true.pools],"_tp"),pop=tmp.nom[true.pools], fullid=tmp.nom[true.pools], hapsize=100,type="TruePool",stringsAsFactors = FALSE)
tmp.ind.CnQTP=cbind(tmp.nom[fake.pools.CnQTP],matrix(unlist(strsplit(tmp.nom[fake.pools.CnQTP],split=",")), ncol=2,byrow=T))
colnames(tmp.ind.CnQTP)=c("FullID","POP","ID")


dum.CnQTP=as.matrix(table(tmp.ind.CnQTP[,2]))
tmp.nomfp.CnQTP=data.frame(nom=paste0(rownames(dum.CnQTP),"_fp"),pop=rownames(dum.CnQTP),hapsize=2*dum.CnQTP, type="FakePool",stringsAsFactors = F)
det.pools.CnQTP=rbind(tmp.nomtp,tmp.nomfp.CnQTP)
det.all.CnQTP=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CnQTP[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CnQTP[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CnQTP))), stringsAsFactors = F)

rm(dum.CnQTP,tmp.nom,tmp.nomfp.CnQTP,tmp.ind.CnQTP,tmp.nomtp)

#############
fake.pools.idx.CnQTP=which(det.pools.CnQTP$type=="FakePool")
nfake.pools.CnQTP=length(fake.pools.idx)
fakepool.ind.idx.CnQTP=list()
for(i in 1:nfake.pools){
tmp.nom.CnQTP=det.pools.CnQTP$pop[fake.pools.idx.CnQTP[i]]
fakepool.ind.idx.CnQTP[[i]]=which(det.all.CnQTP$pop==tmp.nom.CnQTP & det.all.CnQTP$seqtype=="IndSeq")
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
tmp.refcount=cbind(tmp.refcount,rowSums(tmp@randomallele.pca[,fakepool.ind.idx.CnQTP[[i]]]))
tmp.readcov=cbind(tmp.readcov,rowSums(tmp@readcoverage[,fakepool.ind.idx.CnQTP[[i]]]))
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
