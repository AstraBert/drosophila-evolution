require(poolfstat)
##recup info: det.all=nom des pools ou individus
## det.pools=true pools (with name extension _tp) + fake pools (with name extension _fp)
true.pools=c(26, 290, 291, 292, 293, 294, 295, 296, 329, 330, 331, 332, 333, 334, 335)
fake.pools.CnOther=c(27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
tmp.nom=read.table("/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles_2.txt",stringsAsFactors = F)
tmp.nom=matrix(unlist(strsplit(as.character(tmp.nom[,1]),split="\\.")),ncol=2,byrow=T)[,1]
tmp.nomtp=data.frame(nom=paste0(tmp.nom[true.pools],"_tp"),pop=tmp.nom[true.pools], fullid=tmp.nom[true.pools], hapsize=100,type="TruePool",stringsAsFactors = FALSE)
tmp.ind.CnOther=cbind(tmp.nom[fake.pools.CnOther],matrix(unlist(strsplit(tmp.nom[fake.pools.CnOther],split=",")), ncol=2,byrow=T))
colnames(tmp.ind.CnOther)=c("FullID","POP","ID")

dum.CnOther=as.matrix(table(tmp.ind.CnOther[,2]))
tmp.nomfp.CnOther=data.frame(nom=paste0(rownames(dum.CnOther),"_fp"),pop=rownames(dum.CnOther),hapsize=2*dum.CnOther, type="FakePool",stringsAsFactors = F)
det.pools.CnOther=rbind(tmp.nomtp,tmp.nomfp.CnOther)
det.all.CnOther=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CnOther[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CnOther[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CnOther))), stringsAsFactors = F)

rm(dum.CnOther,tmp.nom,tmp.nomfp.CnOther,tmp.ind.CnOther,tmp.nomtp)

#############
fake.pools.idx.CnOther=which(det.pools.CnOther$type=="FakePool")
nfake.pools.CnOther=length(fake.pools.idx)
fakepool.ind.idx.CnOther=list()
for(i in 1:nfake.pools){
tmp.nom.CnOther=det.pools.CnOther$pop[fake.pools.idx.CnOther[i]]
fakepool.ind.idx.CnOther[[i]]=which(det.all.CnOther$pop==tmp.nom.CnOther & det.all.CnOther$seqtype=="IndSeq")
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
tmp.refcount=cbind(tmp.refcount,rowSums(tmp@randomallele.pca[,fakepool.ind.idx.CnOther[[i]]]))
tmp.readcov=cbind(tmp.readcov,rowSums(tmp@readcoverage[,fakepool.ind.idx.CnOther[[i]]]))
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
}