library("poolfstat")

allpops.readcount <-vcf2pooldata(vcf.file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_all.vcf.gz",poolsizes=rep(80,335),min.maf="0.01")
setwd("/media/DD8T_1/DSU/Data_lewald/res_vcf") #directory containing the vcf files
require(poolfstat)
##recup info: det.all=nom des pools ou individus
## det.pools=true pools (with name extension _tp) + fake pools (with name extension _fp)
true.pools=c(26, 290, 291, 292, 293, 294, 295, 296, 329, 330, 331, 332, 333, 334, 335)
fake.pools.DGN=c(127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289)
fake.pools.CNXJ=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
fake.pools.CnOther=c(27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
fake.pools.CnQTP=c(77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)
fake.pools.ISR=c(297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328)
tmp.nom=read.table("/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles_2.txt",stringsAsFactors = F)
tmp.nom=matrix(unlist(strsplit(as.character(tmp.nom[,1]),split="\\.")),ncol=2,byrow=T)[,1]
tmp.nomtp=data.frame(nom=paste0(tmp.nom[true.pools],"_tp"),pop=tmp.nom[true.pools], fullid=tmp.nom[true.pools], hapsize=100,type="TruePool",stringsAsFactors = FALSE)
tmp.ind.DGN=cbind(tmp.nom[fake.pools.DGN],matrix(unlist(strsplit(tmp.nom[fake.pools.DGN],split=",")), ncol=2,byrow=T))
tmp.ind.CNXJ=cbind(tmp.nom[fake.pools.CNXJ],matrix(unlist(strsplit(tmp.nom[fake.pools.CNXJ],split=",")), ncol=2,byrow=T))
tmp.ind.CnOther=cbind(tmp.nom[fake.pools.CnOther],matrix(unlist(strsplit(tmp.nom[fake.pools.CnOther],split=",")), ncol=2,byrow=T))
tmp.ind.CnQTP=cbind(tmp.nom[fake.pools.CnQTP],matrix(unlist(strsplit(tmp.nom[fake.pools.CnQTP],split=",")), ncol=2,byrow=T))
tmp.ind.ISR=cbind(tmp.nom[fake.pools.ISR],matrix(unlist(strsplit(tmp.nom[fake.pools.ISR],split=",")), ncol=2,byrow=T))
colnames(tmp.ind.DGN)=c("FullID","POP","ID")
colnames(tmp.ind.CNXJ)=c("FullID","POP","ID")
colnames(tmp.ind.CnOther)=c("FullID","POP","ID")
colnames(tmp.ind.CnQTP)=c("FullID","POP","ID")
colnames(tmp.ind.ISR)=c("FullID","POP","ID")

dum.DGN=as.matrix(table(tmp.ind.DGN[,2]))
tmp.nomfp.DGN=data.frame(nom=paste0(rownames(dum.DGN),"_fp"),pop=rownames(dum.DGN),hapsize=2*dum.DGN, type="FakePool",stringsAsFactors = F)
det.pools.DGN=rbind(tmp.nomtp,tmp.nomfp.DGN)
det.all.DGN=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.DGN[,3]), pop=c(tmp.nom[true.pools],tmp.ind.DGN[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.DGN))), stringsAsFactors = F)

dum.CNXJ=as.matrix(table(tmp.ind.CNXJ[,2]))
tmp.nomfp.CNXJ=data.frame(nom=paste0(rownames(dum.CNXJ),"_fp"),pop=rownames(dum.CNXJ),hapsize=2*dum.CNXJ, type="FakePool",stringsAsFactors = F)
det.pools.CNXJ=rbind(tmp.nomtp,tmp.nomfp.CNXJ)
det.all.CNXJ=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CNXJ[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CNXJ[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CNXJ))), stringsAsFactors = F)

dum.CnOther=as.matrix(table(tmp.ind.CnOther[,2]))
tmp.nomfp.CnOther=data.frame(nom=paste0(rownames(dum.CnOther),"_fp"),pop=rownames(dum.CnOther),hapsize=2*dum.CnOther, type="FakePool",stringsAsFactors = F)
det.pools.CnOther=rbind(tmp.nomtp,tmp.nomfp.CnOther)
det.all.CnOther=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CnOther[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CnOther[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CnOther))), stringsAsFactors = F)

dum.CnQTP=as.matrix(table(tmp.ind.CnQTP[,2]))
tmp.nomfp.CnQTP=data.frame(nom=paste0(rownames(dum.CnQTP),"_fp"),pop=rownames(dum.CnQTP),hapsize=2*dum.CnQTP, type="FakePool",stringsAsFactors = F)
det.pools.CnQTP=rbind(tmp.nomtp,tmp.nomfp.CnQTP)
det.all.CnQTP=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.CnQTP[,3]), pop=c(tmp.nom[true.pools],tmp.ind.CnQTP[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.CnQTP))), stringsAsFactors = F)

dum.ISR=as.matrix(table(tmp.ind.ISR[,2]))
tmp.nomfp.ISR=data.frame(nom=paste0(rownames(dum.ISR),"_fp"),pop=rownames(dum.ISR),hapsize=2*dum.ISR, type="FakePool",stringsAsFactors = F)
det.pools.ISR=rbind(tmp.nomtp,tmp.nomfp.ISR)
det.all.ISR=data.frame(fullid=c(tmp.nomtp$fullid,tmp.ind.ISR[,3]), pop=c(tmp.nom[true.pools],tmp.ind.ISR[,2]), seqtype=c(rep("PoolSeq",length(true.pools)),rep("IndSeq",nrow(tmp.ind.ISR))), stringsAsFactors = F)

rm(dum,tmp.nom,tmp.nomfp,tmp.ind,tmp.nomtp)

#############
fake.pools.idx.DGN=which(det.pools.DGN$type=="FakePool")
nfake.pools.DGN=length(fake.pools.idx)
fakepool.ind.idx.DGN=list()
for(i in 1:nfake.pools){
tmp.nom.DGN=det.pools.DGN$pop[fake.pools.idx.DGN[i]]
fakepool.ind.idx.DGN[[i]]=which(det.all.DGN$pop==tmp.nom.DGN & det.all.DGN$seqtype=="IndSeq")
}
#######
vcf.files=list.files(path="./results/",pattern="bcftools_all.vcf.gz")
nctg.gg=0 ; nctg.pp=0
for(f in vcf.files){#loop over each contig vcf
tmp<-vcf2pooldata(f,poolsizes = rep(10,287),nlines.per.readblock = 1e4,min.maf=0.05)

#creation de fake pool
nctg.pp=nctg.pp+1
tmp.refcount=tmp@refallele.readcount[,pools.pos.idx]
tmp.readcov=tmp@readcoverage[,pools.pos.idx]
for(i in 1:nfake.pools){
tmp.refcount=cbind(tmp.refcount,rowSums(tmp@refallele.readcount[,fakepool.ind.idx[[i]]]))
tmp.readcov=cbind(tmp.readcov,rowSums(tmp@readcoverage[,fakepool.ind.idx[[i]]]))
}
if(nctg.pp==1){
all.pools=tmp
all.pools@refallele.readcount=tmp.refcount
all.pools@readcoverage=tmp.readcov
all.pools@npools=nrow(det.pools)
all.pools@poolsizes=det.pools$hapsize
all.pools@poolnames=det.pools$nom
}else{
all.pools@refallele.readcount=rbind(all.pools@refallele.readcount,tmp.refcount)
all.pools@readcoverage=rbind(all.pools@readcoverage,tmp.readcov)
all.pools@nsnp=all.pools@nsnp+tmp@nsnp
all.pools@snp.info=rbind(all.pools@snp.info,tmp@snp.info)
}
}
rm(tmp)
}