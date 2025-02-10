require(poolfstat)

# script.R
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    stop("At least one argument must be supplied", call.=FALSE)
}

# Access arguments by index
vcf_file <- args[1]
writing_dir <- args[2]

vcfdata <- vcf2pooldata(vcf.file=vcf_file, min.maf="0.05", poolsizes=c(40,50,80,80,80))

pooldata2diyabc(
    vcfdata,
    writing.dir = writing_dir,
    prefix = "NoInvData",
    diyabc.mrc = 5,
)


