#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(VariantAnnotation)
library(rhdf5)
source("/usr/src/app/data/R/mutations.R")
load("/usr/src/app/data/R/constants.RData")
genomePath <- "/usr/src/app/data/R/genome.fa.gz" # path to the reference genome
genomePath <- "/data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X.fa" #"genome.fa.gz" # path to the reference genome

# load functions and ranges please adjust!
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}
print("INPUT")
print(args)

outputPath <- args[length(args)] # where the output tensor is saved to

print(outputPath)

readVcfSave <- function(path) {
    vcf <- readVcf(path)
    vcf <- vcf[seqnames(vcf) %in% c(1:22, "X", "Y")] # filt(vcf) == "PASS" &
    elementLengths <- elementNROWS(alt(vcf))
    vcf <- vcf[elementLengths==1]
}

processVcf <- function(vcf) {
    print("Getting Trinculeotide context")
    tnc <- getTrinucleotideContext(vcf, genomePath)
    print("Getting Subs")
    sub <- getTrinucleotideSubs(vcf, tnc)
    print("Getting Trx")
    ts <- getStrandOrientation(vcf, TS)
    print("Getting Rep")
    rs <- getStrandOrientation(vcf, RT)
    print("Getting epi")
    ep <- getChromatinState(vcf, EPI)
    print("Getting miuc")
    nu <- getNucleosomeState(vcf, NUC)
    print("Getting clu")
    cl <- getClustering(vcf)
    print("Getting tab")
    t <- table(ts=ts, rs=rs, ep=ep, nu=nu, cl=cl, sub=sub)
    t[,,,,,SUB]
}

df2vcf <- function(df) {
  c <- df$chr
  s <- df$pos
  e <- as.integer(s+sapply(as.character(df$ref), nchar, simplify="array")-1)
  g <- GRanges(c, IRanges(s, e), "*")
  genome(g) <- "hg19"
  v <- VCF(rowRanges=g)
  ref(v) <- DNAStringSet(df$ref)
  alt(v) <- DNAStringSetList(lapply(df$alt, function(x) x))
  toNCBI(v)
}

print("loading vcf ...")
vcf <- lapply(args[1:length(args)-1], readVcfSave)

print("processing ...")
snvTensor <- sapply(vcf, function(x) processVcf(x[isSNV(x)]), simplify="array")
indelTable <- sapply(vcf, function(x) getIndels(x[isIndel(x)]), simplify="array")
mnvTable <- sapply(vcf, function(x) getMNV(x), simplify="array")

print("Trying to save ...")
h5createFile(outputPath)

print("Writing")
h5write(snvTensor, "SNVR", file=outputPath)
h5write(indelTable, "INDELS", file=outputPath)
h5write(mnvTable, "MNV", file=outputPath)
