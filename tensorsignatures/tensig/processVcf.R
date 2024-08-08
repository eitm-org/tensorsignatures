#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(VariantAnnotation)
library(rhdf5)
library(rtracklayer)
# PLEASE adjust paths
source("/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/mutations.R")
load("/home/xchen@okta-oci.eitm.org/projects/tensorsignatures/tensig/constants.RData")

genomePath <- "/data/xchen/refs/GRCh38/GRCh38.primary_assembly.genome_X.fa" #"genome.fa.gz" # path to the reference genome

liftOverAnnots <- function(gr){
  chain_path = '~/projects/tensorsignatures/gh38_granges/hg19ToHg38.over.chain'
  ch = rtracklayer::import.chain(chain_path)
  seqlevelsStyle(gr) = "UCSC"  # necessary
  gr = liftOver(gr, ch)
  unlist(gr)
}


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}
readVcfSave <- function(path) {
    # may need some user modifcation
    vcf <- VariantAnnotation::readVcf(path)
    vcf <- vcf[seqnames(vcf) %in% paste0('chr', c(1:22, "X", "Y"))] # filt(vcf) == "PASS" &
    elementLengths <- elementNROWS(alt(vcf))
    vcf <- vcf[elementLengths==1]
}

processVcf <- function(vcf, ts, epi, nuc, rt) {
    
    # extracts the trinucleotide context of each single base substitution
    tnc <- getTrinucleotideContext(vcf, genomePath)
    # annotates the mutation substitution
    sub <- getTrinucleotideSubs(vcf, tnc)
    # extracts transcription directionality for each substitution
    ts <- getStrandOrientation(vcf, ts)
    # extracts replication directionality for each substitution
    rs <- getStrandOrientation(vcf, rt)
    # extracts epigenetic state for each substitution
    ep <- getChromatinState(vcf, epi)
    # extracts the nucleosome directionality for each substitution
    nu <- getNucleosomeState(vcf, nuc)
    # extracts clustering state
    cl <- getClustering(vcf)
    t <- table(ts=ts, rs=rs, ep=ep, nu=nu, cl=cl, sub=sub)
    t[,,,,,SUB]
}

# output file is specified in the last argument
outputPath <- args[length(args)]
# path = '/data/scratch/xchen/STATE_vcfs_filter3_funcotated_region_filtered/chrY/EIBS-00218_09421_1_1_1_20230714_vs_HG002_06119_15_20230714.ann.filtered3.chrom.vcf.gz'
# vcf = readVcfSave(path)
inputDir = '/data/scratch/xchen/STATE_vcfs_f3_region_filtered_funcotated'
inputDir = args[1]

# inputPaths = list.files(inputDir, pattern='*.chrom.vcf$', full.names = T, recursive = TRUE)
inputPaths = list.files(inputDir, pattern='*.filtered[1,3].vcf$', full.names = T, recursive = TRUE)
# print(length(inputPaths))

inputPaths = sort(list.files(inputDir, pattern=args[2], full.names = T, recursive = TRUE))
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

TS38 = liftOverAnnots(TS)
EPI38 = liftOverAnnots(EPI)
NUC38 = liftOverAnnots(NUC)
RT38 = liftOverAnnots(RT)
c = 1
for (inputPaths_chunk in chunk(inputPaths, 3)){
  print(inputPaths_chunk)
  vcf <- lapply(inputPaths_chunk, readVcfSave)
  snvTensor <- sapply(vcf, function(x) processVcf(x[isSNV(x)], TS38, EPI38, NUC38, RT38), simplify="array")
  indelTable <- sapply(vcf, function(x) getIndels(x[isIndel(x)]), simplify="array")
  mnvTable <- sapply(vcf, function(x) getMNV(x), simplify="array")

  print("Trying to save ...")
  h5createFile(outputPath)
  print("Writing")
  h5write(snvTensor, "SNV", file=paste0(outputPath, '_chunk', c, '.h5'))
  h5write(indelTable, "INDELS", file=paste0(outputPath, '_chunk', c, '.h5'))
  h5write(mnvTable, "MNV", file=paste0(outputPath, '_chunk', c, '.h5'))
  c = c + 1
}