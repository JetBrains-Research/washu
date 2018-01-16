files <- list.files()

geneResults <- files[grepl("\\.genes\\.results", files)]
rnames <- as.character(read.table(geneResults[1], sep="\t", header=1, stringsAsFactors=F)[, 1])

sampleName <- function(resultsFile) {
        bname <- basename(resultsFile)
            woExt <- strsplit(bname, "\\.")[[1]][1]
                return(woExt)
}

sampleNames <- sapply(geneResults, sampleName)

allCounts <- lapply(geneResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 5, drop=F]
})
allCounts <- do.call(cbind, allCounts)
rownames(allCounts) <- rnames
colnames(allCounts) <- sampleNames



allTpms <- lapply(geneResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 6, drop=F]
})
allTpms <- do.call(cbind, allTpms)
rownames(allTpms) <- rnames
colnames(allTpms) <- sampleNames



allFpkms <- lapply(geneResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 7, drop=F]
})
allFpkms <- do.call(cbind, allFpkms)
rownames(allFpkms) <- rnames
colnames(allFpkms) <- sampleNames

write.table(allCounts, "genome_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(allTpms, "genome_tpms.tsv", sep="\t", quote=F, col.names=NA)
write.table(allFpkms, "genome_fpkms.tsv", sep="\t", quote=F, col.names=NA)


isoformResults <- files[grepl("\\.isoforms\\.results", files)]
rnames <- as.character(read.table(isoformResults[1], sep="\t", header=1, stringsAsFactors=F)[, 1])
sampleNames <- sapply(isoformResults, sampleName)


allCounts <- lapply(isoformResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 5, drop=F]
})
allCounts <- do.call(cbind, allCounts)
rownames(allCounts) <- rnames
colnames(allCounts) <- sampleNames



allTpms <- lapply(isoformResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 6, drop=F]
})
allTpms <- do.call(cbind, allTpms)
rownames(allTpms) <- rnames
colnames(allTpms) <- sampleNames



allFpkms <- lapply(isoformResults, function(geneFile) {
        read.table(geneFile, sep="\t", header=1, stringsAsFactors=F)[, 7, drop=F]
})
allFpkms <- do.call(cbind, allFpkms)
rownames(allFpkms) <- rnames
colnames(allFpkms) <- sampleNames

write.table(allCounts, "transcriptome_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(allTpms, "transcriptome_tpms.tsv", sep="\t", quote=F, col.names=NA)
write.table(allFpkms, "transcriptome_fpkms.tsv", sep="\t", quote=F, col.names=NA)


