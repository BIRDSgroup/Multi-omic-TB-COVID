library("dada2")
for(f_name in c("batch1","batch2","batch3","batch4","batch5"))
{
    path <- paste0(f_name,"/")
    print(path)

    fnFs <- sort(list.files(path, pattern="_R1_P.fastq.gz", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern="_R2_P.fastq.gz", full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
   

    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    print("filter_done")
     
    #trunbeg set trimLeft = 0,
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=0,
                           maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                           compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

    dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
    seqtab <- makeSequenceTable(mergers)
    dim(seqtab)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE,minFoldParentOverAbundance=1)
    dim(seqtab.nochim)
    sum(seqtab.nochim)/sum(seqtab)
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    head(track)
    save.image(paste0(path,"workspace_mfP_default.RData"))

    write.table(t(seqtab.nochim), paste0(path,"seqtab-nochim_mFP_default.txt"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
    uniquesToFasta(seqtab.nochim, fout=paste0(path,"rep-seqs_mFP_default.fna"), ids=colnames(seqtab.nochim))
}

