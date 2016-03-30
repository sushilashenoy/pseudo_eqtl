

gene2refseq.url <- 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz'
gene2refseq.cols <- c('taxid', 'entrez', 'status', 'refseq', 'refseq.gi', 'prot_acc', 'prot.gi',
                      'genome.acc', 'genome.gi', 'start', 'end', 'orientation', 'assembly',
                      'peptide.acc', 'peptide.gi', 'symbol')

if ( ! file.exists('gene2refseq.gz') ) {
  require('curl')
  curl_download(gene2refseq.url, 'gene2refseq.gz')
  system('gzip -dc gene2refseq.gz | grep "^9606" | gzip - > Homo_sapiens.gene2refseq.gz')
} else {
  message('Using previously downloaded gene2refseq.gz, file age is ',
          format(Sys.time() - file.mtime('gene2refseq.gz')))
}

gene2refseq <- read.table('Homo_sapiens.gene2refseq.gz', sep='\t', col.names=gene2refseq.cols)


# Identify the most used assembly and ignore others
use.assembly <- names(which.max(table(gene2refseq$assembly)))
entrez2refseq <- subset(gene2refseq, assembly==use.assembly)

# Save table for looking up entrez genes from refseq IDs
refseq2entrez <- unique(entrez2refseq[, c('refseq', 'entrez')])
save(refseq2entrez, file='refseq2entrez.RData')

# Identify entries corresponding to full chromosomes
is.genomic <- substring(entrez2refseq$genome.acc, 1, 2)=='NC'
# Get rid of some columns we don't really care about
ignore.cols <- c('taxid', 'refseq.gi', 'prot_acc', 'prot.gi', 'genome.gi',
                 'assembly', 'peptide.acc', 'peptide.gi')
keep.cols <- !colnames(entrez2refseq) %in% ignore.cols
entrez2genome <- subset(entrez2refseq, is.genomic)[, keep.cols]

# Replace chromosome refseq accessions with chromosome identifiers
chrom.acc <- sort(unique(entrez2genome$genome.acc))
chromosomes <- c(1:22, 'X', 'Y')
entrez2genome$genome.acc <- chromosomes[match(entrez2genome$genome.acc, chrom.acc)]
colnames(entrez2genome)[colnames(entrez2genome)=='genome.acc'] <- 'chromosome'

# There is duplicate info for many genes, some are real duplicates, others are
# alternative transcripts that map to the same genome location, the latter can
# be removed

# First rank duplicate entries based on status
status.rank <- c('VALIDATED', 'REVIEWED', 'PROVISIONAL', 'PREDICTED', 'MODEL', 'INFERRED', 'SUPPRESSED')
entrez2genome <- entrez2genome[with(entrez2genome, order(entrez, match(status, status.rank))), ]

# Remove duplicates where entrez/chromosome/start/end/orient are identical
remove.duplicate <- with(entrez2genome, duplicated(paste(entrez, chromosome, start, end, orientation)))
entrez2genome <- subset(entrez2genome, !remove.duplicate)

# Now identify & deal with remaining duplicates
duplicate.idx <- which(duplicated(entrez2genome$entrez))
original.idx <- match(entrez2genome$entrez[duplicate.idx], entrez2genome$entrez)

# If the status of the duplicate doesn't match, keep the original only.
status.mismatch <- entrez2genome$status[duplicate.idx] != entrez2genome$status[original.idx]
if ( any(status.mismatch) )
  entrez2genome <- entrez2genome[-duplicate.idx[status.mismatch], ]

# Remaining duplicates are genes which are actually duplicated on the genome,
# either on X & Y chromosome and others have proximal duplicates on the same
# autoosome. For now leaving them in the table, but they will be ignored when
# we use match() to locate an entrez gene ID.


save.cols <- c('entrez', 'chromosome', 'start', 'end', 'orientation', 'symbol')
entrez2genome <- entrez2genome[, save.cols]

save(entrez2genome, file='entrez2genome.RData')
