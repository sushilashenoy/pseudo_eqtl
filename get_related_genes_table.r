gene.groups.url <- 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_group.gz'

tmp <- tempfile()
require('curl')
curl_download(gene.groups.url, tmp)

all.related.genes <- read.table(tmp, sep='\t', col.names=c('taxid', 'id', 'relationship', 'ptaxid', 'partner'))

related.genes <- subset(all.related.genes, taxid=='9606' & ptaxid=='9606')[, c('id', 'relationship', 'partner')]

save(related.genes, file='related_genes.RData')
