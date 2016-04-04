gene.groups.url <- 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_group.gz'

if ( ! file.exists('Homo_sapiens.gene_groups.gz') ) {
  require('curl')
  message('Downloading file...')
  curl_download(gene.groups.url, 'gene_groups.gz')
  system('gzip -dc gene_groups.gz | grep "^9606" | gzip - > Homo_sapiens.gene_groups.gz')
} else {
  message('Using previously downloaded Homo_sapiens.gene_groups.gz, file age is ',
          format(Sys.time() - file.mtime('Homo_sapiens.gene_groups.gz')))
}

gene.groups.cols <- c('taxid', 'id', 'relationship', 'ptaxid', 'partner')
all.related.genes <- read.table('Homo_sapiens.gene_groups.gz', sep='\t', col.names=gene.groups.cols,
                                stringsAsFactors=FALSE)
related.genes <- subset(all.related.genes, ptaxid=='9606')[, c('id', 'relationship', 'partner')]

save(related.genes, file='related_genes.RData')
