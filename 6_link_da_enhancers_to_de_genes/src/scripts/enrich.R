# library("biomaRt")
library("clusterProfiler")
library("enrichplot")
library(org.Hs.eg.db)


args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=6) {
  stop("gene_list, gsea table, kegg table, gsea figure, kegg figure filenames and tmp_dir must be given", call.=FALSE)
}

gene_file = args[1]
gseaout_file = args[2]
keggout_file = args[3]
gseafigout_file = args[4]
keggfigout_file = args[5]
tmp_dir = args[6]

# get the genes and store as a vector
genes = read.table(gene_file, header=FALSE)

# Sys.setenv(BIOMART_CACHE = tmp_dir)
# print(biomartCacheInfo())

# # convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
# # mart <- useMart("ensembl","hsapiens_gene_ensembl")
# mart = useEnsembl(biomart='ensembl', dataset=ensembl_dataset, mirror = "uswest")
# entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), "hgnc_symbol", genes, mart)

# use enrichGO for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html

df = as.data.frame(org.Hs.egSYMBOL)
go_gene_list = unique(sort(df$symbol))

goenrich = enrichGO(
    gene=genes$V1,
    universe=go_gene_list,
    keyType = "SYMBOL",
    OrgDb='org.Hs.eg.db',
    pAdjustMethod="BH",
    pvalueCutoff=0.05,
    ont="ALL"
)

# save to file .. 
write.table(goenrich, file=gseaout_file, sep=",", row.names=TRUE, col.names=TRUE)

pdf(gseafigout_file)
myplot <- dotplot(goenrich, showCategory=20, font.size=8) + ggtitle("dotplot for GO enrichment")
print(myplot)
dev.off()



# # use kegg for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
# kenrich = enrichKEGG(
#     gene=genes,
#     organism="hsa",
#     keyType = "SYMBOL",
#     pvalueCutoff=0.05,
#     pAdjustMethod="BH"
# )

# print(length(kenrich))

# if (length(kenrich)>1){
#     pdf(keggfigout_file)
#     myplot <- dotplot(kenrich, showCategory=25, font.size=8) + ggtitle("dotplot for KEGG")
#     print(myplot)
#     dev.off()
#     # save to file .. 
#     write.table(kenrich, file=keggout_file, sep=",", row.names=TRUE, col.names=TRUE)
# }


