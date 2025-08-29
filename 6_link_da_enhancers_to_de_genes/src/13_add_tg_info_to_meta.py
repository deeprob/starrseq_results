import pandas as pd

def add_info(meta_activity_map, abc_file, nearest_file, save_file):
    # parse activity map and only store required info
    meta_activity_df = pd.read_csv(meta_activity_map)
    libraries = ["CC", "ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    relevant_columns = ["chrom_coord"] + libraries
    lib_peaks = [f"{lib}_peak" for lib in libraries]
    lib_fc = [f"{lib}_log2FoldChange" for lib in libraries[1:]]
    lib_pval = [f"{lib}_pvalue" for lib in libraries[1:]]
    lib_padj = [f"{lib}_padj" for lib in libraries[1:]]
    relevant_columns = relevant_columns + lib_peaks + lib_fc + lib_pval + lib_padj
    meta_activity_df = meta_activity_df.loc[:, relevant_columns]

    # parse abc file
    abc_df = pd.read_csv(abc_file, sep="\t")
    # select strict threshold
    strict_abc_df = abc_df.loc[abc_df["ABC.Score"]>0.1]
    # add chrom coord
    strict_abc_df["chrom_coord"] = strict_abc_df.chr + "_" + strict_abc_df.start.astype(str) + "_" + strict_abc_df.end.astype(str)
    # aggregate target gene info
    strict_abc_df = strict_abc_df.groupby("chrom_coord").agg({"TargetGene": lambda x: "|".join(x)})
    strict_abc_df = strict_abc_df.rename(columns={"TargetGene": "abc_target"})

    # parse nearest gene file
    nearest_df = pd.read_csv(nearest_file, sep="\t", header=None)
    # feature must be within 5000 bp
    nearest_df = nearest_df.loc[nearest_df[9]<=5000]
    # rename columns
    nearest_df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "nearest_target", "nearest_gene_symbol", "strand", "distance"]
    # add chrom coord
    nearest_df["chrom_coord"] = nearest_df.chr1 + "_" + nearest_df.start1.astype(str) + "_" + nearest_df.end1.astype(str)
    # aggregate target gene info
    nearest_df = nearest_df.groupby("chrom_coord").agg({"nearest_target": lambda x: "|".join(x)})

    # add info to activity map
    meta_activity_df["abc_gene"] = meta_activity_df.chrom_coord.map(strict_abc_df.squeeze().to_dict())
    meta_activity_df["nearest_gene"] = meta_activity_df.chrom_coord.map(nearest_df.squeeze().to_dict())

    # save activity map
    meta_activity_df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    meta_activity_map = "/data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/data/meta_activity_map.csv"
    abc_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/predictions/activity_from_starr_with_hic/EnhancerPredictions.txt"
    nearest_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/predictions/nearest_genes/closest.bed"
    save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/meta_enhancer_gene.csv"

    add_info(meta_activity_map, abc_file, nearest_file, save_file)
