import os
import subprocess
import pandas as pd
import time


CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

########################
# gsea kegg enrichment #
########################

def run_enrichment_helper(gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir):
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/scripts/enrich.sh", 
        gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir
        ]
    subprocess.run(cmd)
    return

def run_enrichment(de_genes_dir, lib, store_dir, pre):
    table_dir = os.path.join(store_dir, lib, "tables")
    os.makedirs(table_dir, exist_ok=True)
    fig_dir = os.path.join(store_dir, lib, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    tmp_dir = os.path.join(store_dir, lib, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    gene_in = os.path.join(de_genes_dir, lib, f'{pre}_sig_genes.txt')
    gseaout_file = os.path.join(table_dir, f'{pre}_gsea.csv')
    keggout_file = os.path.join(table_dir, f'{pre}_kegg.csv')
    gseafigout_file = os.path.join(fig_dir, f'{pre}_gsea.pdf')
    keggfigout_file = os.path.join(fig_dir, f'{pre}_kegg.pdf')
    run_enrichment_helper(gene_in, gseaout_file, keggout_file, gseafigout_file, keggfigout_file, tmp_dir)
    return

def save_array_to_text(arr, save_file):
    with open(save_file, "w") as f:
        for gene in arr:
            f.write(f"{gene}\n")
    return

def read_tf_df(filename, libn):
    colnames = ["chrom", "start", "end", "chrom_coord", f"{libn}_act", "strand", f"{libn}_log2FoldChange_act", f"gene_name", f"{libn}_log2FoldChange_exp", f"{libn}_padj_act", f"{libn}_padj_exp", "CC_peak", f"{libn}_peak"]
    return pd.read_csv(filename, sep="\t", header=None, names=colnames)


if __name__=="__main__":
    sda_sde_table_dir = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks"
    libraries = ["ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    for lib_name in libraries:
        dl_lib_file = os.path.join(sda_sde_table_dir, lib_name, "direct_loss.bed")
        dl_lib_df = read_tf_df(dl_lib_file, lib_name)
        idg_lib_file = os.path.join(sda_sde_table_dir, lib_name, "indirect_gained.bed")
        idg_lib_df = read_tf_df(idg_lib_file, lib_name)

        # direct loss enrichment
        sig_gene_file = os.path.join(sda_sde_table_dir, lib_name, f'direct_loss_sig_genes.txt')
        save_array_to_text(dl_lib_df.gene_name.unique(), sig_gene_file)
        run_enrichment(sda_sde_table_dir, lib_name, sda_sde_table_dir, 'direct_loss')

        # indirect gained enrichment
        sig_gene_file = os.path.join(sda_sde_table_dir, lib_name, f'indirect_gained_sig_genes.txt')
        save_array_to_text(idg_lib_df.gene_name.unique(), sig_gene_file)
        run_enrichment(sda_sde_table_dir, lib_name, sda_sde_table_dir, 'indirect_gained')
