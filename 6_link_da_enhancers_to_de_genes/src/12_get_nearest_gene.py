import os
import pybedtools

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
pybedtools.helpers.set_tempdir(CURRENT_DIR_PATH)


def get_closest_gene(master_file, gtf_file, save_file):
    a = pybedtools.BedTool(master_file)
    # get the closest feature in 'other.bed' on the same strand
    gtf_bed = pybedtools.BedTool(gtf_file).sort()
    b = a.closest(gtf_bed, d=True)
    b.moveto(save_file)
    return


if __name__ == "__main__":
    master_file = "/data5/deepro/starrseq/papers/results/0_compare_activity_ko_vs_wt/data/window/IN/master_filtered_windows.bed"
    gtf_file = "/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/genome_annot/parsed_gtf.tsv"
    save_file = "/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/predictions/nearest_genes/closest.bed"

    get_closest_gene(master_file, gtf_file, save_file)
    pybedtools.helpers.cleanup(verbose=False, remove_all=False)
