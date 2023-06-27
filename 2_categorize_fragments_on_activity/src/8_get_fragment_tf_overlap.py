import os
import pandas as pd
import pybedtools

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath("__file__"))
pybedtools.helpers.set_tempdir(CURRENT_DIR_PATH)

def get_chip_files(chip_dir):
    tf_chip_files = []
    tf_chip_names = []

    for tf_dir in os.scandir(chip_dir):
        tf_name = tf_dir.name
        for encode_dir in os.scandir(tf_dir.path):
            for encode_file in os.scandir(encode_dir.path):
                if encode_file.name.endswith("bed.gz"):
                    tf_chip_files.append(encode_file.path)
                    tf_chip_names.append(tf_name)
    return tf_chip_files, tf_chip_names


if __name__ == "__main__":
    chip_dir = "/data5/deepro/starrseq/papers/reproducibility/0_in-house_dataset/data/encode/hek293/chip/tf"
    master_fragments = "/data5/deepro/starrseq/papers/results/1_compare_activity_ko_vs_wt/data/window/IN/master_filtered_windows.bed"
    save_bed = "/data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/data/hot_sites/tf_overlap_fragments.bed"

    tf_chip_files, tf_chip_names = get_chip_files(chip_dir)
    master_bed = pybedtools.BedTool(master_fragments) 
    tf_overlap = master_bed.intersect(a=master_fragments, b=tf_chip_files, f=0.5, wa=True, wb=True, names=tf_chip_names).to_dataframe(disable_auto_names=True, header=None).loc[:, :3]
    tf_overlap.columns = ["chrm", "start", "end", "tf"]
    tf_overlap = tf_overlap.groupby(["chrm", "start", "end"]).aggregate(lambda x: ",".join(x)).reset_index()
    tf_overlap["number"] = tf_overlap.tf.str.split(",").apply(len)
    tf_overlap.to_csv(save_bed, index=False, header=None, sep="\t")
    pybedtools.helpers.cleanup(remove_all=True)
