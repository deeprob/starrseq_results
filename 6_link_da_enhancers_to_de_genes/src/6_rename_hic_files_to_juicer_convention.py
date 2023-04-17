# manually rename files according to the juicer convention
import os
import re
import logging


def rename_file(file):
    if file.name.endswith(".fastq"):
        dirname = os.path.dirname(file.path)
        basename = file.name
        newname = re.sub(r'_([1-2])', r'_R\1', basename)
        oldname = os.path.join(dirname, basename)
        newname = os.path.join(dirname, newname)
        logging.info(f"Renaming {oldname} to {newname}!")
        os.rename(oldname, newname)
    return


if __name__ == "__main__":
    renaming_dir = "/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/data/hic/hek293t/fastq"
    logging.basicConfig(filename=os.path.join(renaming_dir, "file_rename.log"), encoding='utf-8', level=logging.INFO)
    for file in os.scandir(renaming_dir):
        rename_file(file)
