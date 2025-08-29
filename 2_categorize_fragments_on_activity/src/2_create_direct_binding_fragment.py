import os
import pandas as pd
import pybedtools
import argparse
import gzip

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath("__file__"))
pybedtools.helpers.set_tempdir(CURRENT_DIR_PATH)


def get_intersecting_fragments(responsive_file, chip_files, save_dir):
    responsive_bed = pybedtools.BedTool(responsive_file)
    tmp_chip_file = os.path.join(save_dir, "tmp_chip.bed")
    mode = "w"
    for file in chip_files:
        with open(tmp_chip_file, mode) as wf:
            with gzip.open(file, "rt") as rf:
                for lines in rf:
                    wf.write(lines)
        mode = "a"
    chip_bed = pybedtools.BedTool(tmp_chip_file).sort().merge()
    # regions in responsive fragments that intersect with chip
    intersect_bed = responsive_bed.intersect(chip_bed, u=True)
    non_intersect_bed = responsive_bed.intersect(chip_bed, v=True)
    intersect_file = os.path.join(save_dir, "direct_binding_responsive.bed")
    non_intersect_file = os.path.join(save_dir, "indirect_binding_responsive.bed")
    intersect_bed.moveto(intersect_file)
    non_intersect_bed.moveto(non_intersect_file)
    pybedtools.helpers.cleanup(remove_all=True)
    os.remove(tmp_chip_file)
    return

def get_chip_files(chip_dir, lib_name):
    chip_files = []
    chip_dir = os.path.join(chip_dir, lib_name)
    for folder in os.listdir(chip_dir):
        for file in os.listdir(os.path.join(chip_dir, folder)):
            if file.endswith("bed.gz"):
                chip_files.append(os.path.join(chip_dir, folder, file))
    return chip_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("response_dir", type=str, help="The response file")
    parser.add_argument("chip_dir", type=str, help="The dir where chip files are stored")
    parser.add_argument("store_dir", type=str, help="Save files dir")
    parser.add_argument("lib_name", type=str, help="The library name")

    cli_args = parser.parse_args()

    chip_files = get_chip_files(cli_args.chip_dir, cli_args.lib_name)
    response_file = os.path.join(cli_args.response_dir, cli_args.lib_name, "responsive.bed")
    save_dir = os.path.join(cli_args.store_dir, cli_args.lib_name)
    
    get_intersecting_fragments(
        response_file,
        chip_files, 
        save_dir
    )

