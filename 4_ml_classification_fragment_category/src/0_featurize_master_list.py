import os
import argparse
import peaker as pk


def featurize_regions(roi_bed, genome_fasta, homer_pwm_motifs, homer_outdir, threads):
    hv = pk.HomerVectorizer(roi_bed, genome_fasta, homer_pwm_motifs, homer_outdir, threads)
    hv.featurize()
    homer_pickle = os.path.join(homer_outdir, "motif_features.pkl")
    hv.store_features_to_pickle(homer_pickle)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq ML classification")
    parser.add_argument("roi_bed", type=str, help="The roi bed file which contains all the regions that need to be featurized")
    parser.add_argument("homer_pwm_motifs", type=str, help="The motifs in homer format used to featurize the regions")
    parser.add_argument("genome_fasta", type=str, help="The path to the genome fasta file whose regions are given")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")    
    parser.add_argument("--threads", type=int, default=64, help="Number of cores to use")

    cli_args = parser.parse_args()

    featurize_regions(cli_args.roi_bed, cli_args.genome_fasta, cli_args.homer_pwm_motifs, cli_args.store_dir, cli_args.threads)
