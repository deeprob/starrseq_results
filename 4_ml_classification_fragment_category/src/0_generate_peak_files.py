import argparse
import utils as ut

CHROM_SIZES = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chrX": 156040895,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr11": 135086622,
    "chr10": 133797422,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr20": 64444167,
    "chr19": 58617616,
    "chrY": 57227415,
    "chr22": 50818468,
    "chr21": 46709983
}

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Generate peak files for STARR-seq data")
    parser.add_argument("peak_file", type=str,
                        help="Path to peak file")
    parser.add_argument("master_file", type=str, 
                        help="Path to master filtered bed file")
    parser.add_argument("save_dir", type=str,
                        help="Directory to save output files")
    parser.add_argument("--window_size", type=int, default=1000,
                        help="Window size for peak generation")
    parser.add_argument("--overlap_frac", type=float, default=0.5,
                        help="Overlap fraction for peak generation")
    
    args = parser.parse_args()
    
    ut.generate_peak_files(args.peak_file, args.window_size, CHROM_SIZES, 
                          args.save_dir, overlap_frac=args.overlap_frac, 
                          target_file=args.master_file)

