import os
import numpy as np
import pandas as pd
import pybedtools
from pybedtools.featurefuncs import greater_than


def make_windows_faster(in_bed, window_size=500, window_stride=50):
    """
    Break the ROIs into fragments of an user defined window size and stride
    """
    window = pybedtools.BedTool().window_maker(b=in_bed, w=window_size, s=window_stride)
    # only keep windows of length greater than w-s-1
    window_new_filtered = window.filter(greater_than, window_size - window_stride + 1).saveas()
    window_df = window_new_filtered.to_dataframe()
    return window_df

def create_centered_peaks(bed_file, window_size, chrom_sizes):
    """
    Create peaks centered around given peak coordinates with a user-defined window size.
    Ensures peaks remain exactly `window_size` bp long, adjusting within chromosome limits.
    Parameters:
    - bed_file (str): Path to the narrowPeak BED file.
    - window_size (int): Desired window size (minimum 500 bp).
    - chrom_sizes (dict): Dictionary of chromosome sizes.
    Returns:
    - DataFrame with fixed-size centered peaks.
    """
    if window_size < 500:
        raise ValueError("Window size cannot be less than 500 bp.")
    # Load the BED file
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(bed_file, sep="\t", header=None, names=cols, usecols=[0, 1, 2])
    # Filter out invalid chromosomes early
    df = df[df["chrom"].isin(chrom_sizes)]
    # Determine offset from center
    df_peaks = pd.DataFrame()
    for position, offset in zip(["center", "left", "right"], [window_size//2, window_size//4, 3 * window_size//4]):
        # Compute initial start and end
        peak_center = (df["start"] + df["end"]) // 2
        new_start = peak_center - offset 
        new_end = new_start + window_size
        # Get chromosome max size
        chrom_max_sizes = df["chrom"].map(chrom_sizes)
        # Shift if out of bounds
        shift_right = np.where(new_start < 1, -new_start+1, 0)  # Extra bp to shift right
        shift_left = np.where(new_end > chrom_max_sizes, new_end - chrom_max_sizes, 0)  # Extra bp to shift left
        # Adjust start and end while keeping size fixed
        new_start = new_start + shift_right - shift_left
        new_end = new_end + shift_right - shift_left
        assert np.all(new_end-new_start==window_size)
        # Create final dataframe
        df_centered = pd.DataFrame({
            "chrom": df["chrom"],
            "start": new_start,
            "end": new_end
        })
        df_peaks = pd.concat([df_peaks, df_centered])
    return df_peaks

def generate_non_overlapping_windows(chrom_sizes, window_size=1000):
    """
    Generate non-overlapping genomic windows of fixed size for given chromosome sizes.
    Args:
    - chrom_sizes (dict): Dictionary mapping chromosome names to their lengths.
    - window_size (int): Length of each window (default=1000).
    Returns:
    - pd.DataFrame: DataFrame containing columns ["chrom", "start", "end"]
    """
    windows = []
    for chrom, size in chrom_sizes.items():
        # Generate start positions at intervals of `window_size`
        starts = np.arange(1, size - window_size + 1, window_size)
        ends = starts + window_size  # Compute end positions
        # If last end exceeds chromosome size, shift the last window
        if ends[-1] < size:
            last_end = size
            last_start = size - window_size
            starts = np.append(starts, [last_start])
            ends = np.append(ends, [last_end])
        # Store results as tuples (chrom, start, end)
        windows.extend(zip([chrom] * len(starts), starts, ends))
    return pd.DataFrame(windows, columns=["chrom", "start", "end"])

def generate_peak_files(peakfile_path, window_size, chrom_sizes, save_dir, overlap_frac=0.5, target_file=""):
    norm_bed_df = create_centered_peaks(peakfile_path, window_size, chrom_sizes)
    if not target_file:
        genome_df = generate_non_overlapping_windows(chrom_sizes, window_size)
    else:
        window_stride = int(0.1*window_size)
        genome_df = make_windows_faster(target_file, window_size, window_stride=window_stride)
    os.makedirs(save_dir, exist_ok=True)
    pybedtools.helpers.set_tempdir(save_dir)
    genome_bed = pybedtools.BedTool.from_dataframe(genome_df).sort()
    peak_bed = pybedtools.BedTool.from_dataframe(norm_bed_df).sort()
    non_peak_bed = genome_bed.intersect(peak_bed, v=True, F=overlap_frac)
    peak_path = peak_bed.saveas(os.path.join(save_dir, "peaks.bed.gz"))
    non_peak_path = non_peak_bed.saveas(os.path.join(save_dir, "not_peaks.bed.gz"))
    pybedtools.helpers.cleanup(verbose=False, remove_all=True)
    return peak_path, non_peak_path
