import pandas as pd
import requests
import os
import pysam


#############
# parse gtf #
#############

def parse_gtf_attributes(attr):
    attr_dict = dict()
    for iattr in attr.split(";"):
        if iattr:
            k,v = iattr.strip().split()
            attr_dict[k] = v.strip('"')
    return pd.Series(data=attr_dict)

def parse_gtf_file(gtf_file, save_file):
    # read gtf
    gtf_df = pd.read_csv(
        gtf_file, 
        sep="\t", 
        comment="#", 
        header=None,
        names=["chrm", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        )
    # only keep genes
    gtf_df = gtf_df.loc[gtf_df.feature=="gene"]
    # parse the attributes column to get gene id and gene name
    gtf_df = gtf_df.merge(gtf_df.attribute.apply(parse_gtf_attributes), left_index=True, right_index=True)
    gtf_df = gtf_df.loc[:, ["chrm", "start", "end", "gene_id", "gene_name", "strand"]].reset_index(drop=True)
    gtf_df.to_csv(save_file, sep="\t", index=False, header=False)
    return


###################
# encode download #
###################

def get_search_dict(url):
    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}
    # GET the search result
    response = requests.get(url, headers=headers)
    # Extract the JSON response as a python dictionary
    search_results = response.json()
    return search_results


def download_file(metadata, storage_dir, filename):
    url = metadata["url"]
    filename = os.path.join(storage_dir, filename)
    if os.path.exists(filename):
        # TODO: change to logger
        print(f"Warning: File {filename} already exists")
    file = open(filename, "wb")
    response = requests.get(url)
    assert response.status_code == 200
    file.write(response.content)
    return filename


def get_file_from_encode(acc, storage_dir, filebase):
    search_file_url = f"https://www.encodeproject.org/files/{acc}"
    search_file_dict = get_search_dict(search_file_url)
    cloud_metadata = search_file_dict["cloud_metadata"]
    os.makedirs(storage_dir,  exist_ok=True)
    download_file(cloud_metadata, storage_dir, filebase)
    return


####################
# expression table #
####################

def create_expression_table_helper(counts_file, gtf_file, store_file):
    expr_df = pd.read_csv(counts_file, sep="\t", engine="python", skipfooter=5).set_index(["gene_id", "gene_name"])
    expr_df["mean_raw_counts"] = expr_df.mean(axis=1)
    gtf_df = pd.read_csv(gtf_file, sep="\t", header=None, names=["chrm", "start", "end", "gene_id", "gene_name", "strand"]).set_index(["gene_id", "gene_name"])
    expr_df = expr_df.merge(gtf_df, left_index=True, right_index=True)
    assert len(expr_df) == len(gtf_df)
    expr_df["length"] = expr_df.end - expr_df.start
    expr_df["expression"] = (expr_df.mean_raw_counts*1e3)/expr_df.length
    expr_df["expression"] = (expr_df["expression"]*1e6)/expr_df["expression"].sum()
    expr_df = expr_df.reset_index()
    expr_df.loc[:, ["gene_id", "expression"]].to_csv(store_file, sep="\t", index=False, header=False)
    return


##################
# activity table #
##################

def get_num_reads_bam_pair1(bam_file):
    # val = pysam.view("-c", bam_file)
    # read mapped in proper pair and first in pair :: https://broadinstitute.github.io/picard/explain-flags.html
    val = pysam.view("-f", "66", "-c", bam_file)
    return int(val.strip())

def get_rpm_normalized_depth_df(depth_file, bam_file):
    depth_df = pd.read_csv(depth_file, sep="\t", header=None)
    depth_df = depth_df.iloc[:, [0, 1, 2, -4]]
    depth_df.columns = ["chr", "start", "end", "reads"]
    bam_reads_pair = get_num_reads_bam_pair1(bam_file)
    depth_df["RPM"] = (depth_df.reads*1e6)/bam_reads_pair
    depth_df["activity_base"] = depth_df["RPM"]
    return depth_df

def create_activity_table_helper(depth_file, bam_file, enhancer_file, store_file):
    depth_df = get_rpm_normalized_depth_df(depth_file, bam_file)
    enhancer_df = pd.read_csv(enhancer_file, sep="\t", header=0, usecols=['chr','start','end','name','class'])
    activity_df = depth_df.merge(enhancer_df, left_on=["chr", "start", "end"], right_on=["chr", "start", "end"]) 
    activity_df.to_csv(store_file, sep="\t", index=False, header=True)
    return
