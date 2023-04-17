import argparse
import utils as ut


def main(
    acc,
    storage_dir,
    filename
    ):
    ut.get_file_from_encode(acc, storage_dir,  filename)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="download files from encode given their accession number")
    parser.add_argument("accession", type=str, help="The encode file accession number for the file to be downloaded")
    parser.add_argument("storage_dir", type=str, help="dir to store results")
    parser.add_argument("filename", type=str, default="", help="Enter file basename to save the file")
    
    args = parser.parse_args()
    main(args.accession, args.storage_dir, args.filename)
