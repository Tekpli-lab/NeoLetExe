import os
import argparse
import pickle
import pycisTopic
import scipy
import pandas
from pycisTopic.cistopic_class import create_cistopic_object

def convert_region_ids(region_ids):
    converted_ids = []
    for i in range(len(region_ids)):
        tmp = region_ids[i].split("-")
        converted_ids.append(f"{tmp[0]}:{tmp[1]}-{tmp[2]}")
    return(converted_ids)


def main():
    parser = argparse.ArgumentParser(description="Create pycisTopic object.")
    parser.add_argument("--input_counts_dir", type=str, required=True, help="Path to the input count matrix (.mtx file).")
    parser.add_argument("--input_metadata", type=str, required=True, help="Path to the input metadata CSV.")
    parser.add_argument("--output_pkl", type=str, required=True, help="Path to the output .pkl file.")

    args = parser.parse_args()

    # Read input files
    counts = scipy.io.mmread(open(os.path.join(args.input_counts_dir, "count_mat.mtx"), "r")).tocsr()
    cell_ids = open(os.path.join(args.input_counts_dir, "colnames.txt"), "r").read().splitlines()
    region_ids = open(os.path.join(args.input_counts_dir, "rownames.txt"), "r").read().splitlines()
    metadata = pandas.read_csv(args.input_metadata, index_col=0)
    metadata = metadata.loc[cell_ids]

    cistopic_object = create_cistopic_object(fragment_matrix=counts, cell_names=cell_ids, region_names=convert_region_ids(region_ids))
    cistopic_object.add_cell_data(metadata)
    cistopic_object.cell_data['patient_number'] = cistopic_object.cell_data['patient_number'].astype(str)

    # Save the cistopic object using the output path specified in the arguments
    pickle.dump(cistopic_object, open(args.output_pkl, "wb"))

if __name__ == "__main__":
    main()
