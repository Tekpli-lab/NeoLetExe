import os
import argparse
import pickle
import pycisTopic
from pycisTopic.clust_vis import harmony
import pandas
import scipy

parser = argparse.ArgumentParser(description="Run harmony and export to .mtx")
parser.add_argument("--output_dir", type=str)
parser.add_argument("--input_object", type=str)
parser.add_argument("--input_model", type=str)


args = parser.parse_args()

output_dir = args.output_dir
input_object = args.input_object
input_model = args.input_model


with open(input_object, "rb") as in_file:
    cistopic_object = pickle.load(in_file)

with open(input_model, "rb") as in_file:
    model = pickle.load(in_file)




cistopic_object.add_LDA_model(model)
harmony(cistopic_object, vars_use = ["orig.ident"])



with open(os.path.join(output_dir, "cell_names.txt"), "w") as cell_names:
    cell_names.writelines([s[:-11]+"\n" for s in model.cell_topic.columns.to_list()])

with open(os.path.join(output_dir, "dim_names.txt"), "w") as dim_names:
    dim_names.writelines([s+"\n" for s in model.cell_topic.index.to_list()])

scipy.io.mmwrite(os.path.join(output_dir, "mat.mtx"), scipy.sparse.coo_matrix(cistopic_object.selected_model.cell_topic_harmony))
