import os
import argparse
import pickle
import pycisTopic
import pandas
from pycisTopic.lda_models import run_cgs_models

def main():
    parser = argparse.ArgumentParser(description="Run pycisTopic model.")
    parser.add_argument("--input_object", type=str, required=True, help="Path to the input cistopic object.")
    parser.add_argument("--output_model", type=str, required=True, help="Path to the output model to save.")
    parser.add_argument("--nTopics", type=int, required=True, help="Number of topics.")
    parser.add_argument("--nIter", type=int, required=True, help="Number of iterations.")

    args = parser.parse_args()

    with open(args.input_object, 'rb') as infile:
        cistopic_object = pickle.load(infile)
    
    model_dir = os.path.dirname(args.output_model)
    os.makedirs(model_dir, exist_ok=True)

    models = run_cgs_models(cistopic_object, n_topics=[args.nTopics], n_iter=args.nIter,
                            n_cpu=5, _temp_dir=os.path.abspath("/tmp"),
                            save_path=model_dir)

if __name__ == "__main__":
    main()
