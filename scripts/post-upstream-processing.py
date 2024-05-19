#!/usr/bin/env python

import argparse
import pandas as pd

def process_count_data(metadata_path, counts_data_path, output_path):
    # load count data and metadata
    metadata_df = pd.read_csv(metadata_path)
    count_data_df = pd.read_csv(counts_data_path, sep="\t")

    # rename columns
    count_data_df.columns = ["geneSymbol"] + [name.split("/")[-1].split(".")[0] for name in count_data_df.columns[1:]]

    # Assume 'sampleAccession' column in metadata contains the sample names in the desired order
    desired_order = metadata_df['sampleAccession'].tolist()

    # Reorder the columns of count data DataFrame according to the desired order
    count_data_df = count_data_df.reindex(columns=["geneSymbol"]+desired_order)

    # Remove duplicates based on GeneSymbol column
    count_data_df = count_data_df.drop_duplicates(keep='first')
    count_data_df = count_data_df.drop_duplicates(subset='geneSymbol', keep='first')

    # Remove missing values
    count_data_df = count_data_df.dropna().reset_index(drop=True)

    # Remove genes with zero values across all samples
    count_data_df = count_data_df.loc[:, (count_data_df != 0).any(axis=0)]

    # Remove genes with negative values
    count_data_df = count_data_df.loc[(count_data_df.iloc[:,1:] >= 0).all(axis=1), :]

    # save count data
    count_data_df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process count data and metadata")
    parser.add_argument("metadata", help="Path to the metadata CSV file")
    parser.add_argument("counts_data", help="Path to the counts data TXT file")
    parser.add_argument("output", help="Path to save the processed count data CSV file")
    args = parser.parse_args()

    process_count_data(args.metadata, args.counts_data, args.output)

if __name__ == "__main__":
    main()
