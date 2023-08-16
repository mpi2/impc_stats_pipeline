import argparse
import csv
import os.path

def chunking(input_file, output_dir_path):
    csvreader = csv.DictReader(open(input_file))

    a = "procedure_stable_id"
    b = "parameter_stable_id"
    c = "phenotyping_center"
    d = "pipeline_stable_id"
    e = "strain_accession_id"
    f = "metadata_group"
    g = "biological_sample_group"

    for row in csvreader:
        filename = "_".join([row[a], row[b], row[c], row[d], row[e], row[f], row[g]])
        filename = os.path.join(output_dir_path, filename + ".csv")
        with open(filename, mode='a') as outfile:
            out_writer = csv.DictWriter(outfile, fieldnames=csvreader.fieldnames)
            if outfile.tell() == 0:
                 out_writer.writeheader()
            out_writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="File input and output path parser")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_dir_path", help="Path to the output directory")
    args = parser.parse_args()
    chunking(args.input_file, args.output_dir_path)
