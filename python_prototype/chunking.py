import argparse
import csv
import math
import os
import os.path
import zipfile


def chunking(input_file, output_dir_path):
    csvreader = csv.DictReader(open(input_file, "rt"))

    a = "procedure_stable_id"
    b = "parameter_stable_id"
    c = "phenotyping_center"
    d = "pipeline_stable_id"
    e = "strain_accession_id"
    f = "metadata_group"
    g = "biological_sample_group"

    all_chunks = {}
    for row in csvreader:
        filename = "_".join([row[a], row[b], row[c], row[d], row[e], row[f], row[g]])
        all_chunks.add(filename)
        filename = os.path.join(output_dir_path, filename + ".csv")
        with open(filename, mode='a') as outfile:
            out_writer = csv.DictWriter(outfile, fieldnames=csvreader.fieldnames)
            if outfile.tell() == 0:
                 out_writer.writeheader()
            out_writer.writerow(row)
    return all_chunks


def divide_chunk(file_ctrl,
                 file_exp,
                 output_dir_path,
                 control_size=1500,
                 n_max=10000,
                 min_colonies_in_chunks=32,
                 chunk_size=24):
    
    if not (os.path.isfile(file_ctrl) and os.path.isfile(file_exp)):
        print("Either control or experimental file do not exist")
        return

    control_data = open(file_ctrl).read()
    n_controls = control_data.count('\n') - 1
    csv_experiment = csv.DictReader(open(file_exp))
    data_dict = {}
    elem_name = os.path.split(file_ctrl)[1].split("_")[:-1]

    for row in csv_experiment:
        zyg = row["zygosity"]
        col_id = row["colony_id"]
        
        data_dict.setdefault(zyg, {})
        data_dict[zyg].setdefault(col_id, [])
        data_dict[zyg][col_id].append(row)
        
    for zygosity, colonies in data_dict.items():
        n_colonies = len(colonies)
        if n_controls < control_size:
            chunks = 1
        elif control_size <= n_controls < n_max:
            if n_colonies < min_colonies_in_chunks:
                chunks = 1
            elif n_colonies >= min_colonies_in_chunks:
                chunks = round(n_colonies / chunk_size)
        elif n_controls >= n_max:
            chunks = n_colonies
                   
        exp_chunks = list(colonies.values())
        chunk_size = math.ceil(len(colonies) / chunks)
        
        chunks_list = [exp_chunks[i:i + chunk_size] for i in range(0, len(colonies), chunk_size)]
        for count, chunk in enumerate(chunks_list):
            outfile_basename = "_".join(elem_name + [zygosity, str(count)])
            outfile_csv = os.path.join(output_dir_path, outfile_basename + ".csv")
            outfile_zip = os.path.join(output_dir_path, outfile_basename + ".zip")

            with open(outfile_csv, mode='w') as outfile:
                # Write experimental to file
                out_writer = csv.DictWriter(outfile, fieldnames=csv_experiment.fieldnames)
                out_writer.writeheader()
                for colony in chunk:
                    for row_expr in colony:
                        out_writer.writerow(row_expr)

            # Compress file with ZIP and remove temporary CSV file
            with zipfile.ZipFile(outfile_zip, "w", compression=zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(outfile_csv, arcname = outfile_basename + ".csv")
            os.remove(outfile_csv)

        # Finally, compress the control file once
        outfile_basename = "_".join(elem_name + [zygosity, "control"])
        outfile_zip = os.path.join(output_dir_path, outfile_basename + ".zip")
        with zipfile.ZipFile(outfile_zip, "w", compression=zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(file_ctrl, arcname = outfile_basename + ".csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="File input and output path parser")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("intermediate_dir_path", help="Path to the intermediate output directory")
    parser.add_argument("output_dir_path", help="Path to the final output directory")
    args = parser.parse_args()
    all_chunks = chunking(args.input_file, args.intermediate_dir_path)
    for chunk in all_chunks:
        file_ctrl = os.path.join(args.intermediate_dir_path, chunk) + "_control.csv"
        file_exp = os.path.join(args.intermediate_dir_path, chunk) + "_experimental.csv"
        divide_chunk(file_ctrl, file_exp, args.output_dir_path)
