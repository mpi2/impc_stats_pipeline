import argparse
import csv
import math
import os
import os.path
import zipfile

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

    csv_control_list = [row for row in csv.DictReader(open(file_ctrl))]
    n_controls = len(csv_control_list)
    csv_experiment = csv.DictReader(open(file_exp))
    data_dict = {}
    elem_name = os.path.split(file_ctrl)[1].split("_")

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
            outfile_name = "_".join(elem_name[:-1] + [zygosity, str(count)])
            outfile_name_csv =  os.path.join(output_dir_path, outfile_name + ".csv")

            with open(outfile_name_csv, mode='w') as outfile:
                out_writer = csv.DictWriter(outfile, fieldnames=csv_experiment.fieldnames)
                out_writer.writeheader()

                # Write control to file
                for row_ctrl in csv_control_list:
                    out_writer.writerow(row_ctrl)
                # Write experimental to file
                for colony in chunk:
                    for row_expr in colony:
                        out_writer.writerow(row_expr)
            outfile_name_zip = os.path.join(output_dir_path, outfile_name + ".zip")
            
            # Compress file with ZIP and remove original CSV file
            with zipfile.ZipFile(outfile_name_zip, "w") as zipf:
                zipf.write(outfile_name_csv, arcname = os.path.split(outfile_name_csv)[1])
            os.remove(outfile_name_csv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="File control and experimental samples parser")
    parser.add_argument("file_ctrl", help="Path to the file with control samples")
    parser.add_argument("file_exp", help="Path to the file with experimental samples")
    parser.add_argument("output_dir_path", help="Path to the output directory")
    args = parser.parse_args()

    divide_chunk(args.file_ctrl, args.file_exp, args.output_dir_path)
