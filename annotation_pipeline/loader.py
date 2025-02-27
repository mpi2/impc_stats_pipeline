import argparse
from pathlib import Path

from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr


# Activate pandas conversion.
pandas2ri.activate()

def main():
    parser = argparse.ArgumentParser(description="Process files and annotate them using R's annotationChooser.")
    parser.add_argument("file", help="Path to the file containing the list of files to process.")
    parser.add_argument("mp_chooser_file", help="Path to the mp_chooser.json.Rdata file.")
    args = parser.parse_args()

    file_list_path = Path(args.file)
    mp_chooser_file = args.mp_chooser_file

    # Load necessary R libraries.
    importr("base")
    importr("jsonlite")
    importr("utils")
    data_table = importr("data.table")
    
    r["load"](mp_chooser_file)

    with open(file_list_path, "r") as f:
        file_list = [line.strip() for line in f]
    total_files = len(file_list)

    # Store StatPackets temporary.
    tmp_dir = Path("tmp")
    if not tmp_dir.exists():
        tmp_dir.mkdir()

    tmplocalfile = tmp_dir / (file_list_path.name + "_.statpackets")
    
    for i, file in enumerate(file_list):
        print(f"\r{i+1}/{total_files}", end="")
        print(f"\n{i+1}/{total_files} ~> {file}")
        file_path = Path(file)
        if file_path.exists() and ( "NotProcessed" in file or "Successful" in file):
            try:
              df = data_table.fread(
                  file=str(file_path),
                  header=False,
                  sep="\t",
                  quote="",
                  stringsAsFactors=False
              )
              
              # Convert R's ncol and nrow to Python integers
              num_cols = int(r["ncol"](df)[0])
              num_rows = int(r["nrow"](df)[0])
              
              if num_cols != 20 or num_rows > 1:
                  print(f"file ignored (!=20 columns): {file}")
                  continue
                  
              # Call R's annotationChooser
              DRrequiredAgeing = importr("DRrequiredAgeing")
              
              rN = DRrequiredAgeing.annotationChooser(
                  statpacket=df,
                  level=0.0001,
                  rrlevel=0.0001,
                  mp_chooser_file=mp_chooser_file
              )
              
              rW = DRrequiredAgeing.annotationChooser(
                  statpacket=rN.rx2("statpacket"),
                  level=0.0001,
                  rrlevel=0.0001,
                  resultKey="Windowed result",
                  TermKey="WMPTERM",
                  mp_chooser_file=mp_chooser_file
              )
              
              statpacket_v20_values = rW.rx2("statpacket").rx2("V20")
              
              with open(tmplocalfile, "a") as outfile:
                  outfile.write("".join(r["as.character"](statpacket_v20_values)) + "\n")

            except Exception as e:
              print(f"Error processing file {file}: {e}")

if __name__ == "__main__":
    main()
