import argparse
import os
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

# Activate pandas conversion
pandas2ri.activate()

def main():
    parser = argparse.ArgumentParser(description="Process files and annotate them using R's annotationChooser.")
    parser.add_argument("file", help="Path to the file containing the list of files to process.")
    parser.add_argument("mp_chooser_file", help="Path to the mp_chooser.json.Rdata file.")
    args = parser.parse_args()

    file_list_path = args.file
    mp_chooser_file = args.mp_chooser_file

    # Load necessary R libraries
    base = importr('base')
    jsonlite = importr('jsonlite')
    utils = importr('utils')
    data_table = importr('data.table')
    
    r['load'](mp_chooser_file)
    a = r['a']

    # Set levels
    level = 0.0001
    rrlevel = 0.0001

    # Start annotation pipeline
        # today = r['format'](r['Sys.time'](), "%d%m%Y")
        # today_str = r['as.character'](today)[0] # Convert R character vector to python string.

    with open(file_list_path, 'r') as f:
        file_list = [line.strip() for line in f]
    lflist = len(file_list)

    # Store StatPackets temporary
    if not os.path.exists("tmp"):
        os.makedirs("tmp")

    tmplocalfile = os.path.join('tmp', os.path.basename(file_list_path) + '_.statpackets')
    
    for i, file in enumerate(file_list):
        print(f"\r{i+1}/{lflist}", end="")
        print(f"\n{i+1}/{lflist} ~> {file}")
        if os.path.exists(file) and ( 'NotProcessed' in file or 'Successful' in file):
            try:
              df = data_table.fread(
                  file=file,
                  header=False,
                  sep='\t',
                  quote="",
                  stringsAsFactors=False
              )
              
            #   if r['ncol'](df) != 20 or r['nrow'](df) > 1:
            #       print(f'file ignored (!=20 columns): {file}')
            #       continue
                  
              # Call R's annotationChooser
              DRrequiredAgeing = importr('DRrequiredAgeing')
              
              rN = DRrequiredAgeing.annotationChooser(
                  statpacket=df,
                  level=level,
                  rrlevel=rrlevel,
                  mp_chooser_file=mp_chooser_file
              )
              
              rW = DRrequiredAgeing.annotationChooser(
                  statpacket=rN.rx2('statpacket'),
                  level=level,
                  rrlevel=rrlevel,
                  resultKey='Windowed result',
                  TermKey='WMPTERM',
                  mp_chooser_file=mp_chooser_file
              )
              
              statpacket_v20_values = rW.rx2('statpacket').rx2('V20')
              
              with open(tmplocalfile, 'a') as outfile:
                  outfile.write("".join(r['as.character'](statpacket_v20_values)) + "\n")

            except Exception as e:
              print(f"Error processing file {file}: {e}")

if __name__ == "__main__":
    main()
