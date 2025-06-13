"""Run annotationChooser on a list of files."""
import argparse
import logging
from pathlib import Path

from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr


# Activate pandas conversion.
pandas2ri.activate()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_statpacket(file, output_file, mp_chooser_file):
    data_table = importr("data.table")
    df = data_table.fread(
        file=file,
        header=False,
        sep="\t",
        quote="",
        stringsAsFactors=False
    )

    # Convert R's ncol and nrow to Python integers.
    num_cols = int(r["ncol"](df)[0])
    num_rows = int(r["nrow"](df)[0])

    if num_cols != 20 or num_rows > 1:
        logging.error("Invalid file format !=20 columns or >1 rows: %s", file)
        return

    # Call R's annotationChooser.
    dr_required_ageing = importr("DRrequiredAgeing")

    rN = dr_required_ageing.annotationChooser(
        statpacket=df,
        level=0.0001,
        rrlevel=0.0001,
        mp_chooser_file=mp_chooser_file
    )

    rW = dr_required_ageing.annotationChooser(
        statpacket=rN.rx2("statpacket"),
        level=0.0001,
        rrlevel=0.0001,
        resultKey="Windowed result",
        TermKey="WMPTERM",
        mp_chooser_file=mp_chooser_file
    )

    statpacket_v20_values = rW.rx2("statpacket").rx2("V20")

    with open(output_file, "a", encoding="utf-8") as outfile:
        outfile.write(
            "".join(r["as.character"](statpacket_v20_values)) + "\n"
        )

def main():
    """Process statpackets using annotationChooser."""
    parser = argparse.ArgumentParser(
        description="Process files and using annotationChooser."
    )
    parser.add_argument(
        "file", help="Path to the file containing the list of files to process."
    )
    parser.add_argument(
        "mp_chooser_file", help="Path to the mp_chooser.json.Rdata file."
    )
    args = parser.parse_args()

    # Log the arguments
    logging.info("Starting loader.R with arguments:")
    logging.info("  file: %s", args.file)
    logging.info("  mp_chooser_file: %s", args.mp_chooser_file)

    file_list_path = Path(args.file)

    # Load necessary R libraries.
    importr("base")
    importr("jsonlite")
    importr("utils")

    with open(file_list_path, "r", encoding="utf-8") as f:
        file_list = [line.strip() for line in f]
    total_files = len(file_list)

    # Store StatPackets temporary.
    output_dir = Path("annotation_pipeline_output")
    if not output_dir.exists():
        output_dir.mkdir()

    output_file = output_dir / (file_list_path.name + "_.statpackets")

    for i, file in enumerate(file_list):
        file_path = Path(file)
        if file_path.exists() and (
            "NotProcessed" in file or "Successful" in file
        ):
            try:
                process_statpacket(file, output_file, args.mp_chooser_file)
                logging.info("File processed successfully %s out of %s: %s", i+1, total_files, file)
            except Exception as e:
                logging.error(f"Error processing {file}:\n%s", e)

if __name__ == "__main__":
    main()
