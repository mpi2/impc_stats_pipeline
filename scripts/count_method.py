#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""This script processes a ZIP file to count occurrences
of applied methods within output_Successful.tsv files in Results_IMPC_SP_Windowed.zip.
The script is executed from the command line with the path to the ZIP file as an argument. 
The results are outputted in separate JSON files for each pattern, summarizing the counts
of the applied methods found in the processed files."""

import json
import re
import sys
import zipfile
from collections import Counter

def count_method(input_zip):
    # Initialize counters
    counters = {p: Counter() for p in ("_IPG_010_001", "_IPG_011_001", "_IPG_012_001")}
    
    # Open the ZIP file
    with zipfile.ZipFile(input_zip, 'r') as zip_file:
        # Iterate through each file in the ZIP
        for file_info in zip_file.infolist():
            # Get the file name
            file_name = file_info.filename
            # Check all patterns
            for pattern, counter in counters.items():
                # Check if the pattern matches the file name and if the file name is "output_Successful.tsv"
                if re.search(pattern, file_name) and file_name.endswith("output_Successful.tsv"):
                    # Extract file content
                    with zip_file.open(file_name) as file:
                        file_content = file.read().decode('utf-8')
                        json_string = file_content.split('\t')[-1]
                        # Convert JSON string to Python dictionary
                        stats_results = json.loads(json_string)
                        # Extract the value of "Applied method"
                        applied_method = stats_results["Result"]["Vector output"]["Normal result"]["Applied method"]
                        # Update the counter
                        counter.update([applied_method])
    
    # Write counters to separate files
    for param, counter in counters.items():
        with open("".join([param[-11:], ".json"]), 'w') as f:
            json.dump(counter, f, indent=4)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process zip file and count applied methods.")
    parser.add_argument("zip_file", help="Path to the zip file")
    args = parser.parse_args()
    count_method(args.zip_file)

