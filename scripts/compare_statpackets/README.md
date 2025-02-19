# Compare Statpackets
`compare.sh` script can be used to normalise and compare JSON files.
`preprocess.py` masks fields that changes between data releases.
It is a prototype of the script, so it is not optimised for large scale usage.

# How to Run
## 1. Install Dependencies
```console
python3 -m venv env
source env/bin/activate
python3 -m pip install -r requirements.txt
``` 
Script is running locally, but not on SLURM, because of delta dependency that can't be installed (requires sudo priveleges).

## 2a. Run for a JSON file
One row of the statpackets is a JSON file.
```console
head -n1 split_index_aa_1.statpackets > 1.json
head -n1 split_index_aa_2.statpackets > 2.json
touch report.html.gz 
bash compare.sh \
  1.json \
  2.json \
  report.html.gz 
```
  
## 2b. Run for Statpackets
One statpackets is a file with 50 rows where each row is a JSON file.
```console
mkdir reports
parallel touch reports/report_{1}.html.gz ::: {1..50}
parallel --colsep ' ' bash compare.sh \
  <(sed -n "{1}p" 21.0.statpackets) \
  <(sed -n "{1}p" 21.3.statpackets) \
  reports/report_{1}.html.gz ::: {1..50}
```
The script will generate 50 reports in the `report` folder. 

## 3. Unzip Report
```console
gzip -dk report.html.gz
```