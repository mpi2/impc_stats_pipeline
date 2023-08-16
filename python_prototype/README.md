# Python prototype for quick preparation of chunks

## 1. Install dependencies
To execute commands you will need: curl, jq, parallel, python3 

## 2. Download the input data
Run the following commands to download all input data:
```bash
num_rows=$(curl --silent 'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*:*&rows=0' | jq '.response.numFound')
curl "https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*:*&wt=csv&rows=${num_rows}" | gzip -c > data.csv.gz
```

## 3. Sorting the file into preliminary chunks 
```bash
python3 step_1_preliminary_chunking.py data.csv.gz results/
```

## 4. Splitting preliminary chunks into final chunks 
```bash
ls results/ | sed -e "s/_control.csv//g" -e "s/_experimental.csv//g" | sort -u > lists_of_chunks.txt
parallel python3 step_2_final_chunking.py {}_control.csv {}_experimental.csv :::: list_of_chunks.txt
```
