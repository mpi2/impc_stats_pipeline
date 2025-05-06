# Python prototype for quick preparation of chunks

## 1. Install dependencies
To execute commands you will need: curl, jq, parallel, python3 

## 2. Get list of all procedure_stable_id chunks
```bash
python3 fetch_procedure_ids.py > procedure_chunks.txt
```

## 3. Download the input data
```bash
mkdir input_procedures
parallel --progress --jobs 4 --colsep "\t" \
  curl --silent --output input_procedures/{1}.csv {2} \
  :::: procedure_chunks.txt
```

## 4. Chunk the input data
```bash
mkdir intermediate results
parallel --progress --colsep "\t" \
  python3 chunking.py input_procedures/{1}.csv intermediate/ results/ \
  :::: procedure_chunks.txt
```
