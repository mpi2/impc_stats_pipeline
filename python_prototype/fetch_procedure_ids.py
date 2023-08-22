import json
from urllib.request import urlopen

url = "https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*:*&rows=0&facet=on&facet.field=procedure_stable_id&facet.limit=-1"
response = urlopen(url)
data = json.load(response)
procedures_and_rows = data["facet_counts"]["facet_fields"]["procedure_stable_id"]
procedures = procedures_and_rows[::2]
row_counts = procedures_and_rows[1::2]
for procedure, row_count in zip(procedures, row_counts):
    fetch_url = f"https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=*:*&fq=procedure_stable_id:{procedure}&rows={row_count}"
    print(f"{procedure}\t{fetch_url}")
