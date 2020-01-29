Exception URL for procedure_stable_ids (not used): https://www.ebi.ac.uk/mi/impc/solr/pipeline/select?q=annotate:false&facet=true&facet.field=procedure_stable_id&rows=0&facet.mincount=1&facet.limit=-1


Exception URL for parameter_stable_ids (USED): https://www.ebi.ac.uk/mi/impc/solr/pipeline/select?q=annotate:false&facet=true&facet.field=parameter_stable_id&rows=0&facet.mincount=1&facet.limit=-1


## already solved. For the category list: write.csv(DRrequiredAgeing::GetPossibleCategories(file = FALSE),file = 'd:AllCts.csv')

For ABR parameters (method list ==RR): http://ves-ebi-d0:8090/mi/impc/dev/solr/experiment/select?q=*%3A*&fq=procedure_group%3A*ABR*&rows=0&fl=parameter_stable_id&wt=json&indent=true&facet=true&facet.field=parameter_stable_id&facet.limit=10000&facet.mincount=1