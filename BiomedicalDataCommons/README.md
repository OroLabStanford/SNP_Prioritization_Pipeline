Scripts that can be used to perform analyses drawing data from Google Biomedical Data Commons. In order to use these scripts a Google Data Commons key needs to be generated. This can be done by following documentation on the [Data Commons API](https://docs.datacommons.org/api/). This key should then be set equal to the global variable API_KEY in each script.

###Update 1/12/22###
Biomedical Data Commons supports use of the python API without use of an API key. See [get_genetic_variant_clinical_significance.py](https://github.com/OroLabStanford/SNP_Prioritization_Pipeline/edit/main/BiomedicalDataCommons/get_genetic_variant_clinical_significance.py) for an example of how to use the datacommons python library to query Biomedical Data Commons.

###Update 1/18/22###
Added [get_genetic_variant_functional_category.py](https://github.com/OroLabStanford/SNP_Prioritization_Pipeline/tree/main/BiomedicalDataCommons) which uses the new functionality of Data Commons python API, which doesn't require use of an API key.
