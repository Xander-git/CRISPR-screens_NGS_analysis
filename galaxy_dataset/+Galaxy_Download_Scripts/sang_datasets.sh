#!/bin/bash

# Inputs 

source ./galaxy_sgRNA_curl_func.sh

category_name="sang"
galaxy_key="dc8571d5abfdeaaf240f10aabf6e6fba"

##################
sample_name="top_lipid_survivor"

cutadapt_dLink="http://localhost:8080/api/histories/6ee95371d25d237b/contents/dataset_collections/634f2991fe1a386d/download"
bowtie2_bam_dLink="http://localhost:8080/api/histories/6ee95371d25d237b/contents/dataset_collections/d343a822bd747ee4/download"

##################

sgRNA_galaxy_curl $category_name $sample_name $galaxy_key $cutadapt_dLink $bowtie2_bam_dLink

##################
sample_name="whole_survivor"

cutadapt_dLink="http://localhost:8080/api/histories/7e2da049d9f5c1fa/contents/dataset_collections/b3a78854daef4a5a/download"
bowtie2_bam_dLink="http://localhost:8080/api/histories/7e2da049d9f5c1fa/contents/dataset_collections/b4a9efd5c50010ba/download"

##################

sgRNA_galaxy_curl $category_name $sample_name $galaxy_key $cutadapt_dLink $bowtie2_bam_dLink
