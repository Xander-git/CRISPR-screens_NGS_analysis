#!/bin/bash

# Inputs 

source ./~galaxy_sgRNA_curl_func.sh

##################

category_name="<category_name>"
galaxy_key="<galaxy_id>"

##################
sample_name="<sample_name>"

reads_url="<reads_collection_link>"
bowtie2_bam_url="<bowtie2_bam_collection_link>"

##################

sgRNA_galaxy_curl $category_name $sample_name $galaxy_key $reads_url $bowtie2_bam_url