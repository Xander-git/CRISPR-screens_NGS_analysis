#!/bin/bash
. import_galaxy_dataset.sh

galaxy_id="<galaxy_api_key"
dataset_name="<dataset_name>"

sample_name="<sample_name>"
readfile_url="<readfile_url>"
bowtie2_url="<bowtie_url>"

import_galaxy "$galaxy_id" "$dataset_name" "$sample_name" "$readfile_url" "$bowtie2_url"