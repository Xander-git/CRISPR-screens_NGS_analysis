#!/bin/bash

import_galaxy () {
    local dataset_label="$2"
    local sample_name="$3"
    local galaxy_id="$1"

    local readfile_url="$4"
    local bowtie2_bam_url="$5"

    readfile_url+="?key=$galaxy_id"
    bowtie2_bam_url+="?key=$galaxy_id"
    echo "$readfile_url"
    echo "$bowtie2_bam_url"

    mkdir -p "../../results/${dataset_label}/${sample_name}/"
    mkdir -p "../../galaxy_dataset/${dataset_label}/${sample_name}/"

    curl --parallel \
        -o "../../galaxy_dataset/${dataset_label}/${sample_name}/readfile.zip" "${readfile_url}" \
        -o "../../galaxy_dataset/${dataset_label}/${sample_name}/bowtie2_bam.zip" "${bowtie2_bam_url}"

    unzip "../../galaxy_dataset/${dataset_label}/${sample_name}/*" -d "../../galaxy_dataset/${dataset_label}/${sample_name}"   
}