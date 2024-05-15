#!/bin/bash

sgRNA_galaxy_curl () {
    local category_label="$1"
    local sample_name="$2"
    local galaxy_id="$3"

    local cutadapt_url="$4"
    local bowtie2_bam_url="$5"

    cutadapt_url+="?key=$galaxy_id"
    bowtie2_bam_url+="?key=$galaxy_id"
    echo "$cutadapt_url"
    echo "$bowtie2_bam_url"

    mkdir -p "../../results/${category_label}/${sample_name}/"
    mkdir -p "../${category_label}/${sample_name}/"

    curl --parallel \
        -o "../${category_label}/${sample_name}/trimmomatic.zip" "${cutadapt_url}" \
        -o "../${category_label}/${sample_name}/bowtie2_bam.zip" "${bowtie2_bam_url}"

    unzip "../${category_label}/${sample_name}/*" -d "../${category_label}/${sample_name}"
    gunzip -r "../${category_label}/${sample_name}/"

    
}