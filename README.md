# CRISPR-screens_NGS_analysis
---
## Quick Sample Analysis

1. Add scripts to matlab path
> addpath("./src")
2. Run  pipeline on collection
---
## Syntax
### Function Parameter Definitions
`<collection_name>`: String with the name of the collection to search for the sample in

`<sample_name>`: String with the name of the sample to search for

`<read_directory>`: String with the name of the directory to use name exact matching on. By convention it should be "cutadapt" for untrimmed reads or "trimmomatic" for trimed reads.

`{<adapter1>,... <adapterN>}`: A cell array with the name of each adapter of interest in the form of a string or character array. 

`<adapter>`: A string with the name of a single adapter

### Pipeline Run
> `NGS_pipeline(<collection_name>, <sample_name>, <read_directory>, {<adapter1>,... <adapterN>});`
Executes the code for step 1-5 on every adapter for a single sample

### Batch Run
> `NGS_batch_pipeline(<collection_name>, <read_directory>, {<adapter1>,... <adapterN>});`

### Step 1: Import fastqsanger
> `step1_import_fastq(<collection_name>, <sample_name>, <read_directory>, {<adapter1>,... <adapterN>});` 
Imports and counts the reads inside a .fastqsanger or .fastqsanger.gz file. The unique sequences and abundance are then stored in a .mat file to motivate matlab import and export. Should be capable of handling a large range of file sizes. Recommended to use .fastqsanger.gz option to save on memory space.

### Step 2: Name Exact Matching
> `step2_nem(<collection_name>, <sample_name>, {<adapter1>,... <adapterN>});`
Imports the data from the .mat variable file created in step 1 and matches it to the guide rna library specified in the settings.

### Step 3: Process Bowtie2
> `step3_bowtie2(<collection_name>, <sample_name>, <adapter>)`
Processes and generates the counts given by the bowtie2 algorithm from galaxy.

### Step 4: Reconcile Counts
> `step4_reconcile(<collection_name>, <sample_name>, <adapter>)`
Compares the NEM and bowtie2 counts of a specific adapter for a given sample and reports the larger of the two as the final count.

### Step 5: Combine Adapters for each Sample
> `step5_combine_adapters(<collection_name>,<sample_name>, {<adapter1>,... <adapterN>})`
Combines the final count for each adapter onto one table for ease of analysis.
---
## NGS_Settings.m

The settings for the NGS_pipeline are stored here. These are vital to the underlying code used. This works similar to class properties.

 > `Params.guide_table_file`
Default: "Cas9_Opt.mat"


> `Params.read_filetype`
Default: ".fastqsanger.gz" 

Must either be ".fastqsanger.gz" or ".fastqsanger". The pipeline is capable of reading unzipped or zipped fastqsanger files. By default this is set to the zipped files in an effort to save storage space at the expense of a slight runtime increase

 > `Params.cpu_cores`
 Default: true

 Can either be true or a specific number. If true, matlab will check how many available cores are available for use and use the maximum amount during parallel processing. Otherwise a specific number will tell the parallel pool to only use that amount. This is useful if sharing cpu cores with other users.
 
 > `Params.start_step`
 Default: 0
 
 Tells the pipeline which step to start at. Useful in cases where the pipeline crashes at a certain step.

 > `Params.parpool_timeout`
Default: 1440

Determines how many minutes the parallel process will run before timing out. If NEM continously times out, it usually means this value is too low. In most cases 1,440 minutes is long enough for NEM to finish. 

 > `Params.throw_exception`
Default: true

Tells the pipeline whether to throw exceptions or continue

 > `Params.post_pipeline_cleanup`
 Default: true

 once the pipeline is done, it removes any unzipped files and deletes the generated .mat files in order to save space


## Filepath Settings
In general, these should not be changed unless you want to change the directory structure. Essential for pipeline algorithms to function properly. All filepaths should end with a path separator.
 > `Params.galaxy_dir`
 Default:"./galaxy_dataset/"

 Filepath to where your galaxy data is stored.

 > `Params.mat_workspace_dir`
 Default:"./src/mat_workspace/"

 > `Params.guide_lib_dir`
 Default:"./guide_library/"

 > `Params.results_dir`
 Default: "./results/"

 > `Params.bam_dir`
 Default: "./bowtie2_bam/"
