# CRISPR-screens_NGS_analysis
## Installation
1. In target directory, clone the repository
```
$ git clone https://github.com/Xander-git/CRISPR-screens_NGS_analysis.git
```

## Quick Sample Analysis using an example dataset

1. From inside the repo's main directory, run matlab from the terminal
> `CRISPR-screens_NGS_analysis$ matlab`

2. Add scripts to matlab path
> `addpath("./src");`

3. (a) Run pipeline on a single sample
> `NGS_pipeline("example_dataset","sample_name", "trimmomatic",{"1665","1666"});`

3. (b) Run pipeline on every sample in a dataset directory using the trimmed reads
> `NGS_batch_pipeline("example_dataset", "trimmomatic",{"1665","1666"});`
---
## Pipeline Syntax

### Pipeline Run
```
NGS_pipeline(<dataset_name>, <sample_name>, <read_directory>, {<adapter1>,... <adapterN>});
```
> Use Cases: When each sample within a dataset has different adapters to target, or if batch pipeline gets cut off 

Executes the code for step 1-5 on every adapter for a single sample

### Batch Run
```
NGS_batch_pipeline(<dataset_name>, <read_directory>, {<adapter1>,... <adapterN>});`
```
> Use Cases: When an entire dataset uses the same adapters

Executes the pipeline run code for every sample within a dataset directory. In the galaxy dataset folder, you should organize your files with a dataset folder containing each individual sample you plan to analyze.The order of folders should look like this.
> galaxy_dataset -> dataset -> samples -> reads_folder & bowtie2_bam -> individual adapter files
---
### Function Parameter Definitions
`<dataset_name>`: String with the name of the target dataset directory within the galaxy_dataset directory

`<sample_name>`: String with the name of the target sample to search for

`<read_directory>`: String with the name of the directory to use naive exact matching on. By convention it should be `"cutadapt"` for untrimmed reads or `"trimmomatic"` for trimmed reads.

`{<adapter1>,... <adapterN>}`: A cell array with the name of each adapter of interest in the form of a string or character array. 

`<adapter>`: A string with the name of a single adapter

---
## Individual Function Calls
> It is recommended that you use the pipeline functions above and just change the starting step in the settings. Using the individual functions is more complex due to the lazy-import implementation. In general, step 1 must always be run first before any others.

### Step 1: Import fastqsanger
```
step1_import_fastq(<dataset_name>, <sample_name>, <read_directory>, {<adapter1>,... <adapterN>});
```
Imports and counts the reads inside a .fastqsanger or .fastqsanger.gz file. The unique sequences and abundance are then stored in a .mat file to motivate matlab import and export. Should be capable of handling a large range of file sizes. Recommended to use .fastqsanger.gz option to save on memory space.

### Step 2: Name Exact Matching
```
step2_nem(<dataset_name>, <sample_name>, {<adapter1>,... <adapterN>});
```

Imports the data from the .mat variable file created in step 1 and matches it to the guide rna library specified in the settings.

### Step 3: Process Bowtie2
```
step3_bowtie2(<dataset_name>, <sample_name>, <adapter>);
```
Processes and generates the counts given by the bowtie2 algorithm from galaxy.

### Step 4: Reconcile Counts
```
step4_reconcile(<dataset_name>, <sample_name>, <adapter>);
```
Compares the NEM and bowtie2 counts of a specific adapter for a given sample and reports the larger of the two as the final count.

### Step 5: Combine Adapters for each Sample
```
step5_combine_adapters(<dataset_name>, <sample_name>, {<adapter1>,... <adapterN>})
```
Combines the final count for each adapter onto one table for ease of analysis.

---
## Helper Functions 
### Convert guide library .csv to a .mat table
```
guidelib_csv2table("<name_of_csv_file>");
```
Converts a csv file with the guide library into a mat file for the pipeline to read. The column order should be the Guide RNA ID, Guide RNA Sequence, Cutting Score, Guide Position, and Library Design Information. See the included CSV file in guide_library for reference. 

### Get bowtie2 start positions
```
get_bowtieStartPos("<fpath_to_bamfile>")
```
A helper function to get the starting positions from bowtie2 given a bam file. The function returns a matlab table with the bowtie2 reference id, the start position, and the corresponding sequence.

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
