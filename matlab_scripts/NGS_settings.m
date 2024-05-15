function Params = NGS_settings()
    Params = struct;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Change pipeline settings below:
    % For directories, 
    %   be sure to add the "/" at the end of the directory name
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pipeline Settings
    Params.read_filetype = ".fastqsanger.gz"; %Tells matlab which file to search for in reads folder. Either ".fastqsanger" or ".fastqsanger.gz"

    
    Params.cpu_cores = true; % true or number
    Params.start_step = 0;
    Params.parpool_timeout = 1440;
    Params.throw_exception = true;
    
    % Set to true if you want the workspace mat file to be deleted after
    % execution of the pipeline to save space. set to false to save .mat
    % files.
    Params.post_pipeline_cleanup = true;

    %% Library Settings
    % Must convert library to .mat table. This speeds up the algorithm
    % If the guide_table_file has not been made yet use"
    %   guide_rna_csv2table("<guide_rna.csv>")
    Params.guide_type_dir = "Cas9/";
    Params.guide_table_file = "Cas9_opt.mat";
    
    %% Background Filepaths (shouldn't need to change)
    Params.galaxy_dir = "./galaxy_dataset/";
    Params.mat_workspace_dir = "./mat_workspace/";
    Params.guide_lib_dir = "./guide_library/";
    Params.results_dir = "./results/";
    
    Params.bam_dir = "bowtie2_bam/";
    Params.untrimmed_dir = "cutadapt/";
    Params.trimmed_dir = "trimmomatic/";
    

end

