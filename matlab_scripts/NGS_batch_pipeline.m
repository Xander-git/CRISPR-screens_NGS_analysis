function [] = NGS_batch_pipeline(collection, read_dir, adapters)
disp("-------------------------------------------------------------------")
fprintf(">> [%s] STARTING EXECUTION(sgRNA_batch_pipelineX)...", datetime('now',Format='default'))
NGS_SETTINGS = NGS_settings();

warning('off',"MATLAB:MKDIR:DirectoryExists")
mkdir("../log");
log = fopen('../log/batch_pipeline_log.txt','w');

collection = string(collection) %#ok<NOPRT> 

folder_info = dir(strcat(NGS_SETTINGS.galaxy_dir,collection,"/"));
sample_idx = find(vertcat(folder_info.isdir));
samples = folder_info(sample_idx); %#ok<FNDSB> 
samples = {samples.name};
samples = samples(3:length(samples)) %#ok<NOPRT> 


batch_start = tic;
for i = 1:length(samples)
    disp("-------------------------------------------------------------------")
    msg = sprintf(">> [%s] Starting sgRNA pipeline on sample %s...\n", ...
        datetime('now',Format='default'), string(samples(i)));
    fprintf(log,msg);
    fprintf(msg)

    [status,msg,err] = NGS_pipeline( collection, samples(i), read_dir, adapters);

    fprintf(log,msg);
    if ~status
        fprintf(log,err.getReport('extended','hyperlinks','off'));
        fclose(log);
        rethrow(err)
    end

end

batch_end = duration(0,0,toc(batch_start));
fprintf(">> [%s] FINISHED EXECUTION(sgRNA_pipeline)\n",datetime('now',Format='default'))
fprintf("Total Elapsed Time for Batch Analyis: %s \n",batch_end)
fclose(log);
end