%% Switch Bowtie Counts with NEM counts
function [status,msg,err] = step4_reconcile(collection,sample_name,adapter)
    NGS_SETTINGS = NGS_settings();
    func_name="step4_reconcile";
    clearvars -except func_name date sample adapter

    try
        %%
        disp("-------------------------------------------------------------------")
        fprintf(">> [%s] STARTING EXECUTION(%s)...\n", ...
            datetime('now',Format='default'),func_name)
        collection = string(collection);
        sample_name = string(sample_name);
        adapter_name = string(adapter);
        
        fpath_mat_data = sprintf("%s%s/%s/%s_%s.mat",...
            NGS_SETTINGS.mat_workspace, collection, sample_name, sample_name, adapter_name);
        load(fpath_mat_data,"GUIDE_RNA_SEQUENCE","NEM","BOWTIE")
        whos GUIDE_RNA_SEQUENCE NEM BOWTIE
        
        disp(">> Starting reconcile...")
        tic
        FINAL_COUNT=max(NEM,BOWTIE);
        disp(">> ...Finished reconcile")

        %%
        disp(">> Creating Final Count Table...")
        ADAPTER_TABLE = table(GUIDE_RNA_SEQUENCE,NEM,BOWTIE,FINAL_COUNT);
        head(ADAPTER_TABLE)
        disp(">> ...Finished Table Creation")

        disp(">> Saving Final Count & Results Table...")
        save(fpath_mat_data,"FINAL_COUNT","ADAPTER_TABLE","-append");
        disp(">> ...Finished Saving Results")
        
        status=true;
        msg = sprintf(">> [%s] ...FINISHED EXECUTION(%s)\n",...
            datetime('now',Format='default'),func_name);
        fprintf(msg)
        err="";
        toc
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",...
            datetime('now',Format='default'),func_name);
        fprintf(msg)
    end
        
end