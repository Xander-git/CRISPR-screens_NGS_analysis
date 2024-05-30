function [status,msg,err] = step5_combine_adapters(collection, sample_name, adapter_names)
% STEP5_COMBINE_ADAPTERS Combines the final counts for each adapter and
% outputs a csv file
%   Detailed explanation goes here
    NGS_SETTINGS = NGS_settings();
    func_name="step5_combine_adapters";

    try
        disp("-------------------------------------------------------------------")
        fprintf(">> [%s] STARTING EXECUTION(%s)...\n",datetime('now',Format='default'),func_name)
        collection = string(collection);
        sample_name = string(sample_name);
        
        fpath_results_collection_dir = strcat(NGS_SETTINGS.results_dir, collection,"/");
        fpath_results_sample_dir = strcat(fpath_results_collection_dir,sample_name,"/");

        fpath_mat_sample_prefix = strcat(NGS_SETTINGS.mat_workspace_dir,collection,"/",sample_name,"/",sample_name,"_");
        fpath_guide_table = NGS_SETTINGS.guide_lib_dir + NGS_SETTINGS.guide_table_file;

        warning('off',"MATLAB:MKDIR:DirectoryExists")
        mkdir(fpath_results_collection_dir);
        mkdir(fpath_results_sample_dir);
        
        load(fpath_guide_table,"guide_table");
        guide_rna_id = guide_table.(1);
        sequence = guide_table.(2);
        cutting_score = guide_table.(3);
        guide_position = guide_table.(4);
        library_design_information = guide_table.(5);

        full_table_names = {'Guide RNA ID', 'Sequence', 'Cutting Score', 'Guide Position', 'Library Design Information'};
        full_count = table(guide_rna_id, sequence,cutting_score, guide_position, library_design_information, 'VariableNames', full_table_names);
        
        for i = 1:length(adapter_names)
            adapter = load(strcat(fpath_mat_sample_prefix,adapter_names(i),".mat"),"FINAL_COUNT","NEM","BOWTIE");
            count = table(adapter.FINAL_COUNT,'VariableNames',[adapter_names(i)]);

            disp(">> Creating Adapter "+adapter_names(i)+" Report Table...")
            adapter_table=table(sequence, adapter.NEM, adapter.BOWTIE, adapter.FINAL_COUNT);
            
            writetable(adapter_table,strcat(fpath_results_sample_dir,"/",sample_name,"_",adapter_names(i),".csv"));

            disp(">> Creating Full Report Table...")
            full_count = [full_count count]; %#ok<AGROW> 
        end
        
        head(full_count)
        
        disp(">> Saving Results")
        save(strcat(fpath_mat_sample_prefix,"full.mat"), "full_count",'-v7.3')
        writetable(full_count,strcat(fpath_results_sample_dir,"/",sample_name,"_full_count.csv"))
        
        if NGS_SETTINGS.post_pipeline_cleanup
            [status, msg] = rmdir((NGS_SETTINGS.mat_workspace_dir + collection + "/" + sample_name),'s') %#ok<ASGLU> 
        end

        status=true;
        msg = sprintf(">> [%s] ...FINISHED EXECUTION(%s)\n",datetime('now',Format='default'),func_name);
        fprintf(msg)
        err="";
        toc
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",datetime('now',Format='default'),func_name);
        fprintf(msg)
    end
end

