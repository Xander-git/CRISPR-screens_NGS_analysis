function [status,msg,err] = step2_nem(collection,sample_name,adapter)
    NGS_SETTINGS = NGS_settings();
    func_name="step2_nem";
    try
        disp("-------------------------------------------------------------------")
        fprintf(">> [%s] STARTING EXECUTION(%s)...\n", ...
            datetime('now',Format='default'),func_name)
        
        sample_name = string(sample_name);
        adapter_name = string(adapter);
        collection = string(collection);
        
        mat_workspace_dir = NGS_SETTINGS.mat_workspace_dir;
        fpath_guide_data = NGS_SETTINGS.guide_lib_dir + NGS_SETTINGS.guide_table_file;

        
        fpath_mat_data = mat_workspace_dir + collection + "/" + sample_name + "/" + strcat(sample_name,"_",adapter_name,".mat");
        
        disp(">> Loading Variables From step1_import_fastq...")
        load(fpath_mat_data,"READ_SEQUENCE","ABUNDANCE")
        disp(">> Loading Guide RNA Library...")
        load(fpath_guide_data,"guide_table")
        
        nem_start = tic;
        
        GUIDE_RNA_SEQUENCE = guide_table.(2);
        lib_sz = length(GUIDE_RNA_SEQUENCE);
        NEM=zeros(lib_sz, 1);
        
        clearvars guide_table
        % start = datetime('now');
        whos GUIDE_RNA_SEQUENCE NEM READ_SEQUENCE ABUNDANCE
        
        disp(strcat(">> Starting NEM with Library_Size:",string(lib_sz)," & Sample_Size:",string(length(READ_SEQUENCE)),"..."))
        parfor i=1:lib_sz
        
            w1 = find( contains(READ_SEQUENCE, GUIDE_RNA_SEQUENCE(i) ) )
            if ~any(w1)
                NEM(i,1)=0;
            else
                NEM(i,1)=sum(ABUNDANCE(w1)) %#ok<PFBNS> 
            end
        end
        disp(strcat(">> ...Finished NEM ",adapter_name))
        

        %%
        disp(">> Previewing NEM Results Table...")
        nem_table = table(GUIDE_RNA_SEQUENCE,NEM);
        head(nem_table)
        
        disp(">> Saving NEM Counts, Guide RNA Sequence, & Results Table...")
        save(fpath_mat_data,"NEM","GUIDE_RNA_SEQUENCE","-append")
        disp(">> ...Finished Saving Results")

        nem_end = duration(0,0,toc(nem_start));
        msg = sprintf(">> Elapsed Time For NEM Counting: %s\n",nem_end);
        fprintf(msg)

        status = true;
        err="";
        msg = sprintf(">> [%s] ...FINISHED EXECUTION(%s)\n",datetime('now',Format='default'),func_name);
        fprintf(msg)
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",datetime('now',Format='default'),func_name);
        fprintf(msg)
    end
end