    %% Input Settings
    function [status, msg, err] = step1_import_fastq(collection, read_dir, sample_name, adapter)
    NGS_SETTINGS = NGS_settings();
    func_name="step1_import_fastq";
    try
        disp("-------------------------------------------------------------------")
        fprintf(">> [%s] STARTING EXECUTION(%s)...\n",...
            datetime('now',Format='default'),func_name)
        
        collection = string(collection);
        sample_name = string(sample_name);
        read_dir = string(read_dir);
        
        galaxy_dir = NGS_SETTINGS.galaxy_dir;
        mat_workspace = NGS_SETTINGS.mat_workspace;

        galaxy_collection_dir = galaxy_dir + collection + "/";
        
        read_fType = NGS_SETTINGS.read_filetype;
        collection_mat_dir = mat_workspace + collection + "/";
        sample_mat_dir = collection_mat_dir + "/" + sample_name %#ok<NOPRT> 
        
        warning('off',"MATLAB:MKDIR:DirectoryExists")
        mkdir(mat_workspace);
        mkdir(collection_mat_dir);
        mkdir(sample_mat_dir);
        
        



        fpath_read_dataset = strcat(galaxy_collection_dir, sample_name, "/",read_dir,"/",adapter, read_fType) %#ok<NOPRT> 
        fpath_output = strcat(sample_mat_dir,"/",sample_name,"_",adapter, ".mat") %#ok<NOPRT> 
        
        %% Check for gzip
        if read_fType==".fastqsanger.gz"
            fprintf(strcat(">> Beginning decompression of ","'",fpath_read_dataset,"' with adapter ",adapter,"...\n"))
            gunzip(fpath_read_dataset);
            fpath_read_dataset = galaxy_collection_dir + sample_name + "/" + read_dir + "/" + adapter + ".fastqsanger";
        end



        %% Sequence Import Start
        
        fprintf(strcat(">> Beginning Import of ","'",fpath_read_dataset,"' with adapter ",adapter,"...\n"))
        tic
        for a=1:1
            [~,sequence] =  fastqread(fpath_read_dataset);
            if ischar(sequence)==1
                seq=sequence;
                sequence=cell(1);
                sequence{1}=seq;
            end
        end
        Seq=sequence';
        [SAMPLE_SEQUENCE, ABUNDANCE]=count_unique(Seq);
        disp(">> ...Finished Import & Abundance Counting")

        sample_table = table(SAMPLE_SEQUENCE, ABUNDANCE);
        head(sample_table)
        whos S1 N1

        disp(">> Saving Sample Sequences & Abundance Counts...")
        save(fpath_output,"SAMPLE_SEQUENCE","ABUNDANCE",'-v7.3') % S: sequence, N: Unique Counts
        
        
        if read_fType==".fastqsanger.gz"
            fprintf(strcat(">> Beginning compression of ","'",fpath_read_dataset,"' with adapter ",adapter,"...\n"))
            gzip(fpath_read_dataset);
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
