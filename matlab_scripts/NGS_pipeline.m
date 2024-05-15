function [status, msg, err] = NGS_pipeline(collection, sample, read_dir,adapters)
    NGS_SETTINGS = NGS_settings();
    func_name="NGS_pipeline";
    disp("-------------------------------------------------------------------")
    
    warning('off',"MATLAB:MKDIR:DirectoryExists")
    mkdir("../log");
    log = fopen("../log/pipeline_log.txt","w");
    

    msg = sprintf(">> [%s] STARTING EXECUTION(%s)...\n",...
        datetime('now',Format='default'),func_name);
    fprintf(log, msg);
    fprintf(msg)

    
    if NGS_SETTINGS.cpu_cores == true
        cpu_cores = gpuDeviceCount("available");
    else
        cpu_cores = NGS_SETTINGS.cpu_cores;
    end    

    starting_step =NGS_SETTINGS.start_step %#ok<NOPRT> 
    parpool_timeout = NGS_SETTINGS.parpool_timeout %#ok<NOPRT> 
    throw_err = NGS_SETTINGS.throw_exception %#ok<NOPRT> 

    read_dir = string(read_dir)
    collection = string(collection) %#ok<NOPRT> 
    sample = string(sample) %#ok<NOPRT> 

    %%    
    try
        sample_start = tic;

        %% import & count
        if starting_step<=1
            step1_start=tic;
            for i=1:length(adapters)
                [status,msg,err] = step1_import_fastq(collection, sample, read_dir, string(adapters(i)));
                fprintf(log,msg);
                if ~status
                    rethrow(err)
                end

            end
            step1_end = duration(0,0,toc(step1_start));
            msg = sprintf("Elapsed Time For step1_import_fastq for each adapter: %s \n",step1_end);
            fprintf(log,msg);
            fprintf(msg)

        end
        
        %% NEM
        if starting_step<=2
            step2_start = tic;

            pipeline_pool = parpool("Processes",cpu_cores,...
                IdleTimeout=parpool_timeout) %#ok<NOPRT> 

            for i = 1:length(adapters)
                [status,msg,err]=step2_nem(collection,sample,adapters(i));
                fprintf(log,msg);
                if ~status
                    delete(pipeline_pool)
                    rethrow(err)
                end
            end
            delete(pipeline_pool)
            step2_end = duration(0,0,toc(step2_start));
            msg = sprintf("Elapsed Time For step2_nem for each adapter: %s \n",step2_end);
            fprintf(log,msg);
            fprintf(msg)
        end
        
        %% bowtie2
        if starting_step<=3
            step3_start=tic;
            for i = 1:length(adapters)
                [status,msg,err]=step3_bowtie2(collection,sample,adapters(i));
                fprintf(log,msg);
                if ~status
                    rethrow(err)
                end
            end 
            step3_end = duration(0,0,toc(step3_start));
            msg = sprintf("Elapsed Time for step3_bowtie2 for each adapter: %s \n",step3_end);
            fprintf(log,msg);
            fprintf(msg)
        end
        

        %% reconcile
        if starting_step<=4
            step4_start=tic;
            for i = 1:length(adapters)
                [status,msg,err]=step4_reconcile(collection,sample,adapters(i));
                fprintf(log,msg);
                if ~status
                    rethrow(err)
                end
            end
            step4_end = duration(0,0,toc(step4_start));
            msg = sprintf("Elapsed Time for step4_reconcile for each adapter: %s \n",step4_end);
            fprintf(log,msg);
            fprintf(msg)
        end
        
        %% Combine Data
        [status,msg,err]=step5_combine_adapters(collection,sample,adapters);
        fprintf(log,msg);
        if ~status
            rethrow(err)
        end
        
        status=true;
        msg = sprintf(">> [%s] ...FINISHED EXECUTION(%s)\n",datetime('now',Format='default'),func_name);
        fprintf(log, msg);
        fprintf(msg)
        err="";
        
        sample_end = duration(0,0,toc(sample_start));
        msg = sprintf("Total Elapsed Time for sample %s: %s \n",sample, sample_end);
        fprintf(log,msg);
        fprintf(msg)
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",datetime('now',Format='default'),func_name);
        fprintf(log,msg);
        fprintf(msg)
        fprintf(log,err.getReport('extended','hyperlinks','off'));
        if throw_err
            rethrow(err)
        end
    end
    fclose(log);
end