function [status, msg, err] = guide_lib_csv2table(guide_rna)
    NGS_SETTINGS = NGS_settings();
    func_name="guide_rna_csv2table";
    try
        status=true;
        err="";

        tic
        fprintf(">> Beginning Import of '"+ NGS_SETTINGS.guide_lib_dir + guide_rna + "'...\n")

        % guide_lib = readtable(guide_lib_path,"VariableNamingRule","preserve");
        guide_table = readtable(NGS_SETTINGS.guide_lib_dir + guide_rna + ".csv");
        head(guide_table)
        save(strcat(NGS_SETTINGS.guide_lib_dir,guide_rna,".mat"),"guide_table")
        toc
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",datetime('now',Format='default'),func_name);
    end
end