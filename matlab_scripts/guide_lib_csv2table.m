function [status, msg, err] = guide_lib_csv2table(guide_rna)
    SETTINGS = settings();
    func_name="guide_rna_csv2table";
    status=false;err="";msg=""; %#ok<NASGU> 
    try
        status=true;
        err="";

        tic
        fprintf(">> Beginning Import of '"+ guide_lib_path + "'...\n")

        % guide_lib = readtable(guide_lib_path,"VariableNamingRule","preserve");
        guide_table = readtable(guide_lib_path+".csv");
        save(strcat(SETTINGS.guide_lib_folder,guide_rna,".mat"),"guide_table")
        toc
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n",datetime('now',Format='default'),func_name);
    end
end