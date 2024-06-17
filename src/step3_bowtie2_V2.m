% This script generates sgRNA counts for a given sample by inexact
% matching. This is done by generating genomic alignment start positions
% for all reads in a given sample. Then these start positions are compared
% to the start positions that we determined for our library of sgRNA
% previously, and if they match, then the read is mapped to the corresponding
% sgRNA.

% While Bowtie itself allows for mismatches while aligning reads to the
% genome, We can also have reads are aligned 1 nt to the left or right of
% the true start position due to sequencing errors. To capture this, we
% allow for some wiggle room in the start position mapping. If the start
% position of a read matches to the start position of an sgRNA to within 1
% nt in either direction, then we map that read to that sgRNA.

% Inputs:
% 1. Examplefile1_Cas9_bam. Obtained as a result of running bowtie on
% 'Examplefile1_Cas9.fastqsanger' against the genome file
% 'CLIB89+control_Cas9.fasta'.
% 2. All_sgRNA_pos_Cas9.mat: variable file containg genomic start positions of
% all Cas9 sgRNA in the library.

% Outputs:
% Bowtie_counts_1.mat: Variable file containing Bowtie counts for a given
% sample.

% *Note*: This script also requires the use of the fast and efficient custom
% 'count_unique' function written by the user Anthony Kendall. This
% function may be downlaoded at:
% https://www.mathworks.com/matlabcentral/fileexchange/23333-determine-and-count-unique-values-of-an-array


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Reading BAM file and storing alignments against all chromosomes

% Inputs are Bam files that are the result of aligning sgRNA from any
% particular sample against the genome using Bowtie. There will again be an
% extra chromosome here as all the nontargeting sgRNA were added as a
% seperate chromosome to allow for the aligment of the entire library.

function [status,msg,err]=step3_bowtie2(collection, sample_name, adapter)
    NGS_SETTINGS = NGS_settings();
    func_name="step3_bowtie2";
    
    try
        disp("-------------------------------------------------------------------")
        fprintf(">> [%s] STARTING EXECUTION(%s)...\n",datetime('now',Format='default'),func_name)
        
        collection = string(collection);
        galaxy_dataset_dir = strcat(NGS_SETTINGS.galaxy_dir, collection, "/");

        sample_name = string(sample_name);
        adapter_name = string(adapter);
        
        fpath_guide_library = NGS_SETTINGS.guide_lib_dir + NGS_SETTINGS.guide_table_file;
        fpath_bam_dataset = strcat(galaxy_dataset_dir, sample_name, "/", NGS_SETTINGS.bam_dir, adapter_name, ".bam"); 
        fpath_mat_data = strcat(NGS_SETTINGS.mat_workspace_dir,collection,"/",sample_name,"/",sample_name,"_",adapter_name,".mat");
        
        %%
        fprintf(strcat(">> Starting bowtie2 counting on sample ",sample_name," with adapter ", adapter_name,"\n"))
        
        bowtie_start = tic;
        bamFilename = fpath_bam_dataset; %This file will be the result of running bowtie on 'Examplefile1_Cas9.fastqsanger' against the genome file 'CLIB89+control_Cas9.fasta'
        info = baminfo(bamFilename,'ScanDictionary',true);
        
        bm_cArr = {};
        for i = 1:length(info.ScannedDictionary)
            try
                bm_tmp = BioMap(bamFilename, 'SelectReference',info.ScannedDictionary{i});
                fprintf(">> Found Bowtie Reference: " + info.ScannedDictionary{i} + "\n");
                bm_cArr = [bm_cArr,{bm_tmp}];
            catch
                fprintf(">> Could not read reference: "+info.ScannedDictionary{i}+" (If reference = unmapped, ignore this warning)\n");
            end
        end
        fprintf(">> Amount of references found: " + length(bm_cArr)+"\n");

        load(fpath_guide_library,"guide_table");
        guide_lib = table2cell(guide_table);

        %% Concatenate Start Positions and Flags
        
        % Start positions of every sgRNA in the current reads file is recorded.
        % Other information that is appended to the start position of an sgRNA is
        % the chromosome it aligned to, and the strand that it was found on. For
        % example: a start position of 2500018 is read this way: 25000|1|8. 8
        % refers to the chromosome, a flag of 0 or 1 mean Top and Bottom strands
        % respectively, and the remaining is the actual start position of the
        % sgRNA. So this sgRNA is found at the 25000 position on the bottom strand
        % in the 8th chromosome. Actions here are repeated for all 6 chromosomes
        % in Y. lipolytica (A,B,C,D,E,F).
        
        chromosomes_cArr={};
        for i = 1:length(bm_cArr)
            msg = ">> Starting alignment sequence import for chromosome number: " + i + "\n";
            fprintf(msg);
            curr_bm = bm_cArr{i};
            start = curr_bm.Start;
            flag = curr_bm.Flag;
            flag(flag==16)=1; % The reverse strand alignment to a chromosome is indicated as 16 in bowtie, is changed to 1 here
            
            fn=zeros(length(start),1);
            for j = 1:length(start)
                as = sprintf('%1d%1d%1d ', start(j), flag(j), i);
                anum=str2num(as);
                fn(j,1)=anum;
            end
            
            chromosome={};
            chromosome(:,1)=num2cell(fn);
            chromosome(:,2)=curr_bm.Sequence;
            msg = ">> Preview of Aligment Matrix: \n";
            fprintf(msg);
            chromosome(1:5,:) %#ok<NOPRT> 
            chromosomes_cArr=[chromosomes_cArr,{chromosome}];

        end
        

        ReadAll=vertcat(chromosomes_cArr{:});
        msg = ">> Finished compiling alignment sequences...final size: " + length(ReadAll) + "\n";
        fprintf(msg);
        
        msg = ">> Beginning final processing procedures...";
        fprintf(msg)


        ReadAllstart=cell2mat(ReadAll(:,1));
        [S,N]=count_unique(ReadAllstart);
        %% Check rows that are within 1 nt of each other in guide list
        
        % For sgRNA that were designed within 1 nt wiggle rooms of each other in the
        % library, we need to be careful not to double count these when allowing a 1
        % nt wiggle room in true start position. Only true start positions are
        % considered for these sgRNA.
        
        Ap=cell2mat(guide_lib(:,4));
        c=1;
        for i=1:length(Ap)-1
            diff=Ap(i+1)-Ap(i);
            ones(c)=diff;
            c=c+1;
        end
        next=find(ones==100);
        next=next';
        c=1;
        for i=1:length(next)
            next2(c,1)=next(i);
            next2(c+1,1)=next(i)+1;
            c=c+2;
        end
        Successive=Ap(next2); %These are sgRNAs that were designed within 1 nt of each other.
        %% Generate start sites of 1 on either side (except for successive sgRNA)
        
        c=1;
        for i=1:length(Ap)
            All_pos(c,1)=Ap(i)-100;
            All_pos(c+1,1)=Ap(i);
            All_pos(c+2,1)=Ap(i)+100;
            c=c+3;
        end
        %% Match and Count start positions
        
        % start positions for all reads are compared to the start positions that we
        % determined for our library of sgRNA previously, and if they match (with a
        % 1 nt wiggle room), then the read is mapped to the corresponding sgRNA.
        % However, For sgRNA within 1 bp of each other, we need to be careful while
        % collapsing, as we could double count. Precautions are taken in collapsing.
        
        for i=1:length(All_pos)
            val=All_pos(i,1);
            ind=find(S==val);
            if isempty(ind)
                ct=0;
            else
                ct=N(ind);
            end
            All_pos(i,2)=ct;
        end
        %% Collapse the counts for each read and map it to actual sgRNA
        
        % If the sgRNA position is in the list successive, then those sgRNA only
        % have their true start site counts counted.
        
        w=2;
        for i=1:length(Ap)
            if ~ismember(Ap(i),Successive)
                Ap(i,2)=sum(All_pos(w-1:w+1,2));
                w=w+3;
            else
                Ap(i,2)=All_pos(w,2);
                w=w+3;
            end
        end
        
        BOWTIE=Ap(:,2);
        %%
        fprintf(strcat(">> Finished bowtie2 counting for sample ",sample_name," with adapter ",adapter_name,"\n"))
        load(fpath_mat_data,"GUIDE_RNA_SEQUENCE");
        bowtie_table = table(GUIDE_RNA_SEQUENCE,BOWTIE);
        head(bowtie_table)
        disp(">> Saving Bowtie2 Counts & Results Table...")
        save(fpath_mat_data,"BOWTIE","-append")
        disp(">> ...Finished Saving Results")

        bowtie_end = duration(0,0,toc(bowtie_start));
        fprintf(">> Elapsed time for bowtie2 counting: %s \n", ...
            bowtie_end)

        status=true;
        msg = sprintf(">> [%s] ...FINISHED EXECUTION(%s)\n", ...
            datetime('now',Format='default'),func_name);
        fprintf(msg)
        err="";
        
        
    catch err
        status=false;
        msg = sprintf(">> [%s] ...Failed to finish executing (%s)\n", ...
            datetime('now',Format='default'),func_name);
        fprintf(msg)
    end
end
% change .mat variable name in line 60
% change input .bam file name in line 47
% change variable name in line 169
% change output variable & file names in line 229