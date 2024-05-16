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
        fpath_mat_data = strcat(NGS_SETTINGS.mat_workspace,collection,"/",sample_name,"/",sample_name,"_",adapter_name,".mat");
        
        %%
        fprintf(strcat(">> Starting bowtie2 counting on sample ",sample_name," with adapter ", adapter_name,"\n"))
        
        bowtie_start = tic;
        bamFilename = fpath_bam_dataset; %This file will be the result of running bowtie on 'Examplefile1_Cas9.fastqsanger' against the genome file 'CLIB89+control_Cas9.fasta'
        info = baminfo(bamFilename,'ScanDictionary',true);
        bmA = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{1});
        bmB = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{2});
        bmC = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{3});
        bmD = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{4});
        bmE = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{5});
        bmF = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{6});
        if length(info.ScannedDictionary)==9
            bmH = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{8});
        else
            bmH = BioMap(bamFilename, 'SelectReference', info.ScannedDictionary{7});
        end
        
        load(fpath_guide_library,"guide_table");
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
        
        Astart=bmA.Start; %Start positions of sgRNA aligning to the A chromosome
        Aflag=bmA.Flag;
        
        Bstart=bmB.Start;
        Bflag=bmB.Flag;
        
        Cstart=bmC.Start;
        Cflag=bmC.Flag;
        
        Dstart=bmD.Start;
        Dflag=bmD.Flag;
        
        Estart=bmE.Start;
        Eflag=bmE.Flag;
        
        Fstart=bmF.Start;
        Fflag=bmF.Flag;
        
        Hstart=bmH.Start;
        Hflag=bmH.Flag;
        
        Aflag(Aflag==16)=1; %The reverse strand alignment to a chromosome (for chexample Chromosome A here) is indicated as the number 16 on Bowtie, that is chnaged to 1 here
        Bflag(Bflag==16)=1;
        Cflag(Cflag==16)=1;
        Dflag(Dflag==16)=1;
        Eflag(Eflag==16)=1;
        Fflag(Fflag==16)=1;
        Hflag(Hflag==16)=1;
        
        for i=1:length(Astart)
            as=sprintf('%1d%1d%1d ',Astart(i),Aflag(i),1); %Astart is start position, Aflag is flag (16->1)
            a1=str2num(as); %str2num is done so allocation is not messed up
            Afn(i,1)=a1;
        end
        for i=1:length(Bstart)
            as=sprintf('%1d%1d%1d ',Bstart(i),Bflag(i),2);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Bfn(i,1)=a1;
        end
        for i=1:length(Cstart)
            as=sprintf('%1d%1d%1d ',Cstart(i),Cflag(i),3);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Cfn(i,1)=a1;
        end
        for i=1:length(Dstart)
            as=sprintf('%1d%1d%1d ',Dstart(i),Dflag(i),4);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Dfn(i,1)=a1;
        end
        for i=1:length(Estart)
            as=sprintf('%1d%1d%1d ',Estart(i),Eflag(i),5);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Efn(i,1)=a1;
        end
        for i=1:length(Fstart)
            as=sprintf('%1d%1d%1d ',Fstart(i),Fflag(i),6);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Ffn(i,1)=a1;
        end
        for i=1:length(Hstart)
            as=sprintf('%1d%1d%1d ',Hstart(i),Hflag(i),8);
            a1=str2num(as); %str2num is done so allocation is not messed up
            Hfn(i,1)=a1;
        end
        
        A(:,1)=num2cell(Afn);
        A(:,2)=bmA.Sequence;
        
        B(:,1)=num2cell(Bfn);
        B(:,2)=bmB.Sequence;
        
        C(:,1)=num2cell(Cfn);
        C(:,2)=bmC.Sequence;
        
        D(:,1)=num2cell(Dfn);
        D(:,2)=bmD.Sequence;
        
        E(:,1)=num2cell(Efn);
        E(:,2)=bmE.Sequence;
        
        F(:,1)=num2cell(Ffn);
        F(:,2)=bmF.Sequence;
        
        H(:,1)=num2cell(Hfn);
        H(:,2)=bmH.Sequence;
        
        ReadAll=vertcat(A,B,C,D,E,F,H);
        ReadAllstart=cell2mat(ReadAll(:,1));
        [S,N]=count_unique(ReadAllstart);
        %% Check rows that are within 1 nt of each other in guide list
        
        % For sgRNA that were designed within 1 nt wiggle rooms of each other in the
        % library, we need to be careful not to double count these when allowing a 1
        % nt wiggle room in true start position. Only true start positions are
        % considered for these sgRNA.
        
        Ap=cell2mat(guide_table(:,4));
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