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


%Original Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
% Improved by: Alexander Nguyen
% Undergraduate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Reading BAM file and storing alignments against all chromosomes

% Inputs are Bam files that are the result of aligning sgRNA from any
% particular sample against the genome using Bowtie. There will again be an
% extra chromosome here as all the nontargeting sgRNA were added as a
% seperate chromosome to allow for the aligment of the entire library.

function chromosome_startPosit=get_bowtieStartPos(fpath_bam_data)
    %%
    tic
    bamFilename = fpath_bam_data; %This file will be the result of running bowtie on 'Examplefile1_Cas9.fastqsanger' against the genome file 'CLIB89+control_Cas9.fasta'
    info = baminfo(bamFilename,'ScanDictionary',true);
    
    bm_cArr = {};
    for i = 1:length(info.ScannedDictionary)
        try
            bm_tmp = BioMap(bamFilename, 'SelectReference',info.ScannedDictionary{i});
            fprintf(">> Found Bowtie Reference: " + info.ScannedDictionary{i} + "\n");
            bm_cArr = [bm_cArr,{bm_tmp}];
        catch
            fprintf(">> Warning: Could not read reference: "+info.ScannedDictionary{i}+" (If reference = unmapped, ignore this warning)\n");
        end
    end
    fprintf(">> Amount of references found: " + length(bm_cArr)+"\n");

    %% Concatenate Start Positions and Flags
    
    % Start positions of every sgRNA in the current reads file is recorded.
    % Other information that is appended to the start position of an sgRNA is
    % the chromosome it aligned to, and the strand that it was found on. For
    % example: a start position of 2500018 is read this way: 25000|1|8. 8
    % refers to the chromosome, a flag of 0 or 1 mean Top and Bottom strands
    % respectively, and the remaining is the actual start position of the
    % sgRNA. So this sgRNA is found at the 25000 position on the bottom strand
    % in the 8th chromosome. 
    
    chromosomes_cArr={};
    for i = 1:length(bm_cArr)
        msg = ">> Starting start position search for chromosome number: " + i + "\n";
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
        bowtie_refID=cell(length(fn), 1); start_posit=cell(length(fn), 1); sequence=cell(length(fn),1);
        bowtie_refID(:)={convertCharsToStrings(info.ScannedDictionary{i})};
        start_posit(:)=num2cell(fn);
        sequence(:)=curr_bm.Sequence;
        chromosome = table(bowtie_refID,start_posit,sequence);
        head(chromosome)
        chromosomes_cArr=[chromosomes_cArr,{chromosome}];
    end
    chromosome_startPosit=vertcat(chromosomes_cArr{:});
    toc
end
