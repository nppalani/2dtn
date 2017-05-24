function [validReads,mapTable] = removesinglets(ecoliref,input_sam,Z)

%input_sam is actually an imported table.
% Z is the minimum number of reads required per alignment position so that
% that position is retained.


% Import SAM columns - fastqid & alignpos 
% sort input SAM alignment table by alignment position - ascending order
sorted_sam = sortrows(input_sam,'alignpos','ascend');

%Find unique alignment positions
[C, ~,~]= unique(sorted_sam.alignpos);

% Count number of times an alignment position is present (number of reads)
align_occurence = countmember(C, sorted_sam.alignpos);

%Combine C & align_occurence to trim positions with less than Z reads

minRead_positions.uniqpos = C;
minRead_positions.counts = align_occurence; 

% Convert structure to table
minRead_positionsTab = struct2table(minRead_positions);

% Sort descending
minRead_positionsTabSort = sortrows(minRead_positionsTab,'counts','descend');

% Here we will cut alignment positions that have less than Z reads

[readcounts,readcounts_index,~] = unique(minRead_positionsTabSort.counts,'stable');

CutoffIndex = find(readcounts==(Z-1));
Cutoffrow = (readcounts_index(CutoffIndex))-1;

% All positions with reads less than Z are now removed.
minRead_positionsTrimmed = minRead_positionsTabSort(1:Cutoffrow,:);

%Sort by position
minRead_positionsTrimmed = sortrows(minRead_positionsTrimmed,'uniqpos','ascend');

[FastqidIndex,~]=ismember(input_sam.alignpos,minRead_positionsTrimmed.uniqpos);

input_samStruct = table2struct(input_sam);

validReads = removeLowReads(FastqidIndex, input_samStruct);

% function removes fastqids for positions that have number of reads < Z
% (specified in input).
        function out_table = removeLowReads(fqidx,in_struct)

        x = length(fqidx);
        k=1;
            for i = 1:x
                if fqidx(i,1) == 1
                    reps(k).fastqid = in_struct(i).fastqid;
                    reps(k).alignpos = in_struct(i).alignpos;
                    k=k+1;
                end
            end
         out_table = struct2table(reps);


        end
%--------------------------------------------------------------------
 % Mapping to e coli gene IDs
%--------------------------------------------------------------------

 % Find unique positions 

 mapPos = unique(validReads.alignpos);
 
 mapTable = position_geneidmap(mapPos,ecoliref);
 
end
