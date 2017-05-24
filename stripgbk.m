function Stripped_struct = stripgbk(gbkstruct)

%Download E. coli genbank file from NCBI and use genbankread to load into a
%structure. Pass that structure to this script to keep only the elements
%mentioned in Stripped_struct. Output is a 

[~,x]= size(gbkstruct.CDS);
k=1;

Stripped_struct = struct('startpos','stoppos','geneName','geneProduct');
for i = 1:x
    if length(gbkstruct.CDS(i).indices) == 2
        
        pos = gbkstruct.CDS(i).indices;
            
        if pos(1) > pos(2) 
            a = pos(1); pos(1) = pos(2); pos(2) = a;
        end
        
        Stripped_struct(k).startpos = pos(1);
        Stripped_struct(k).stoppos = pos(2);
        Stripped_struct(k).geneName = gbkstruct.CDS(i).gene;
        Stripped_struct(k).geneProduct = gbkstruct.CDS(i).product;
        k=k+1;
    end
    
end
