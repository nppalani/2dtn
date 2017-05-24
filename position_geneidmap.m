function out_table = position_geneidmap(J,reftable)

%J comes from removesinglets.m ; reftable comes from stripgbk.m

samtrim_height = length(J);
reftable_height = length(reftable);
k=1;


for i = 1:samtrim_height
    alignpos = J(i);
    for gh = 1:reftable_height
        
        
        startpos = reftable(gh).startpos;
        endpos = reftable(gh).stoppos;
        
        if (alignpos >= startpos) && (alignpos < endpos)
            out_struct(k).alignpos = J(i);
            out_struct(k).genename = reftable(gh).geneName;
             out_struct(k).geneproduct = reftable(gh).geneProduct;
             
             k=k+1;
        end
    end
end

out_table = struct2table(out_struct);

end
