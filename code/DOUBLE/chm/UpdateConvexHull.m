function UpdateConvexHull(newPt,dims)
    % Given a new extreme point, compute all possible HP with the new EP.
    global chull;
    global ePoints;
    global tol_mem;
    global n_dec_c;
    
    for i=1:size(chull,2)
        pts = chull{2,i};
	% if the new EP is already a member of the set of EPs
         if ismembertol(newPt(dims).',pts(dims,:).',tol_mem,'ByRows', true, 'DataScale',1)
             continue   
         end
         if round(chull{1,i}{1}*newPt,n_dec_c) <= round(chull{1,i}{2},n_dec_c) 
              continue
         end
        for j=1:size(pts,2)
            v = pts;
            v(:,j) = newPt;
            [h,h0] = GetHyperplane(v,dims);
            if HPinChull(h,h0,v,chull,dims,1) || HPinChull(-h,-h0,v,chull,dims,1)
               continue
            end      
            eh = round(h*ePoints,n_dec_c);
             if max(eh) <= round(h0,n_dec_c) 
                 chull = [chull, {{h,h0}; v;1}];
            else
                 if min(eh) >= round(h0,n_dec_c)
                    chull = [chull, {{-h,-h0}; v;1}];
                end
            end
        end
    end
    % Remove HP laying inside the hull
    to_remove = [];
    for j=1:size(chull,2) 
       ec = round(chull{1,j}{1}*ePoints,n_dec_c);
       h0 = round(chull{1,j}{2},n_dec_c);
       if min(ec) < h0 && max(ec) > h0
            to_remove = [to_remove j];
       end
    end
    chull(:,to_remove) = [];
end
