function [h_unique] = getListHPs(dims)
     global tol_mem;
     global n_dec_p
     global chull;
     ht=[];
     for i=1:length(chull)
          ht=[ht; chull{1,i}{1}(dims) chull{1,i}{2}];
     end
     h_unique =round(uniquetol(ht,tol_mem, 'ByRows',true,'DataScale',1),n_dec_p);
end
