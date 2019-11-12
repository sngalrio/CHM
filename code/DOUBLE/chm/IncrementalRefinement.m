function IncrementalRefinement(polyt)
    % Initial convex hull is refined by maximizing/minimizing the hyperplanes
    % containing the extreme points until all the facets of the projection are terminal.
    global chull;
    global ePoints;
    global n_dec_c;
    global n_dec_p;
    global tol_mem;
    global tol_zero;

    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};

    while sum(cell2mat(chull(end,:))) ~= 0
        i = 0;
        while i < size(chull,2)
            i = i+1;
            h = chull{1,i}{1};
            h0 = chull{1,i}{2};
            % maximize HP
            [xopt,~,sol_flag] = cplexlp(-h,A,b,Aeq,beq,lb,ub);
            if sol_flag ~= 1
                error('Error. \nNo feasible solution found for HP in CHULL.')
            end

            hx = round(h*xopt,n_dec_p);
            % If HP is terminal
            if abs(hx - h0) < tol_zero
                chull{end,i} = 0;
            else
            % if HP not terminal, compute new EP
                ep = ExtremePoint(h, hx,-1,polyt);
                if ~ismembertol(ep(dims).',ePoints(dims,:).',tol_mem,'ByRows', true,'DataScale',1)
                    ePoints = [ePoints, ep];
                    % update CH with new EP
                    UpdateConvexHull(ep,dims);
                end
            end
        end
	% Remove HP laying inside the hull
    	to_remove = [] ;
        for j=1:size(chull,2)
           ec = round(chull{1,j}{1}*ePoints,n_dec_c);
           h0 = round(chull{1,j}{2},n_dec_c);
           if min(ec) < h0  && max(ec) > h0
                to_remove = [to_remove j];
           end
        end
        chull(:,to_remove) = [];
    end
end
