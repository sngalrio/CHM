function pnts = InitialPoints(polyt)
    global tol_mem;
    global tol_zero;
    global n_dec_p;
    
    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
    nA = max(size(A,2),size(Aeq,2));
    h = zeros(1,nA);
    h(dims(1)) = 1;
    % maximum
    [xopt,~,sol_flag] = cplexlp(-h,A,b,Aeq,beq,lb,ub);
    if sol_flag ~= 1
    	error('Error. \nNo feasible solution found.')
    end
    xopt=round(xopt,n_dec_p);
    pnts = ExtremePoint(h,h*xopt,-1,polyt);
    % minimum
    [xopt,~,sol_flag] = cplexlp(h,A,b,Aeq,beq,lb,ub);
    if sol_flag ~= 1
    	error('Error. \nNo feasible solution found.')
    end
    xopt=round(xopt,n_dec_p);
    ep = ExtremePoint(h,h*xopt,1,polyt);
    if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1)
        pnts = [pnts, ep];
    end   
    while size(pnts,2) <= length(dims)
        [h,h0] = GetHyperplane(pnts,dims);
        [xopt,~,sol_flag] = cplexlp(h,A,b,Aeq,beq,lb,ub);
        if sol_flag ~= 1
            error('Error. \nNo feasible solution found.')
        end
         xopt=round(xopt,n_dec_p);
         hx =h*xopt;
        if abs(hx - h0) > tol_zero
            ep = ExtremePoint(h,hx,1,polyt);
            if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1)
                 pnts = [pnts, ep];
            end
        else
            [xopt,~,sol_flag] = cplexlp(-h,A,b,Aeq,beq,lb,ub);
            if sol_flag ~= 1
                error('Error. \nNo feasible solution found.')
            end
            xopt=round(xopt,n_dec_p);
            ep = ExtremePoint(h,h*xopt,-1,polyt);
            if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1)
                 pnts = [pnts, ep];
            end
        end
    end
end

