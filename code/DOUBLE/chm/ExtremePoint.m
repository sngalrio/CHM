function xopt = ExtremePoint(hyrp,hyp0,opti,polyt)
    % compute an extreme point for a given projection
    global n_dec_p;
    global tol_lp;
    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
    nA = max(size(A,2),size(Aeq,2));
    A=[A;hyrp;-hyrp];
    b=[b;hyp0;-hyp0];
    h = zeros(1,nA); 
    h(dims) = ones(1,length(dims));
    [xopt,~,sol_flag,output] = cplexlp(opti*h,A,b,Aeq,beq,lb,ub);
    if sol_flag ~= 1
        'adding tolerance to EXTREME POINT OPTIMIZATION. '
        [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
        A=[A;round(hyrp,n_dec_p);-round(hyrp,n_dec_p)];
        b=[b;hyp0  + tol_lp;-hyp0 + tol_lp];
        [xopt,~,sol_flag,output] = cplexlp(opti*h,A,b,Aeq,beq,lb,ub);
        if sol_flag ~= 1
            output
            error('Error. \nNo feasible solution found EXTREME POINT.')
        end
    end
    xopt = round(xopt,n_dec_p);
end
