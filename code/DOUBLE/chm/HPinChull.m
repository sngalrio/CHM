function flag = HPinChull(h,h0,v,array,dims,flag_points)
    global tol_zero;
    flag = 0;
    for i=1:size(array,2)
        % if h0 matches
        if abs(h0 - array{1,i}{2}) < tol_zero
            %if h matches
            if all(abs(h(dims) - array{1,i}{1}(dims)) < tol_zero)
                if flag_points
                    ps1 = sort(reshape(v(dims,:),1,[]));
                    ps2 = sort(reshape(array{2,i}(dims,:),1,[]));
		    % if generating points match
                    if all(abs(ps1 - ps2) < tol_zero)
                        flag = 1;
                        break
                    end
                else
                    flag = 1;
                    break
                end
            end
        end
    end
end
