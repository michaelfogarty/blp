function meanval(theta2, res::Results, data::Data)
    # compute the market shares using the BLP contraction
    @unpack theta_i, theta_j, x2, s_jt = data

    # whats the difference between theta2 and theta2w?
    # we need to calculate the difference b/t current estimates (theta2) and
    # existing estimates

    if max(abs.(theta2-res.theta2)) < 0.01;
        tol = 1e-9;
        flag = 0;
    else
        tol = 1e-6;
        flag = 1;
    end

    # here we need to turn the theta2 vector into a matrix
    theta2w = Array(sparse(theta_i,theta_j,theta2));
    expmu = exp(mufunc(x2, theta2w, data));
    norm = 1;
    avgnorm = 1;
    
    i=0;
    mvalold = res.mval;

    while norm > tol*10^(flag*floor(i/50)) & avgnorm > 1e-3*tol*10^(flag*floor(i/50))
        # the following two lines are equivalent; however, the latter saves on the number of exponents
        # computed at therefore speeds up the computation by 5-10%
        #  mval = mvalold + log(s_jt) - log(mktsh(mvalold,expmu));
        mval = mvalold .* s_jt ./ mktsh(mvalold, expmu); 
        t = abs.(mval-mvalold);
        norm = max(t);
        avgnorm = mean(t);
        mvalold = mval;
        i +=  1; #update loop counter
    end

    if flag == 1 && max(isnan.(mval)) < 1 # update second condition here
        res.mval = mval
        res.theta2 = theta2 # do we really want to be saving this here? worth considering
    end
    
    return log(mval)
     

end