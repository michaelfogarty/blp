function jacob(mval, theta2, data::Data)
    # compute jacobian matrix of the implicit function that defines mean utility
    @unpack ns, theta_i, theta_j, cdid, cdindex, vfull = data

    theta2w = Array(sparse(theta_i,theta_j,theta2));

    expmu = exp(mufunc(theta2w, data));
    shares = ind_sh(mval, expmu);
    # clear expmu

    n, K = size(x2)
    J = size(theta2w,2) - 1;
    f1 = zeros(size(cdid,1), K*(J + 1));

    # computing (∂ s_jt)/(∂σ)
    for i = 1:K
	    xv = (x2[:,i]*ones(1,ns)) .* vfull[:, ns*(i-1)+1:ns*i];    
	    temp = cumsum(xv .* shares);
	    sum1 = temp[cdindex, :];
	    sum1[2:size(sum1,1), :] = diff(sum1);
	    f1[:,i] = mean((shares .* (xv-sum1[cdid,:]))')';
	    #clear xv temp sum1
    end



end