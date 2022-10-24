function mufunc(theta2w::Array{Float64,2}, data::Data)
    # compute nonlinear portion of utility (μ_ijt)
    @unpack x2, ns, vfull, dfull = data

    n,k = size(x2, 1), size(x2,2);
    j = size(theta2w,2)-1;
    μ = zeros(n, ns);

    for i in 1:ns
        v_i = vfull[:,i:ns:k*ns];
        d_i = dfull[:,i:ns:j*ns];
        μ[:,i] = (x2 .* v_i * theta2w[:,1]) + x2 .* (d_i * theta2w[:,2:j+1]') * ones(k,1);
    end

    return μ
end