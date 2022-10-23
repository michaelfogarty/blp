
function gmmobjg(theta2, res::Results, data::Data)
    @unpack invA, theta_i, theta_j, x1, IV = data #unpack data objects from data struct
    delta = meanval(theta2) #calls the meanval function 
    
    if maximum((delta .== NaN)) ==1
        f = 1e10
    else
        temp1 = x1'*IV
        temp2 = delta'*IV
        theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2
        res.theta1 = theta1 # save theta1 in results struct as well
        gmmresid = delta - x1*theta1
        res.gmmresid = gmmresid #this is my addition: we want to save our gmm residuals in the results struct
        temp1 = gmmresid'*IV
        f1 = temp1*invA*temp1'
        f = f1
        if nargout > 1
            # so I think that gmmobjg needs to be single valued 
            # for optimization purpposes, but we can put df into
            # the results struct without problem
            # also nargout is matlab-ese that we need to
            # get rid of here
            @unpack mval = res
            temp = jacob(mval, theta2, data)'
            df = 2*temp*IV*invA*IV'*gmmresid
        end
    end
    f, df
end