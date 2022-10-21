
function gmmobjg(theta2, res::Results, data::Data)
    @unpack invA, theta1, theti, thetj, x1, IV = data #unpack data objects from data struct
    delta = meanval(theta2) #calls the meanval function 
    
    if maximum((delta .== NaN)) ==1
        f = 1e10
    else
        temp1 = x1'*IV
        temp2 = delta'*IV
        theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2
        gmmresid = delta - x1*theta1
        res.gmmresid = gmmresid #this is my addition: we want to save our gmm residuals in the results struct
        temp1 = gmmresid'*IV
        f1 = temp1*invA*temp1'
        f = f1
        if nargout > 1
            @unpack mvalold = res
            temp = jacob(mvalold, theta2)'
            df = 2*temp*IV*invA*IV'*gmmresid
        end
    end
    f, df
end