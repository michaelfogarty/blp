function cr_dumm(long_id)
    #this function creates a set of dummies for each of the values defined by long_id

    b = sort(long_id)
    b1 = [1; diff(b)] #fortunately, Julia's diff function is the same as Matlab's!
    b2 = b[b1 .> 0] #I believe we need to use element-wise comparison here
    f = sparse(zeros(length(long_id[:,1]), length(b2[:,1])))
    for i=1:length(b2[:,1])
        f[:,i] = sparse((long_id .== b2[i])) #I believe we need to compare element-wise
    end
    f #return that bad boi
end