###need to import mvalold, ps2, gmmresid from structs
###To do: we need to create results struct to hold these things! 

#also: need to have a data struct holding x1 and IV (and potentially other stuff)!

#function accepts theta2 and results struct
function var_cov(theta2, res::Results, data::Data)
    @unpack mvalold, ps2, gmmresid = res
    @unpack x1, IV = data

    N,Z = size(x1,1), size(IV, 2)
    temp = jacob(mval, theta2, data) #call the jacobian function and see what's up
    a = [x1 temp]' * IV
    IVres = IV .* (gmmresid * ones(Z))
    b = IVres'*IVres 

    f = inv(a*invA*a')*a*invA*b*invA*a'*inv(a*invA*a') #I am assuming we want this boi returned
    return f #return :)
end