###need to import mvalold, ps2, gmmresid from structs
###To do: we need to create results struct to hold these things! 

#also: need to have a data struct holding x1 and IV (and potentially other stuff)!

#function accepts theta2 and results struct
function var_cov(theta2, res::Results)
    @unpack mvalold, ps2, gmmresid = res
    @unpack x1 = data
    N = length(x1[:,1]) #guessing we want the number of observations in x1
    Z = length(IV[:,2]) #I'm guessing we want the number of instruments
    temp = jacob(mvalold, theta2) #call the jacobian function and see what's up
    a = [x1 temp]' * IV
    IVres = IV .* (gmmresid * ones(Z))
    b = IVres'*IVres 

    f = inv(a*invA*a')*a*invA*b*invA*a'*inv(a*invA*a') #I am assuming we want this boi returned
    f #return :)
end