function mktsh(mval, expmu, data::Data)
    # compute market share for each product
    @unpack ns = data

    f = sum((ind_sh(mval, expmu, data))')/ns;
    f = f'

    return f
    
end