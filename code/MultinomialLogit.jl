function OwnershpMatrix(firm)
    # Compute ownership matrix given vector of firms that own each product

    n = length(firm)
    Ω = zeros(n,n)

    for i in 1:n, j in 1:n
        Ω[i,j] = firm[i] == firm[j]
    end

    return Ω
end

function ComputeElasticites(α::Float64, shares, market)
    # compute matrix of own and cross price elasticites given estimate of α
    # valid for multinomial logit demand system w/o random coefficients
    # shares is a vector of market shares
    # market gives a vector of market IDs

    n = length(shares)
    H = zeros(n,n)

    for i in 1:n, j in 1:n
        if market[i] == market[j]
        # only compute elasticities for goods in same market
            if i != j
            # cross price elasticity
                H[i,j] = α * shares[i] * shares[j]
            elseif i == j
            # own price elasticity
                H[i,j] = -α* shares[i] *(1.0 - shares[j])
            end
        end
    end

    return H
end

function SummarizeMarkups(α::Float64, cereal::DataFrame)
    
    # compute omega matrix 
    Ωstar = OwnershpMatrix(cereal.firm)
    H = ComputeElasticites(α, cereal.share, cereal.market)
    Ω = Ωstar .* H

    markup = - inv(Ω) * cereal.share
    mc = cereal.price .- markup
    margin = (cereal.price .- mc) ./ mc

    means = [mean(x) for x in [markup, mc, margin]]
    medians = [median(x) for x in [markup, mc, margin]]
    stdevs = [std(x) for x in [markup, mc, margin]]

    return [means,
            medians,
            stdevs]
end

function SimulateMerger(own, mc, p0, s0, α, market; tol = 1e-8, err = 100, max_iter = 10000)
    s = s0 # initial vector of shares
    i = 0

    while err > tol && i < max_iter
        i += 1

        # update Ω based on last loop's shares
        H = ComputeElasticites(α, s, market)
        Ω = own .* H

        # calculate new price based on updated Ω
        μ = -inv(Ω) * s # compute markup
        p = mc .+ μ     # compute new prices
        s = -Ω*(p0-mc)        # compute market shares implied by markups

        err = maximum(abs.(p0-p))
        p0 = p

    end

    println(i, " iterations")

    return p0, s, Ω
end

