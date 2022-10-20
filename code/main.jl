using CSV, DataFrames, Pipe, GLM, LinearAlgebra, Statistics

# read in data
cereal = CSV.read("data/cereal.csv", DataFrame)
demographics = CSV.read("data/demographics.csv", DataFrame)

# fix error in ID variable
add(x,y) = x .+ y
transform!(cereal, [:id, :quarter] => add => :id)

# calculate outside option market shares
cereal =
    @pipe cereal |>
    groupby(_, [:city, :year, :quarter]) |>
    combine(_, :share => sum => :outside_share) |>
    leftjoin(_, cereal, on = [:city, :year, :quarter])

# create unique market (city x year x quarter) IDs 
cereal.market = groupby(cereal, [:city, :year, :quarter]).groups;


# calculate equm markups
α = 30.0; # use placeholder value for price sensitivity estimate

Ωstar = OwnershpMatrix(cereal.firm)
H = ComputeElasticites(α, cereal.share, cereal.market)

Ω = Ωstar .* H

markup = -inv(Ω) * cereal.share
mc = cereal.price .- markup
margin = (cereal.price .- mc) ./ mc

SummarizeMarkups(30.0, cereal)

# generate vector of post-merger ownership
# Post-Nabisco (3,6)
# GM -Quaker (2,4)

function NewOwnership(firms, firm1, firm2)
    # generate ownership matrix that merges two firms

    for (i, firm) in enumerate(firms)
        if firm == firm2
            firms[i] = firm1
        end
    end

    return firms
end


Ωstar = OwnershpMatrix(cereal.firm) # post-nabisco merger
H = ComputeElasticites(30.0, cereal.share, cereal.market)
Ω = Ωstar .* H

# markup
μ = - inv(Ω) * cereal.share
mc = cereal.price .- μ
m = (cereal.price .- mc) ./ cereal.price

# Post-nabisco merger
merged_firms_post_nabisco = NewOwnership(cereal.firm, 3, 6)
Ω_merged = OwnershpMatrix(merged_firms_post_nabisco)


p_post_nabisco, s_post_nabisco, Ω_post_nabisco =  SimulateMerger(Ω_merged, mc, cereal.price, cereal.share, α, cereal.market)

p_post_nabisco .- cereal.price

markup


cereal.share


SummarizeMarkups(α, cereal)