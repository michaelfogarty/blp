## Function to estimate BLP random coefficients logit model
using CSV, Parameters, DataFrames, Kronecker, Optim, SparseArrays

include("cr_dumm.jl")

cr_dim(cereal.brand)

cereal = CSV.read("data/cereal.csv", DataFrame)
demographics = CSV.read("data/demographics.csv", DataFrame)

# Data cleaning
add(x,y) = x .+ y
transform!(cereal, [:id, :quarter] => add => :id)

# calculate outside option market shares
complement(x) = 1.0 .- x

cereal =
    @pipe cereal |>
    groupby(_, [:city, :year, :quarter]) |>
    combine(_, :share => sum => :outside_share) |>
    transform!(_, :outside_share => complement => :outside_share) |>
    leftjoin(_, cereal, on = [:city, :year, :quarter])
    

# create unique market (city x year x quarter) IDs 
cereal.market = groupby(cereal, [:city, :year, :quarter]).groups;

# recreate objects from Matlab code

@with_kw struct Data
    ns::Int64 = 20          # number of simulated individauls per market
    nmkt::Int64 = 94
    nbrn::Int64 = 24
    cdid::Vector{Int64} = kronecker(1:nmkt, ones(nbrn,1))
    cdindex::Vector{Int64} = collect(nbrn:nbrn:(nbrn*nmkt))

    # data objects

end

mutable struct Results

end

iv_vars = [string("z", x) for x in 1:20]

# iv variables - need to add brand dummies
iv = Matrix(cereal[:, iv_vars])

convert(Matrix{Int64}, Array(cr_dim(cereal.brand))) == convert(Matrix{Int64}, make_dummy_matrix(cereal.brand))