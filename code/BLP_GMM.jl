## Function to estimate BLP random coefficients logit model
using CSV, Parameters, DataFrames, Kronecker, Optim, SparseArrays, Random, Distributions, Pipe
Random.seed!(104) # set seet for RNG repeatability

include("cr_dumm.jl")


cereal = CSV.read("data/cereal.csv", DataFrame)
demographics = CSV.read("data/demog.csv", DataFrame, header = false)
v = CSV.read("data/v.csv", DataFrame)

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


@with_kw struct Data
    # Indexes
    ns::Int64 = 20          # number of simulated individauls per market
    nmkt::Int64 = 94        # number of markets
    nbrn::Int64 = 24        # number of brands
    cdid::Vector{Int64} = repeat(1:nmkt, inner = nbrn)      # market id
    cdindex::Vector{Int64} = collect(nbrn:nbrn:(nbrn*nmkt)) # index of last obs for each market

    # data matricies
    iv_vars::Vector{String} = [string("z", x) for x in 1:20]
    iv::Array{Float64,2} = Matrix(cereal[:, iv_vars]);
    x1::Array{Float64, 2} = hcat(cereal.price, Array(cr_dumm(cereal.firmbr)));
    x2::Array{Float64, 2} = hcat(ones(size(cereal,1), 1), cereal.price, cereal.sugar, cereal.mushy)

    # market shares
    s_jt::Array{Float64,1} = cereal.share

    # IV variables
    IV::Array{Float64, 2} = hcat(iv, x1[:, 2:end]);
    invA::Array{Float64, 2} = inv(IV' * IV);

    # demographics and random draws
    v::Array{Float64, 2} = Matrix(v[:, [string("v", x) for x in 1:80]])
    vfull::Array{Float64, 2} = repeat(v, inner = (nbrn, 1), outer = (1,1))
    d::Array{Float64, 2} = Matrix(demographics)
    dfull::Array{Float64, 2} = repeat(d, inner = (nbrn, 1), outer = (1,1))

    # parameters in the θ matrix - Initial values
    theta2w::SparseMatrixCSC = sparse([0.3772    3.0888         0    1.1859         0;
                                  1.8480   16.5980    -.6590         0   11.6245;
                                  -0.0035  -0.1925         0    0.0296         0;
                                  0.0810    1.4684         0   -1.5143         0]);
    #theta2w = sparse(theta2w)
    A::Tuple{Vector{Int64}, Vector{Int64}, Vector{Float64}} = findnz(theta2w)
    theta_i::Vector{Int64} = A[1]
    theta_j::Vector{Int64} = A[2]
    theta2::Vector{Float64} = A[3]


end

mutable struct Results
    # Initial values of theta2
    theta2::Array{Float64,1} # nonlinear coefficients
    theta1::Array{Float64,1} # coefficients of linear part of the model

    mval::Array{Float64, 1} # 
    gmmresid::Array{Float64,1} # gmma residuals
end



function Initialize(data::Data)
    # initialize results struct
    @unpack ns, nbrn, nmkt = data

    # Initial parameter values
    theta2 = data.theta2
    theta1 = zeros(nbrn + 1) # price + number of brand dummies

    mval = data.s_jt
    gmmresid = zeros(nbrn*nmkt)

    Results(theta2, theta1, mval, gmmresid)
end

include("meanval.jl")
include("mufunc.jl")
include("mktsh.jl")
include("ind_sh.jl")

data = Data();
res = Initialize(data);
theta2 = randn(13)

meanval(theta2, res, data)