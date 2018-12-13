using Expectations, Distributions 

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test, Compat.Random
else
    using Test, Random, Statistics, LinearAlgebra
end

@testset "Iterable distributions" begin include("iterable.jl") end
