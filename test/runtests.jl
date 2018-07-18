using Expectations, Distributions, QuantEcon, Compat
using Compat.LinearAlgebra

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "Iterable distributions" begin include("iterable.jl") end 