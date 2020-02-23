function expectation(m::UnivariateMixture; kwargs...) # uses the defaults for each dist. otherwise, you can manually construct expectations and glob them together
    dists = components(m)
    mixtureweights = probs(m)
    expectations = [expectation(d; kwargs...) for d in dists]
    return MixtureExpectation(expectations, mixtureweights)
end

function (e::MixtureExpectation)(f::Function; kwargs...)
    return dot(e.mixtureweights, [E(f) for E in e.expectations])
end

weights(e::MixtureExpectation) = e.mixtureweights
expectations(e::MixtureExpectation) = e.expectations
# *(e::MixtureExpectation, h::AbstractArray) =
# *(r::Real, e::IterableExpectation) =  IterableExpectation(nodes(e), r * weights(e)) # Necessary because, for example, multiplying UnitRange * 2 = StepRange
# expectation(f::Function, m::UnivariateMixture) = 
