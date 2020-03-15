# uses the defaults for each dist.
# otherwise, you can manually construct expectations and glob them together
function expectation(m::UnivariateMixture; kwargs...)
    expectations = [expectation(d; kwargs...) for d in components(m)]
    return MixtureExpectation(expectations, probs(m))
end

function (e::MixtureExpectation)(f::Function; kwargs...)
    return dot(e.mixtureweights, [E(f; kwargs...) for E in e.expectations])
end

weights(e::MixtureExpectation) = e.mixtureweights
expectations(e::MixtureExpectation) = e.expectations
expectation(f::Function, m::UnivariateMixture; kwargs...) = dot(probs(m), [expectation(f, dist; kwargs...) for dist in components(m)])
*(r::Real, e::MixtureExpectation) = MixtureExpectation(r*expectations(e), weights(e))
# *(e::MixtureExpectation, h::AbstractArray) = dot(weights(e), [E*h for E in expectations(e)])
import Base.+
+(e1::IterableExpectation, e2::IterableExpectation) = MixtureExpectation([e1, e2], [1.0, 1.0])
+(expectations...) = MixtureExpectation(expectations, ones(length(expectations)))
