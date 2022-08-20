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
*(r::Real, e::MixtureExpectation) = MixtureExpectation(r * expectations(e), weights(e))
# *(e::MixtureExpectation, h::AbstractArray) = dot(weights(e), [E*h for E in expectations(e)])
import Base.+
+(expectations::Expectation...) = MixtureExpectation(expectations, ones(length(expectations)))

import Base.*

# Right-multiplying an expectation by something.
"""
    *(e::MixtureExpectation, h::AbstractArray) = dot(map(x -> x * h, expectations(e)), weights(e))
Implements the right-application of an `MixtureExpectation` by a vector of values on its nodes.
"""
*(e::MixtureExpectation, h::AbstractArray) = dot(map(x -> x * h, expectations(e)), weights(e))

# Left-multiplying an expectation by a scalar.
"""
    *(r::Real, e::MixtureExpectation) =  IterableExpectation(nodes(e), r * weights(e))
Implements left-multiplication of an `IterableExpectation` by a real scalar.
"""
*(r::Real, e::MixtureExpectation) = MixtureExpectation(r * expectations(e), weights(e)) # Necessary because, for example, multiplying UnitRange * 2 = StepRange
