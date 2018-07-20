using Documenter, Expectations 

# Compile the raw documentation. 
makedocs()

# Push the documentation to the server.
deploydocs(
    repo = "github.com/econtoolkit/Expectations.jl.git",
    julia = "nightly"
)