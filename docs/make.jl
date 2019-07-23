using Documenter, Expectations 

# Compile the raw documentation. 
makedocs(sitename = "Expectations.jl")

# Push the documentation to the server.
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/QuantEcon/Expectations.jl.git",
)
