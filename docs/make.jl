using Documenter, Expectations 

# Compile the raw documentation. 
makedocs()

# Push the documentation to the server.
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/econtoolkit/Expectations.jl.git",
    julia = "nightly",
    osname = "linux"
)