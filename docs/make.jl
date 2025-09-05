using Documenter
using MomentumConservedExactDiagonalization

makedocs(
    sitename = "MomentumConservedExactDiagonalization.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Zou-Bo.github.io/MomentumConservedExactDiagonalization.jl",
        assets = String[],
    ),
    modules = [MomentumConservedExactDiagonalization],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
        "Tutorials" => "tutorials.md"
    ],
    repo = "https://github.com/Zou-Bo/MomentumConservedExactDiagonalization.jl/blob/{commit}{path}#L{line}",
    devbranch = "main",
)

deploydocs(
    repo = "github.com/Zou-Bo/MomentumConservedExactDiagonalization.jl",
    devbranch = "main",
    push_preview = true,
)