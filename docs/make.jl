using Correlation, Documenter

makedocs(;
    modules=[Correlation],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://gitlab.invenia.ca/invenia/Correlation.jl/blob/{commit}{path}#L{line}",
    sitename="Correlation.jl",
    authors="Eric Davies",
    assets=[
        "assets/invenia.css",
    ],
    checkdocs = :exports,
    strict = true,
)
