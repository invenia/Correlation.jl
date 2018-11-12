using Documenter, Correlation

makedocs(;
    modules=[Correlation],
    format=:html,
    pages=[
        "Home" => "index.md",
    ],
    repo="https://gitlab.invenia.ca/invenia/Correlation.jl/blob/{commit}{path}#L{line}",
    sitename="Correlation.jl",
    authors="Eric Davies",
    assets=[
        "assets/invenia.css",
        "assets/logo.png",
    ],
    strict = true,
    checkdocs = :none,
    html_prettyurls = false,
)
