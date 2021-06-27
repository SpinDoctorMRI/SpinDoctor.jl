using SpinDoctor
using Documenter

DocMeta.setdocmeta!(SpinDoctor, :DocTestSetup, :(using SpinDoctor); recursive = true)

makedocs(;
    modules = [SpinDoctor],
    authors = "Syver Døving Agdestein <syverda@icloud.com> and contributors",
    repo = "https://github.com/agdestein/SpinDoctor.jl/blob/{commit}{path}#{line}",
    sitename = "SpinDoctor.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://agdestein.github.io/SpinDoctor.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/agdestein/SpinDoctor.jl", devbranch = "main")
