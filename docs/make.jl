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
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Examples" => [
            "Solve BTPDE" => "examples/solve_btpde.md",
            "Custom gradients" => "examples/custom_gradient.md",
            "Compare ADCs" => "examples/compute_adc.md",
            "High Angular Resolution" => "examples/hardi.md",
        ],
        "Theory" => [
            "Theory" => "theory/theory.md",
            "Bloch-Torrey PDE" => "theory/btpde.md",
            "Apparent Diffusion Coefficient" => "theory/adc.md",
            "Matrix Formalism" => "theory/matrix_formalism.md",
            "Discretization" => "theory/discretization.md",
        ],
        "Neurons" => "neurons.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/agdestein/SpinDoctor.jl", devbranch = "main")
