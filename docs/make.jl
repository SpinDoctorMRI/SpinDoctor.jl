using SpinDoctor
using Literate
using Documenter

DocMeta.setdocmeta!(SpinDoctor, :DocTestSetup, :(using SpinDoctor); recursive = true)

# Generate examples
examples = [
    "Solve BTPDE" => "solve_btpde",
    "Custom gradients" => "custom_gradients",
    "Compare ADCs" => "compare_adcs",
    "Matrix Formalism" => "matrix_formalism",
    "High Angular Resolution" => "hardi",
]
output = "generated"
for e ∈ examples
    e = joinpath(@__DIR__, "..", "examples", "$(e.second).jl")
    o = joinpath(@__DIR__, "src", output)
    Literate.markdown(e, o)
    # Literate.notebook(e, o)
    # Literate.script(e, o)
end

makedocs(;
    modules = [SpinDoctor],
    authors = "Syver Døving Agdestein and contributors",
    repo = "https://github.com/SpinDoctorMRI/SpinDoctor.jl/blob/{commit}{path}#{line}",
    sitename = "SpinDoctor.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://spindoctormri.github.io/SpinDoctor.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        # "Examples" => [
        #     "Solve BTPDE" => "examples/solve_btpde.md",
        #     "Custom gradients" => "examples/custom_gradients.md",
        #     "Compare ADCs" => "examples/compare_adcs.md",
        #     "Matrix Formalism" => "examples/matrix_formalism.md",
        #     "High Angular Resolution" => "examples/hardi.md",
        # ],
        "Examples" => [e.first => joinpath(output, e.second * ".md") for e ∈ examples],
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

deploydocs(; repo = "github.com/SpinDoctorMRI/SpinDoctor.jl", devbranch = "main")
