# Render the images of the examples gallery (also shown in the README) with CairoMakie:
#
#     julia --project=examples docs/gallery.jl [example ...]
#
# Each example is run, and the n-th figure it shows (the value of a statement that is not
# an assignment) is saved in `docs/src/public/gallery`.

using CairoMakie

images = [
    "quickstart" => 2,
    "geometries" => 9,
    "white_matter" => 1,
    "periodic_tissue" => 1,
    "neuron" => 3,
    "compare_adcs" => 2,
    "matrix_formalism" => 2,
    "custom_gradients" => 2,
    "validation" => 2,
]
names = isempty(ARGS) ? first.(images) : ARGS
outdir = joinpath(@__DIR__, "src", "public", "gallery")
mkpath(outdir)

for (name, index) in filter(i -> first(i) in names, images)
    file = joinpath(@__DIR__, "..", "examples", "$name.jl")
    mod = Module(Symbol(name))
    Core.eval(mod, :(include(x) = Base.include($mod, x)))
    code = Meta.parseall(read(file, String); filename = file)
    nfigure = 0
    for expr in code.args
        expr isa LineNumberNode && continue
        value = Core.eval(mod, expr)
        if value isa Makie.FigureLike && !Meta.isexpr(expr, :(=))
            nfigure += 1
            nfigure == index && save(joinpath(outdir, "$name.png"), value; px_per_unit = 1)
        end
    end
    @info "Saved the image of $name"
end
