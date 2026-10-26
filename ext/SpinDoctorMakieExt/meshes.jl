# Meshes and fields on the compartments

# Default colors of the compartments
compartment_colors(n) = [Makie.wong_colors()[mod1(i, 7)] for i = 1:n]

# Axis for a mesh in 2D or 3D, with equal scaling of the axes
function mesh_axis(position, dim; title = "", decorations = true)
    if dim == 2
        ax = Axis(position; title, aspect = DataAspect())
    else
        ax = Axis3(position; title, aspect = :data)
    end
    decorations || (hidedecorations!(ax); hidespines!(ax))
    ax
end

# The faces to draw for each compartment: the elements in 2D, and the boundary facets in
# 3D, each facet once (an interface is drawn with the first of its compartments shown)
function visible_faces(femesh::FEMesh{T,dim}, compartments) where {T,dim}
    dim == 2 && return femesh.elements
    ncompartment, nboundary = size(femesh.facets)
    faces = [zeros(Int, 3, 0) for _ = 1:ncompartment]
    for j = 1:nboundary
        i = findfirst(i -> !isempty(femesh.facets[i, j]), compartments)
        isnothing(i) && continue
        icmpt = compartments[i]
        faces[icmpt] = [faces[icmpt] femesh.facets[icmpt, j]]
    end
    faces
end

# All the shown compartments as one mesh (so that the triangles are sorted by depth in
# CairoMakie), with the point indices of each compartment in the global field
function merged_mesh(femesh, compartments)
    faces = visible_faces(femesh, compartments)
    offsets = cumsum([0; size.(femesh.points, 2)])
    points = reduce(hcat, femesh.points[compartments])'
    shown = cumsum([0; size.(femesh.points[compartments], 2)])
    f = reduce(hcat, [faces[i] .+ shown[k] for (k, i) in enumerate(compartments)])'
    indices = reduce(vcat, [(offsets[i]+1):offsets[i+1] for i in compartments])
    points, f, indices
end

function plot_mesh!(
    ax,
    femesh::FEMesh{T,dim};
    compartments = 1:length(femesh.points),
    colors = compartment_colors(length(femesh.points)),
    alpha = 1.0,
    strokewidth = dim == 2 ? 0.5 : 0.0,
    strokecolor = (:black, 0.5),
) where {T,dim}
    alphas = alpha isa Number ? fill(alpha, length(femesh.points)) : alpha
    points, faces, _ = merged_mesh(femesh, compartments)
    color = reduce(
        vcat,
        map(compartments) do i
            c = Makie.to_color(colors[i])
            fill(
                Makie.RGBAf(c.r, c.g, c.b, alphas[i] * c.alpha),
                size(femesh.points[i], 2),
            )
        end,
    )
    if strokewidth > 0
        poly!(ax, points, faces; color, strokewidth, strokecolor)
    else
        mesh!(ax, points, faces; color)
    end
end

function plot_mesh(
    femesh::FEMesh{T,dim};
    title = "",
    decorations = true,
    size = (600, 500),
    kwargs...,
) where {T,dim}
    fig = Figure(; size)
    ax = mesh_axis(fig[1, 1], dim; title, decorations)
    plot_mesh!(ax, femesh; kwargs...)
    fig
end

# Values of a complex field to plot
function field_values(u, part)
    part == :abs && return abs.(u)
    part == :real && return real.(u)
    part == :imag && return imag.(u)
    part == :phase && return angle.(u)
    error("Unknown part $part (use :abs, :real, :imag or :phase)")
end

field_label(part) =
    Dict(:abs => "|M|", :real => "Re M", :imag => "Im M", :phase => "arg M")[part]

function default_colorrange(values, part)
    part == :phase && return (-pi, pi)
    vmax = maximum(abs, values)
    part == :abs ? (zero(vmax), vmax) : (-vmax, vmax)
end

default_colormap(part) =
    part == :abs ? :viridis : part == :phase ? :cyclic_mrybm_35_75_c68_n256 : :RdBu

function plot_field!(
    ax,
    femesh::FEMesh{T,dim},
    u;
    compartments = 1:length(femesh.points),
    part = :abs,
    colormap = default_colormap(part),
    colorrange = nothing,
) where {T,dim}
    points, faces, indices = merged_mesh(femesh, compartments)
    uobs = u isa Observable ? u : Observable(u)
    color = lift(u -> field_values(u[indices], part), uobs)
    isnothing(colorrange) && (colorrange = default_colorrange(color[], part))
    mesh!(ax, points, faces; color, colormap, colorrange, shading = NoShading)
end

function plot_field(
    femesh::FEMesh{T,dim},
    u;
    part = :abs,
    title = "",
    decorations = true,
    size = (650, 500),
    kwargs...,
) where {T,dim}
    fig = Figure(; size)
    ax = mesh_axis(fig[1, 1], dim; title, decorations)
    plot = plot_field!(ax, femesh, u; part, kwargs...)
    Colorbar(fig[1, 2], plot; label = field_label(part))
    fig
end
