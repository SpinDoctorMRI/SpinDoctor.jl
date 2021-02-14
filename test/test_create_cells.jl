using SpinDoctor, Test

# Test create cylinders
cellsetup = CellSetup(
    name = "test",
    shape = "cylinder",
    ncell = 10,
    rmin = 1.,
    rmax = 5.,
    dmin = 0.2,
    dmax = 0.3,
    height = 3.,
)
cells = create_cells(cellsetup)
@test cells isa NamedTuple
@test size(cells.radii) == (1, cellsetup.ncell)
@test size(cells.centers) == (2, cellsetup.ncell)
@test all(cellsetup.rmin .≤ cells.radii .≤ cellsetup.rmax)


# Test create spheres
cellsetup = CellSetup(
    name = "test",
    shape = "sphere",
    ncell = 20,
    rmin = 2.,
    rmax = 7.,
    dmin = 0.2,
    dmax = 0.3,
)
cells = create_cells(cellsetup)
@test cells isa NamedTuple
@test size(cells.radii) == (1, cellsetup.ncell)
@test size(cells.centers) == (3, cellsetup.ncell)
@test all(cellsetup.rmin .≤ cells.radii .≤ cellsetup.rmax)