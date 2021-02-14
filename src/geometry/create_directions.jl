""" Create ndir directions unformly distributed in the unit circle or sphere. """
function create_directions(experiment::ExperimentSetup)
    @unpack ndir, flat_dirs, direction = experiment

    if ndir == 1
        # Take the one gradient direction, and normalize
        directions = direction / norm(direction)
    elseif flat_dirs
        # Create directions (unformly distributed on unit circle)
        θ = 2π * collect(0:ndir-1)' / ndir
        directions = [
            cos.(θ)
            sin.(θ)
            zeros(1, ndir)
        ]
    else
        # Create Fibonacci distributed directions
        directions = create_fibonacci_sphere(ndir)
    end
end