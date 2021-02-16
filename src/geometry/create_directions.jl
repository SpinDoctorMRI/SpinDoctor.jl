""" Create ndirection directions unformly distributed in the unit circle or sphere. """
function create_directions(experiment::ExperimentSetup)
    @unpack ndirection, flat_dirs, direction = experiment

    if ndirection == 1
        # Take the one gradient direction, and normalize
        directions = direction / norm(direction)
    elseif flat_dirs
        # Create directions (unformly distributed on unit circle)
        θ = 2π * collect(0:ndirection-1)' / ndirection
        directions = [
            cos.(θ)
            sin.(θ)
            zeros(1, ndirection)
        ]
    else
        # Create Fibonacci distributed directions
        directions = create_fibonacci_sphere(ndirection)
    end
end
