# MotionPlanningOptimization.jl

MotionPlanningOptimization.jl is a Julia package that uses mathematical optimization techniques to calculate an optimal path for a robot or vehicle to move from an initial position to a final position with a given set of obstacles.

## Installation

Please ensure you have a MOSEK license before using this library.

```sh
import Pkg
Pkg.add("MotionPlanningOptimization")
```

You will need these additional dependencies to run the notebook.jl example.
```sh
import Pkg
Pkg.add("GeometryBasics")
Pkg.add("GLMakie")
```

## Usage

```sh
using MotionPlanningOptimization

# number of trajectory segments
N = 4
# initial position
X₀ = Dict()
X₀[:x], X₀[:y], X₀[:z] = -5.0, -5.0, -5.0
# initial velocity
X₀′ = Dict()
X₀′[:x], X₀′[:y], X₀′[:z] = 1e0, 1e0, 1e0
# initial acceleration
X₀″ = Dict()
X₀″[:x], X₀″[:y], X₀″[:z] = 1, 1, 1
# final position
X₁ = Dict()
X₁[:x], X₁[:y], X₁[:z] = 22.5, -7.5, 0.0
# final velocity
X₁′ = Dict()
X₁′[:x], X₁′[:y], X₁′[:z] = 1e0, 1e0, 1e0
# final acceleration
X₁″ = Dict()
X₁″[:x], X₁″[:y], X₁″[:z] = 1, 1, 1
initial_state = State(X₀, X₀′, X₀″)
final_state = State(X₁, X₁′, X₁″)
Mxl, Mxu, Myl, Myu, Mzl, Mzu = -10.0, 25.0, -10.0, 20.0, -10.0, 30.0
# set safe regions
num_boxes = 1
boxes = Array{Box}(UndefInitializer(), num_boxes)
box₁ = Box(-10.00, 0.00, -10.00, 20.00, -10.00, 10.00)
boxes[1] = box₁

# solve
params = OptimizationParameters(Mxl, Mxu, Myl, Myu, Mzl, Mzu)
model, solve_time, optimal_path = solve_for_optimal_path(N, initial_state, final_state, boxes, params=params)
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please apply at https://opencollective.com/alexander-leong if you would like to sponsor this work.

## License

[MIT](https://choosealicense.com/licenses/mit/)