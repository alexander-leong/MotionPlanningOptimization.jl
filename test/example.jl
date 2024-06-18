# Example usage of trajopt.jl
# MIT License
# Alexander Leong

using GLMakie
using GeometryBasics

#
# Setup problem
#

# set obstacles
obstacle₁ = Box(-10.00, 0.00, 0.00, 20.00, 10.00, 30.00)
obstacle₂ = Box(0.00, 15.00, -10.00, 5.00, -10.00, 20.00)
obstacle₃ = Box(15.00, 25.00, -10.00, 5.00, 10.00, 20.00)

# set safe regions using box constraints
boxes = Array{Box}(UndefInitializer(), 7)
box₁ = Box(-10.00, 0.00, -10.00, 20.00, -10.00, 10.00)
box₂ = Box(0.00, 25.00, 5.00, 20.00, -10.00, 30.00)
box₃ = Box(-10.00, 25.00, 5.00, 20.00, -10.00, 10.00)
box₄ = Box(-10.00, 0.00, -10.00, 0.00, -10.00, 30.00)
box₅ = Box(15.00, 25.00, -10.00, 20.00, -10.00, 10.00)
box₆ = Box(0.00, 25.00, -10.00, 20.00, 20.00, 30.00)
box₇ = Box(-10.00, 25.00, -10.00, 5.00, 20.00, 30.00)
boxes[1] = box₁
boxes[2] = box₂
boxes[3] = box₃
boxes[4] = box₄
boxes[5] = box₅
boxes[6] = box₆
boxes[7] = box₇

X₀ = Dict()
X₀[:x], X₀[:y], X₀[:z] = -5.0, -5.0, -5.0
X₀′ = Dict()
X₀′[:x], X₀′[:y], X₀′[:z] = 1e0, 1e0, 1e0
X₀″ = Dict()
X₀″[:x], X₀″[:y], X₀″[:z] = 1, 1, 1
X₁ = Dict()
X₁[:x], X₁[:y], X₁[:z] = 22.5, -7.5, 0.0
X₁′ = Dict()
X₁′[:x], X₁′[:y], X₁′[:z] = 1e0, 1e0, 1e0
X₁″ = Dict()
X₁″[:x], X₁″[:y], X₁″[:z] = 1, 1, 1
initial_state = State(X₀, X₀′, X₀″)
final_state = State(X₁, X₁′, X₁″)

N = 4 # number of trajectory segments
Mxl = -10.0
Mxu = 25.0
Myl = -10.0
Myu = 20.0
Mzl = -10.0
Mzu = 30.0
params = OptimizationParameters(Mxl, Mxu, Myl, Myu, Mzl, Mzu)
model, solve_time, optimal_path = solve_for_optimal_path(N, initial_state, final_state, boxes, params=params)

#
# Render result with Makie.jl
#

fig = Figure()
ax = Axis3(
    fig[1, 1],
    aspect = :data,
    viewmode = :fit
    )
ax.azimuth = pi/4
ax.elevation = pi/4

theta = 0
C = [
    1          0           0
    0 cos(theta) -sin(theta)
    0 sin(theta)  cos(theta)
]

function render_box(b::Box, color=:blue)
  offset_x=0
  offset_y=0
  offset_z=0
  points = hcat(
      C * [b.xu+offset_x, b.yu+offset_y, b.zu+offset_z],
      C * [b.xl+offset_x, b.yu+offset_y, b.zu+offset_z],
      C * [b.xl+offset_x, b.yl+offset_y, b.zu+offset_z],
      C * [b.xu+offset_x, b.yl+offset_y, b.zu+offset_z],
      C * [b.xu+offset_x, b.yu+offset_y, b.zl+offset_z],
      C * [b.xl+offset_x, b.yu+offset_y, b.zl+offset_z],
      C * [b.xl+offset_x, b.yl+offset_y, b.zl+offset_z],
      C * [b.xu+offset_x, b.yl+offset_y, b.zl+offset_z]
  )
  
  faces = [
      1 2 3
      3 4 1
      1 2 6
      1 5 6
      1 4 5
      4 5 8
      3 4 7
      4 7 8
      5 6 7
      5 7 8
      2 3 6
      3 6 7
  ]
  
  mesh!(ax, points, faces, alpha=0.25, color=color, transparency=true, visible=true)
end

# show obstacles in scene
render_box(obstacle₁)
render_box(obstacle₂)
render_box(obstacle₃)

# uncomment to render safe regions
# render_box(box₁, :yellow)
# render_box(box₂, :yellow)
# render_box(box₃, :yellow)
# render_box(box₄, :yellow)
# render_box(box₅, :yellow)
# render_box(box₆, :yellow)
# render_box(box₇, :yellow)

# show trajectory
scatterlines!(optimal_path.x, optimal_path.y, optimal_path.z, markersize=5, markercolor = :orange)

fig