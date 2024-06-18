module TrajOpt

using DynamicPolynomials
using Hypatia
using JuMP
using LinearAlgebra
using MosekTools
using Pajarito
using SumOfSquares

struct Box
  xl::Float64
  xu::Float64
  yl::Float64
  yu::Float64
  zl::Float64
  zu::Float64
end

struct State
  X::Dict{Symbol, Float64}
  X′::Dict{Symbol, Float64}
  X″::Dict{Symbol, Float64}
end

struct BigM_Parameters
  Mxl::Float64
  Mxu::Float64
  Myl::Float64
  Myu::Float64
  Mzl::Float64
  Mzu::Float64
end

struct OptimalPath
  x::Array{Float64}
  y::Array{Float64}
  z::Array{Float64}
  x′::Array{Float64}
  y′::Array{Float64}
  z′::Array{Float64}
  x″::Array{Float64}
  y″::Array{Float64}
  z″::Array{Float64}
  x‴::Array{Float64}
  y‴::Array{Float64}
  z‴::Array{Float64}
end

export solve_for_optimal_path

"""
    solve_for_optimal_path(N, initial, final, constraints, max_iterations, params)

Returns the optimal path for a given number of trajectory segments N, an initial state, final state and box constraints defining safe regions.

# Example
```julia-repl
julia> N = 4
julia> # initial position
julia> X₀ = Dict()
julia> X₀[:x], X₀[:y], X₀[:z] = -5.0, -5.0, -5.0
julia> # initial velocity
julia> X₀′ = Dict()
julia> X₀′[:x], X₀′[:y], X₀′[:z] = 1e0, 1e0, 1e0
julia> # initial acceleration
julia> X₀″ = Dict()
julia> X₀″[:x], X₀″[:y], X₀″[:z] = 1, 1, 1
julia> # final position
julia> X₁ = Dict()
julia> X₁[:x], X₁[:y], X₁[:z] = 22.5, -7.5, 0.0
julia> # final velocity
julia> X₁′ = Dict()
julia> X₁′[:x], X₁′[:y], X₁′[:z] = 1e0, 1e0, 1e0
julia> # final acceleration
julia> X₁″ = Dict()
julia> X₁″[:x], X₁″[:y], X₁″[:z] = 1, 1, 1
julia> initial_state = State(X₀, X₀′, X₀″)
julia> final_state = State(X₁, X₁′, X₁″)
julia> Mxl, Mxu, Myl, Myu, Mzl, Mzu = -10.0, 25.0, -10.0, 20.0, -10.0, 30.0
julia> num_boxes = 1
julia> boxes = Array{Box}(UndefInitializer(), num_boxes)
julia> box₁ = Box(-10.00, 0.00, -10.00, 20.00, -10.00, 10.00)
julia> boxes[1] = box₁
julia> params = OptimizationParameters(Mxl, Mxu, Myl, Myu, Mzl, Mzu)
julia> model, solve_time, optimal_path = solve_for_optimal_path(N, initial_state, final_state, boxes, params=params)
````
"""
function solve_for_optimal_path(N::Int,
    initial::State,
    final::State,
    constraints::Array{Box},
    max_iterations::Int=38400,
    params::BigM_Parameters=nothing)
  boxes = constraints
  X₀ = initial.X
  X₀′ = initial.X′
  X₀″ = initial.X″
  X₁ = final.X
  X₁′ = final.X′
  X₁″ = final.X″

  model = SOSModel(optimizer_with_attributes(Pajarito.Optimizer,
  "conic_solver" => optimizer_with_attributes(Hypatia.Optimizer),
  "oa_solver" => optimizer_with_attributes(Mosek.Optimizer)))
  set_attribute(model, "iteration_limit", max_iterations)

  r = 3 # order of polynomial
  @polyvar(t)
  Z = monomials([t], 0:r)
  @variable(model, H[1:N, boxes], Bin)

  p = Dict()
  # Reformulate disjunctive constraint using big M method,
  # see https://optimization.cbe.cornell.edu/index.php?title=Disjunctive_inequalities#Big-M_Method
  # or https://en.wikipedia.org/wiki/Big_M_method
  if params != nothing
    Mxl = params.Mxl
    Mxu = params.Mxu
    Myl = params.Myl
    Myu = params.Myu
    Mzl = params.Mzl
    Mzu = params.Mzu
  else
    Mxl = minimum(x -> x.xl, constraints)
    Mxu = maximum(x -> x.xu, constraints)
    Myl = minimum(x -> x.yl, constraints)
    Myu = maximum(x -> x.yu, constraints)
    Mzl = minimum(x -> x.zl, constraints)
    Mzu = maximum(x -> x.zu, constraints)
  end

  # setup the polynomials so they satisfy the box constraints and each polynomial can be a member of only a single segment
  T = collect(0:1:N)
  for j in 1:N
    @constraint(model, sum(H[j, box] for box in boxes) == 1)
    p[(:x, j)] = @variable(model, [1:1], SumOfSquares.Poly(Z))
    p[(:y, j)] = @variable(model, [1:1], SumOfSquares.Poly(Z))
    p[(:z, j)] = @variable(model, [1:1], SumOfSquares.Poly(Z))
    S = @set t >= T[j] && t <= T[j+1]
    for box in boxes
      xl, xu, yl, yu, zl, zu = box.xl, box.xu, box.yl, box.yu, box.zl, box.zu
      @constraint(model, p[(:x, j)][1] >= Mxl + (xl-Mxl)*H[j, box], domain = S)
      @constraint(model, p[(:x, j)][1] <= Mxu + (xu-Mxu)*H[j, box], domain = S)
      @constraint(model, p[(:y, j)][1] >= Myl + (yl-Myl)*H[j, box], domain = S)
      @constraint(model, p[(:y, j)][1] <= Myu + (yu-Myu)*H[j, box], domain = S)
      @constraint(model, p[(:z, j)][1] >= Mzl + (zl-Mzl)*H[j, box], domain = S)
      @constraint(model, p[(:z, j)][1] <= Mzu + (zu-Mzu)*H[j, box], domain = S)
    end
  end

  # formulate optimization problem
  for ax in (:x, :y, :z)
    @constraint(model, p[(ax, 1)][1]([0]) == X₀[ax])
    @constraint(model, differentiate(p[ax, 1][1], t, 1)([0]) == X₀′[ax])
    @constraint(model, differentiate(p[ax, 1][1], t, 2)([0]) == X₀″[ax])
    for j in 1:N-1
      @constraint(model, p[(ax, j)][1]([T[j+1]]) == p[(ax, j+1)][1]([T[j+1]]))
      @constraint(model, differentiate(p[(ax, j)][1], t, 1)([T[j+1]]) == differentiate(p[(ax, j+1)][1], t, 1)([T[j+1]]))
    end
    @constraint(model, p[(ax, N)][1]([T[N]]) == X₁[ax])
    @constraint(model, differentiate(p[ax, N][1], t, 1)([T[N]]) == X₁′[ax])
    @constraint(model, differentiate(p[ax, N][1], t, 2)([T[N]]) == X₁″[ax])
  end

  @variable(model, γ[keys(p)] >= 0)
  for (key, val) in p
    # bound the absolute jerk (third derivative of position with respect to time)
    @constraint(model, γ[key] >= differentiate(val[1], t, 3))
    @constraint(model, -γ[key] <= -differentiate(val[1], t, 3))
  end

  @objective(model, Min, sum(γ))
  # solve the problem
  solve_time = @timed optimize!(model)
  if termination_status(model) != OPTIMAL
    return model, nothing, [], [], [], [], [], [], [], []
  end

  x = Array{Float32}([])
  y = Array{Float32}([])
  z = Array{Float32}([])
  x′ = Array{Float32}([])
  y′ = Array{Float32}([])
  z′ = Array{Float32}([])
  x″ = Array{Float32}([])
  y″ = Array{Float32}([])
  z″ = Array{Float32}([])
  x‴ = Array{Float32}([])
  y‴ = Array{Float32}([])
  z‴ = Array{Float32}([])

  # collect values
  for n in 2:N
    for i in T[n-1]:(T[n]-T[n-1])/20:T[n]
            if i == T[n]
                break
            end
      xₜ, yₜ, zₜ = value.(p[(:x, n-1)])[1](i), value.(p[(:y, n-1)])[1](i), value.(p[(:z, n-1)])[1](i)
      push!(x, xₜ)
      push!(y, yₜ)
      push!(z, zₜ)
      xₜ′, yₜ′, zₜ′ = differentiate(value.(p[(:x, n-1)])[1], t, 1)([i]), differentiate(value.(p[(:y, n-1)])[1], t, 1)([i]), differentiate(value.(p[(:z, n-1)])[1], t, 1)([i])
      push!(x′, xₜ′)
      push!(y′, yₜ′)
      push!(z′, zₜ′)
      xₜ″, yₜ″, zₜ″ = differentiate(value.(p[(:x, n-1)])[1], t, 2)([i]), differentiate(value.(p[(:y, n-1)])[1], t, 2)([i]), differentiate(value.(p[(:z, n-1)])[1], t, 2)([i])
      push!(x″, xₜ″)
      push!(y″, yₜ″)
      push!(z″, zₜ″)
      xₜ‴, yₜ‴, zₜ‴ = differentiate(value.(p[(:x, n-1)])[1], t, 3)([i]), differentiate(value.(p[(:y, n-1)])[1], t, 3)([i]), differentiate(value.(p[(:z, n-1)])[1], t, 3)([i])
      push!(x‴, xₜ‴)
      push!(y‴, yₜ‴)
      push!(z‴, zₜ‴)
    end
  end

  optimal_path = OptimalPath(x, y, z, x′, y′, z′, x″, y″, z″, x‴, y‴, z‴)
  return model, solve_time.time, optimal_path
end

end