```@meta
CurrentModule = FerriteProblems
```

# FerriteProblems

Documentation for [FerriteProblems](https://github.com/KnutAM/FerriteProblems.jl).

When using the `FESolvers.jl` package together with `Ferrite.jl`, 
the user has to specify a `problem` to be solved. 
The purpose of `FESolvers.jl` is to keep this as general as possible, 
and not require users to structure their `problem` structs in a certain way. 

This package brings together the 
[`FESolvers.jl`](https://github.com/KnutAM/FESolvers.jl)
package with 
[`Ferrite.jl`](https://github.com/Ferrite-FEM/Ferrite.jl),
as well as the supporting packages 
[`FerriteAssembly.jl`](https://github.com/KnutAM/FerriteAssembly.jl) and 
[`FerriteNeumann.jl`](https://github.com/KnutAM/FerriteNeumann.jl). 
There is also preliminary support for organizing your simulations 
by saving both setup and results as `.jld2` using [`JLD2.jl`](https://github.com/JuliaIO/JLD2.jl)

All unregistered dependencies, including `FerriteProblems.jl` itself, 
is available in the [`knutamregistry`](https://github.com/KnutAM/knutamregistry)

## Simple workflow
The easiest way to get started is to just follow the examples, but in brief the workflow
for a simple setup is

1. `def = FEDefinition(;dh=dh, ch=ch, material=material, cellvalues=cellvalues)`
2. `problem = FerriteProblem(def)`
3. `solver = QuasiStaticSolver(NewtonSolver(), FixedTimeStepper(;num_steps=10,Δt=0.1))`
4. `solve_problem!(problem, solver)`

Note that the documentation of `Ferrite.jl`, `FESolvers.jl`, 
`FerriteAssembly.jl`, and `FerriteNeumann.jl` should be considered as well. 