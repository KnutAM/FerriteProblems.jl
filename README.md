# FerriteProblems

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteProblems.jl/dev/)
[![Build Status](https://github.com/KnutAM/FerriteProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteProblems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteProblems.jl)

**Warning:** *This package is experimental and breaking changes are expected*

When using the `FESolvers.jl` package together with `Ferrite.jl`, the user has to specify 
a `problem` to be solved. The purpose of `FESolvers.jl` is to keep this as general as possible, 
and not require users to structure their `problem` structs in a certain way. 

The `FerriteProblems.jl` package brings together 
* [`FESolvers.jl`](https://github.com/KnutAM/FESolvers.jl)
* [`Ferrite.jl`](https://github.com/Ferrite-FEM/Ferrite.jl)
* [`FerriteAssembly.jl`](https://github.com/KnutAM/FerriteAssembly.jl)
* [`FerriteNeumann.jl`](https://github.com/KnutAM/FerriteNeumann.jl)


## Installation
Those packages not in the general julia registry are available in [knutamregistry](https://github.com/KnutAM/knutamregistry)
After adding this registry, the `FerriteProblems.jl` can be installed as 
```julia
using Pkg
Pkg.add("FerriteProblems")
```

If the registry is not added, it is possible to install using 
```julia
using Pkg
Pkg.add(url="https://github.com/KnutAM/FerriteProblems.jl")
```
But then it is first required to add all of the above packages 
(that are not registered) manually
