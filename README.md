# FerriteProblems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KnutAM.github.io/FerriteProblems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteProblems.jl/dev/)
[![Build Status](https://github.com/KnutAM/FerriteProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteProblems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteProblems.jl)

When using the `FESolvers.jl` package together with `Ferrite.jl`, the user has to specify 
a `problem` to be solved. The purpose of `FESolvers.jl` is to keep this as general as possible, 
and not require users to structure their `problem` structs in a certain way. 

This package brings together the `FESolvers.jl` package with `Ferrite.jl`, as well as 
the supporting packages `FerriteAssembly.jl` and `FerriteNeumann.jl`. 

**Warning:** *This package is experimental and breaking changes are expected*
