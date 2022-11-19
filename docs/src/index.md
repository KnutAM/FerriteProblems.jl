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
There is also extensive support for organizing your simulations 
by saving both setup and results as `.jld2` using [`JLD2.jl`](https://github.com/JuliaIO/JLD2.jl)

The workflow requires you to create its `FerriteProblem` type which 
can then be solved by `FESolvers.jl`. 
The easiest way to get started is to just follow the examples. 
Note that the documentation of `Ferrite.jl`, `FerriteAssembly.jl`, 
and `FerriteNeumann.jl` should be considered as well. 

```julia
struct FerriteProblem{POST,DEF,BUF,IOT}
    def::DEF
    post::POST
    buf::BUF
    io::IOT
end
```
with `FerriteProblem(def, post=nothing, io=nothing)` as constructor. 
The four parts have distinct tasks:
* `def::FEDefinition` is responsible for the full problem definition. 
  I.e., given `def`, the full simulation should be possible to replicate, 
  given the same solver from `FESolvers.jl`
* `post` contains all information related to the postprocessing of each step. 
  This typically varies a lot between simulations, 
  and is the first type parameter to allow easy 
  dispatch on `problem`s with different `post`s.  
* `buf::FEBuffer` contains all buffer values, 
  these are not necessary (nor desirable) to save, 
  and can be recreated each time the constructor is called.
* `io::FerriteIO`: This field enable file handling to allow easy saving 
  and retrieving of results from a simulation using JLD2 files.
