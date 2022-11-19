```@meta
CurrentModule = FerriteProblems
```

# API
## Main types and running simulations
```@docs
FerriteProblem
FEDefinition
FerriteIO
safesolve
FP.close_problem
```

## Access functions
```@docs
FP.getdh
FP.getch
FP.getnh
FP.getcv
FP.getmaterial
FP.getbodyload
FP.getjacobian
FP.getunknowns
FP.getresidual
FP.getneumannforce
FP.getoldunknowns
FP.getstate
FP.getoldstate
FP.gettime
FP.getoldtime
```

## Special cases
```@docs
FP.allocate_material_cache
```

## Saving and loading data
```@docs
FESolvers.postprocess!
FP.addstep!
FP.gettimedata
FP.savedofdata!
FP.getdofdata
FP.savenodedata!
FP.getnodedata
FP.savecelldata!
FP.getcelldata
FP.saveipdata!
FP.getipdata
FP.saveglobaldata!
FP.getglobaldata
FP.getdef
FP.getpost
```