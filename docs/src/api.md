```@meta
CurrentModule = FerriteProblems
```

# API
The main types below are exported. 
Remaining functions are not exported to avoid polluting the 
name space. *Tip:* To simplify calling the following functions
it is possible to write `import FerriteProblems as FP` as is done
in the examples. 

## Main types
```@docs
FerriteProblem
FEDefinition
FerriteIO
```

## Setting up simulation
```@docs
FerriteProblems.cellbuffertype
FerriteProblems.allocate_material_cache
```

## Convergence criteria
In normal usage, the following convergence criteria can be used
```@docs
FerriteProblems.ConvergenceCriterion
FerriteProblems.AbsoluteResidual
FerriteProblems.RelativeResidualElementScaling
```

To create custom convergence criteria, the following functions 
may require overloading. 
```@docs
FerriteProblems.TolScaling
FerriteProblems.make_assemscaling
FESolvers.calculate_convergence_measure
```

## Access functions
```@docs
FerriteProblems.getdh
FerriteProblems.getch
FerriteProblems.getnh
FerriteProblems.getcv
FerriteProblems.getmaterial
FerriteProblems.getbodyload
FerriteProblems.getjacobian
FerriteProblems.getunknowns
FerriteProblems.getresidual
FerriteProblems.getneumannforce
FerriteProblems.getoldunknowns
FerriteProblems.getstate
FerriteProblems.getoldstate
FerriteProblems.gettime
FerriteProblems.getoldtime
```

## Saving and loading data
```@docs
FESolvers.postprocess!
FerriteProblems.close_postprocessing
FerriteProblems.addstep!
FerriteProblems.gettimedata
FerriteProblems.savedofdata!
FerriteProblems.getdofdata
FerriteProblems.savenodedata!
FerriteProblems.getnodedata
FerriteProblems.savecelldata!
FerriteProblems.getcelldata
FerriteProblems.saveipdata!
FerriteProblems.getipdata
FerriteProblems.saveglobaldata!
FerriteProblems.getglobaldata
FerriteProblems.getdef
FerriteProblems.getpost
```