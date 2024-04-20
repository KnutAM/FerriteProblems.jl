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
FESolvers.calculate_convergence_measure
FerriteProblems.make_assemscaling
```

## Access functions
```@docs
FerriteProblems.get_dofhandler
FerriteProblems.get_constrainthandler
FerriteProblems.get_loadhandler
FerriteProblems.get_material
FerriteProblems.getjacobian
FerriteProblems.getunknowns
FerriteProblems.getresidual
FerriteProblems.get_external_force
FerriteProblems.get_old_unknowns
FerriteProblems.get_state
FerriteProblems.get_old_state
FerriteProblems.get_time
FerriteProblems.get_old_time
```

## Saving and loading data
```@docs
FESolvers.handle_notconverged!(::Any, ::FerriteProblem, ::Any)
FESolvers.postprocess!(::Any, ::FerriteProblem, ::Any)
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