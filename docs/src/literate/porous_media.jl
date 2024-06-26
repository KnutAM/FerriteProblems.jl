# # Porous media 
# 
# Porous media is a two-phase material, consisting of solid parts and a liquid occupying
# the pores inbetween. 
# Using the porous media theory, we can model such a material without explicitly 
# resolving the microstructure, but by considering the interactions between the 
# solid and liquid. In this example, we will additionally consider larger linear 
# elastic solid aggregates that are impermeable. Hence, there is no liquids in 
# these particles and the only unknown variable is the displacement field `:u`. 
# In the porous media, denoted the matrix, we have both the displacement field,
# `:u`, as well as the liquid pressure, `:p`, as unknown. The computational domain
# is shown below (*the outdated figure doesn't show correct boundary conditions*)
#
# ```@raw html
# <table><tbody><tr height="300px"><td style="text-align: left;">
# ```
# ![Computational domain](porous_media/domain.svg)
# ```@raw html
# </td><td>
# ```
# ![Pressure evolution](porous_media/pressure.gif)
# ```@raw html
# </td><td>
# ```
# ![Pressure legend](porous_media/pressure_legend.png)
# ```@raw html
# </td><td>
# ```
# ![u2 evolution](porous_media/u2.gif)
# ```@raw html
# </td><td>
# ```
# ![u2 legend](porous_media/u2_legend.png)
# ```@raw html
# </tr><tr><td>
# ```
# Computational domain
# ```@raw html
# </td><td>
# ```
# Pressure evolution
# ```@raw html
# </td><td></td><td>
# ```
# Vertical displacements¨
# ```@raw html
# </td><td>
# ```
# ```@raw html
# </td></tr></tbody></table>
# ```
#
#
# ## Theory of porous media
# The strong forms are given as
# ```math
# \begin{aligned}
# \boldsymbol{\sigma}(\boldsymbol{\epsilon}, p) \cdot \boldsymbol{\nabla} &= \boldsymbol{0} \\
# \dot{\Phi}(\boldsymbol{\epsilon}, p) + \boldsymbol{w}(p) \cdot \boldsymbol{\nabla} &= 0
# \end{aligned}
# ```
# where 
# ``\boldsymbol{\epsilon} = \left[\boldsymbol{u}\otimes\boldsymbol{\nabla}\right]^\mathrm{sym}`` 
# The constitutive relationships are 
# ```math
# \begin{aligned}
# \boldsymbol{\sigma} &= \boldsymbol{\mathsf{E}}:\boldsymbol{\epsilon} - \alpha p \boldsymbol{I} \\
# \boldsymbol{w} &= - k \boldsymbol{\nabla} p \\
# \Phi &= \phi + \alpha \mathrm{tr}(\boldsymbol{\epsilon}) + \beta p
# \end{aligned}
# ``` 
# with 
# ``\boldsymbol{\mathsf{E}}=2G \boldsymbol{\mathsf{I}}^\mathrm{dev} + 3K \boldsymbol{I}\otimes\boldsymbol{I}``.
# The material parameters are then the 
# shear modulus, ``G``, 
# bulk modulus, ``K``, 
# permeability, ``k``,  
# Biot's coefficient, ``\alpha``, and
# liquid compressibility, ``\beta``.
# The porosity, ``\phi``, doesn't enter into the equations 
# (A different porosity leads to different skeleton stiffness and permeability).
#
# 
# The variational (weak) form can then be derived for the variations ``\boldsymbol{\delta u}``
# and ``\delta p`` as
# ```math
# \begin{aligned}
# \int_\Omega \left[\left[\boldsymbol{\delta u}\otimes\boldsymbol{\nabla}\right]^\mathrm{sym}:
# \boldsymbol{\mathsf{E}}:\boldsymbol{\epsilon} - \boldsymbol{\delta u} \cdot \boldsymbol{\nabla} \alpha p\right] \mathrm{d}\Omega 
# &= \int_\Gamma \boldsymbol{\delta u} \cdot \boldsymbol{t} \mathrm{d} \Gamma \\
# \int_\Omega \left[\delta p \left[\alpha \dot{\boldsymbol{u}} \cdot \boldsymbol{\nabla} + \beta \dot{p}\right] + 
# \boldsymbol{\nabla}(\delta p) \cdot [k \boldsymbol{\nabla}]\right] \mathrm{d}\Omega 
# &= \int_\Gamma \delta p w_\mathrm{n} \mathrm{d} \Gamma 
# \end{aligned}
# ```
# where ``\boldsymbol{t}=\boldsymbol{n}\cdot\boldsymbol{\sigma}`` is the traction and 
# ``w_\mathrm{n} = \boldsymbol{n}\cdot\boldsymbol{w}`` is the normal flux.  
# 
# ### Finite element form
# Discretizing in space using finite elements, we obtain the vector equation 
# ``r_i = f_i^\mathrm{int} - f_{i}^\mathrm{ext}`` where ``f^\mathrm{ext}`` are the external 
# "forces", and ``f_i^\mathrm{int}`` are the internal "forces". We split this into the 
# displacement part ``r_i^\mathrm{u} = f_i^\mathrm{int,u} - f_{i}^\mathrm{ext,u}`` and 
# pressure part ``r_i^\mathrm{p} = f_i^\mathrm{int,p} - f_{i}^\mathrm{ext,p}``
# to obtain the discretized equation system 
# ```math
# \begin{aligned}
# f_i^\mathrm{int,u} &= \int_\Omega [\boldsymbol{\delta N}^\mathrm{u}_i\otimes\boldsymbol{\nabla}]^\mathrm{sym} : \boldsymbol{\mathsf{E}} : [\boldsymbol{u}\otimes\boldsymbol{\nabla}]^\mathrm{sym} \ 
# - [\boldsymbol{\delta N}^\mathrm{u}_i \cdot \boldsymbol{\nabla}] \alpha p \mathrm{d}\Omega 
# &= \int_\Gamma \boldsymbol{\delta N}^\mathrm{u}_i \cdot \boldsymbol{t} \mathrm{d} \Gamma \\
# f_i^\mathrm{int,p} &= \int_\Omega \delta N_i^\mathrm{p} [\alpha [\dot{\boldsymbol{u}}\cdot\boldsymbol{\nabla}]  + \beta\dot{p}] + \boldsymbol{\nabla}(\delta N_i^\mathrm{p}) \cdot [k \boldsymbol{\nabla}(p)] \mathrm{d}\Omega 
# &= \int_\Gamma \delta N_i^\mathrm{p} w_\mathrm{n} \mathrm{d} \Gamma
# \end{aligned}
# ```
# Approximating the time-derivatives, ``\dot{\boldsymbol{u}}\approx \left[\boldsymbol{u}-{}^n\boldsymbol{u}\right]/\Delta t``
# and ``\dot{p}\approx \left[p-{}^np\right]/\Delta t``, we can implement the finite element equations in the residual form 
# ``r_i(\boldsymbol{a}(t), t) = 0`` where the vector ``\boldsymbol{a}`` contains all unknown displacements ``u_i`` and pressures ``p_i``. 
# We use automatic differentiation to get the jacobian.  

# ## Implementation
# We now solve the problem step by step. The full program with fewer comments is found in 
#md # the final [section](@ref porous-media-plain-program)
# 
# Required packages
using Ferrite, FerriteMeshParser, Tensors
using FerriteAssembly, FerriteProblems, FESolvers
using MaterialModelsBase
import FerriteProblems as FP
import MaterialModelsBase as MMB
import FerriteAssembly.ExampleElements: ElasticPlaneStrain, PoroElasticPlaneStrain
# ## Physics
# Both the elastic material and a poroelastic material is available in
# `FerriteAssembly.ExampleElements`:
elastic_material() = ElasticPlaneStrain(;E=2.e3, ν=0.3)
poroelastic_material() = PoroElasticPlaneStrain(;E=2.e3, ν=0.3, k=0.05, α=1.0, β=1/2e3)

# ## Problem definition
# ### Mesh import
# In this example, we import the mesh from the Abaqus input file, 
# [`porous_media_0p75.inp`](porous_media/porous_media_0p75.inp) using `FerriteMeshParser`'s 
# `get_ferrite_grid` function. 
# (A finer mesh, [`porous_media_0p25.inp`](porous_media/porous_media_0p25.inp), is also available)
# We then create one cellset for each phase (solid and porous)
# for each element type. These 4 sets will later be used in their own `FieldHandler`
function get_grid()
    ## Import grid from abaqus mesh
    grid = get_ferrite_grid(joinpath(@__DIR__, "porous_media", "porous_media_0p75.inp"))

    ## Create cellsets for each fieldhandler
    addcellset!(grid, "solid3", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS3")))
    addcellset!(grid, "solid4", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS4R")))
    addcellset!(grid, "porous3", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS3")))
    addcellset!(grid, "porous4", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS4R")))

    ## Create faceset for the sides and top 
    addfaceset!(grid, "sides", x->(first(x) < eps() || first(x) ≈ 5.0))
    addfaceset!(grid, "top", x->(last(x) ≈ 10.0))
    return grid
end

# ### Problem setup 
# Define the finite element interpolation, integration, and boundary conditions. 
function create_definition(;t_rise=0.1, p_max=100.0)
    grid = get_grid()

    ## Define interpolations
    ipu_quad = Lagrange{RefQuadrilateral,2}()^2
    ipu_tri  = Lagrange{RefTriangle,2}()^2
    ipp_quad = Lagrange{RefQuadrilateral,1}()
    ipp_tri  = Lagrange{RefTriangle,1}()

    ## Quadrature rules
    qr_quad = QuadratureRule{RefQuadrilateral}(2)
    qr_tri  = QuadratureRule{RefTriangle}(2)

    ## CellValues
    cvu_quad = CellValues(qr_quad, ipu_quad)
    cvu_tri = CellValues(qr_tri, ipu_tri)
    cvp_quad = CellValues(qr_quad, ipp_quad)
    cvp_tri = CellValues(qr_tri, ipp_tri)

    ## Setup the DofHandler
    dh = DofHandler(grid)
    ## Solid quads
    sdh_solid_quad = SubDofHandler(dh, getcellset(grid,"solid4"))
    add!(sdh_solid_quad, :u, ipu_quad)
    ## Solid triangles
    sdh_solid_tri = SubDofHandler(dh, getcellset(grid,"solid3"))
    add!(sdh_solid_tri, :u, ipu_tri)
    ## Porous quads
    sdh_porous_quad = SubDofHandler(dh, getcellset(grid, "porous4"))
    add!(sdh_porous_quad, :u, ipu_quad)
    add!(sdh_porous_quad, :p, ipp_quad)
    ## Porous triangles
    sdh_porous_tri = SubDofHandler(dh, getcellset(grid, "porous3"))
    add!(sdh_porous_tri, :u, ipu_tri)
    add!(sdh_porous_tri, :p, ipp_tri)

    close!(dh)

    ## Setup each domain
    domains = Dict{String,DomainSpec}()
    ## Solid domain with Triangle elements, quadratic displacement interpolation
    domains["solid3"] = DomainSpec(sdh_solid_tri, elastic_material(), cvu_tri)
    
    ## Solid domain with Quadrilateral elements, quadratic displacement interpolation
    domains["solid4"] = DomainSpec(sdh_solid_quad, elastic_material(), cvu_quad)

    ## Porous domain with Triangle elements
    domains["porous3"] = DomainSpec(sdh_porous_tri, poroelastic_material(), (u=cvu_tri, p=cvp_tri))

    ## Porous domain with Quadrilateral elements
    domains["porous4"] = DomainSpec(sdh_porous_quad, poroelastic_material(), (u=cvu_quad, p=cvp_quad))

    ## Add boundary conditions
    ch = ConstraintHandler(dh);
    ## Fix bottom in y and sides in x
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), Returns(0.0), [2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "sides"), Returns(0.0), [1]))
    ## Zero pressure on top surface
    add!(ch, Dirichlet(:p, getfaceset(grid, "top"), Returns(0.0)))
    close!(ch)

    ## Add Neumann boundary conditions - normal traction on top
    lh = LoadHandler(dh);
    add!(lh, Neumann(:u, 2, getfaceset(grid, "top"), (x,t,n) -> -n*clamp(t/t_rise,0,1)*p_max))

    return FEDefinition(domains; ch, lh)
end;

# ## Postprocessing
struct PM_PostProcess{PVD}
    pvd::PVD
    filestem::String
end
function PM_PostProcess(filestem="porous_media")
    pvd = paraview_collection("$filestem.pvd")
    return PM_PostProcess(pvd, filestem)
end

function FESolvers.postprocess!(post::PM_PostProcess, p, solver)
    step = FESolvers.get_step(solver)
    vtk_grid("$(post.filestem)-$step", FP.get_dofhandler(p)) do vtk
        vtk_point_data(vtk, FP.get_dofhandler(p), FP.getunknowns(p))
        vtk_save(vtk)
        post.pvd[step] = vtk 
    end
end

FP.close_postprocessing(post::PM_PostProcess, args...) = vtk_save(post.pvd);

# ## Solving
# We solve the problem by using linearly increasing time steps
problem = FerriteProblem(create_definition(), PM_PostProcess())
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(map(x->x^2, range(0, 1, 41))))
solve_problem!(problem, solver)


#md # ## [Plain program](@id porous-media-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`porous_media.jl`](porous_media.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```