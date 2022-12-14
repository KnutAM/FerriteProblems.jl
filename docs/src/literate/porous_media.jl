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
# is shown below
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
using FerriteAssembly, FerriteProblems, FerriteNeumann, FESolvers
using MaterialModelsBase
import FerriteProblems as FP
import MaterialModelsBase as MMB

# ## Physics
# ### Elastic material 
# For the elastic material, we just define the material following 
# the `MaterialModelsBase.jl` interface:
struct Elastic{T} <: AbstractMaterial
    G::T
    K::T
    E4::SymmetricTensor{4,2,T,9}
end
function Elastic(;E=2.e3, ??=0.3)
    G = E / 2(1 + ??)
    K = E / 3(1 - 2??)
    I2 = one(SymmetricTensor{2,2})
    I4vol = I2???I2
    I4dev = minorsymmetric(otimesu(I2,I2)) - I4vol / 3
    E4 = 2G*I4dev + K*I4vol
    return Elastic(G, K, E4)
end;

function MMB.material_response(m::Elastic, ??, args...; kwargs...)
    ?? = m.E4 ??? ??
    return ??, m.E4, NoMaterialState()
end;

# ### Poroelastic material
# For the poroelastic material, we reuse the elastic part from above, 
# but add additionally required properties 
struct PoroElastic{E<:Elastic,T}
    elastic::E
    k::T    # [mm^4/Ns] Permeability
    ??::T    # [-] Biot's coefficient
    ??::T    # [1/MPa] Liquid bulk modulus
end 
function PoroElastic(;elastic=Elastic(), k=0.05, ??=1.0, ??=1/2e3)
    return PoroElastic(elastic, k, ??, ??)
end;

# And then we also have to define the specific element routine 
# which we do by defining the residual and using autodiff to calculate the tangent
function FerriteAssembly.element_residual!(re, state, ae, material::PoroElastic, cv::NamedTuple, dh_fh, ??t, buffer)
    ## Setup cellvalues and give easier names
    cv_u = cv[:u]
    cv_p = cv[:p]
    num_u = getnbasefunctions(cv_u)
    num_p = getnbasefunctions(cv_p)

    ## Assign views to the matrix and vector parts
    ae_old = FerriteAssembly.get_aeold(buffer)
    udofs = dof_range(dh_fh, :u)
    pdofs = dof_range(dh_fh, :p)
    ru = @view re[udofs]
    rp = @view re[pdofs]
    ue = @view ae[udofs]
    pe = @view ae[pdofs]
    ue_old = @view ae_old[udofs]
    pe_old = @view ae_old[pdofs]

    ## Assemble stiffness and force vectors
    for q_point in 1:getnquadpoints(cv_u)   
        ## Calculate variables in the current quadrature point 
        d?? = getdetJdV(cv_u, q_point)
        ?? = function_symmetric_gradient(cv_u, q_point, ue)
        ??_old = function_symmetric_gradient(cv_u, q_point, ue_old)
        p = function_value(cv_p, q_point, pe)
        p_old = function_value(cv_p, q_point, pe_old)
        ???p = function_gradient(cv_p, q_point, pe)
        pdot = (p-p_old)/??t 
        div_udot = (tr(??)-tr(??_old))/??t
        ??eff = material.elastic.E4 ??? ??

        ## Assemble residual contributions
        for i??? in 1:num_u
            ?????Nu = shape_symmetric_gradient(cv_u, q_point, i???)
            div_??Nu = shape_divergence(cv_u, q_point, i???)
            ru[i???] += (?????Nu ??? ??eff - div_??Nu*material.??*p)*d??
        end
        for i??? in 1:num_p
            ??Np = shape_value(cv_p, q_point, i???)
            ?????Np = shape_gradient(cv_p, q_point, i???)
            rp[i???] += (??Np*(material.??*div_udot + material.??*pdot) + (?????Np ??? ???p)*material.k) * d??
        end
    end
end

# ## Problem definition
# ### Mesh import
# In this example, we import the mesh from the Abaqus input file, 
# [`porous_media_0p75.inp`](porous_media_0p75.inp) using `FerriteMeshParser`'s 
# `get_ferrite_grid` function. 
# (A finer mesh, [`porous_media_0p25.inp`](porous_media_0p25.inp), is also available)
# We then create one cellset for each phase (solid and porous)
# for each element type. These 4 sets will later be used in their own `FieldHandler`
function get_grid()
    ## Import grid from abaqus mesh
    grid = get_ferrite_grid(joinpath(@__DIR__, "porous_media_0p75.inp"))

    ## Create cellsets for each fieldhandler
    addcellset!(grid, "solid3", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS3")))
    addcellset!(grid, "solid4", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS4R")))
    addcellset!(grid, "porous3", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS3")))
    addcellset!(grid, "porous4", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS4R")))

    ## Create faceset for the sides and top 
    addfaceset!(grid, "sides", x->(first(x) < eps() || first(x) ??? 5.0))
    addfaceset!(grid, "top", x->(last(x) ??? 10.0))
    return grid
end

# ### Problem setup 
# Define the finite element interpolation, integration, and boundary conditions. 
function create_definition(;t_rise=0.1, p_max=100.0)

    grid = get_grid()

    ## Setup the interpolation and integration rules
    dim=Ferrite.getdim(grid)    
    ip3_lin = Lagrange{dim, RefTetrahedron, 1}()
    ip4_lin = Lagrange{dim, RefCube, 1}()
    ip3_quad = Lagrange{dim, RefTetrahedron, 2}()
    ip4_quad = Lagrange{dim, RefCube, 2}()
    qr3 = QuadratureRule{dim, RefTetrahedron}(1)
    qr4 = QuadratureRule{dim, RefCube}(2)

    ## Setup the MixedDofHandler
    dh = MixedDofHandler(grid)
    push!(dh, FieldHandler([Field(:u, ip3_lin, dim)], getcellset(grid,"solid3")))
    push!(dh, FieldHandler([Field(:u, ip4_lin, dim)], getcellset(grid,"solid4")))
    push!(dh, FieldHandler([Field(:u, ip3_quad, dim), Field(:p, ip3_lin, 1)], getcellset(grid,"porous3")))
    push!(dh, FieldHandler([Field(:u, ip4_quad, dim), Field(:p, ip4_lin, 1)], getcellset(grid,"porous4")))
    close!(dh)

    ## Setup cellvalues with the same order as the FieldHandlers in the dh
    ## - Linear displacement elements in the solid domain
    ## - Taylor hood (quadratic displacement, linear pressure) and linear geometry in porous domain
    cv = ( CellVectorValues(qr3, ip3_lin), 
           CellVectorValues(qr4, ip4_lin),
           (u=CellVectorValues(qr3, ip3_quad, ip3_lin), p=CellScalarValues(qr3, ip3_lin)),
           (u=CellVectorValues(qr4, ip4_quad, ip4_lin), p=CellScalarValues(qr4, ip4_lin)) )

    ## Add boundary conditions
    ## Use `Ferrite.jl` PR427 (temporarily included in FerriteProblems.jl) 
    ## to make Dirichlet conditions easier and more general
    ch = ConstraintHandler(dh);
    ## Fix bottom in y and sides in x
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> zero(Vec{1}), [2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "sides"), (x,t) -> zero(Vec{1}), [1]))
    ## Zero pressure on top surface
    add!(ch, Dirichlet(:p, getfaceset(grid, "top"), (x,t) -> 0.0))
    close!(ch)

    ## Add Neumann boundary conditions - normal traction on top
    nh = NeumannHandler(dh);
    add!(nh, Neumann(:u, 2, getfaceset(grid, "top"), (x,t,n) -> -n*clamp(t/t_rise,0,1)*p_max))

    ## We then need one material per fieldhandler:
    materials = (Elastic(), Elastic(), PoroElastic(), PoroElastic())

    return FEDefinition(;dh=dh, ch=ch, nh=nh, cv=cv, m=materials, cc=FP.RelativeResidualElementScaling())
end

# ## Postprocessing
struct PostProcess{PVD}
    pvd::PVD
    filestem::String
end
function PostProcess(filestem="porous_media")
    pvd = paraview_collection("$filestem.pvd")
    return PostProcess(pvd, filestem)
end

function FESolvers.postprocess!(post::PostProcess, p, step, solver)
    vtk_grid("$(post.filestem)-$step", FP.getdh(p)) do vtk
        vtk_point_data(vtk, FP.getdh(p), FP.getunknowns(p))
        vtk_save(vtk)
        post.pvd[step] = vtk 
    end
end

FP.close_postprocessing(post::PostProcess, args...) = vtk_save(post.pvd)

# ## Solving
# Without the real PR427 there will be a lot of warnings. when creating the definition.
# As a temporary solution, we just disable warnings during definition creation
using Logging
Logging.disable_logging(Logging.Warn)
problem = FerriteProblem(create_definition(), PostProcess())
Logging.disable_logging(Logging.LogLevel(-1))
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