using Ferrite, FerriteMeshParser, Tensors
using FerriteAssembly, FerriteProblems, FESolvers
using MaterialModelsBase
import FerriteProblems as FP
import MaterialModelsBase as MMB
import FerriteAssembly.ExampleElements: ElasticPlaneStrain, PoroElasticPlaneStrain

elastic_material() = ElasticPlaneStrain(;E=2.e3, ν=0.3)
poroelastic_material() = PoroElasticPlaneStrain(;E=2.e3, ν=0.3, k=0.05, α=1.0, β=1/2e3)

function get_grid()
    # Import grid from abaqus mesh
    grid = get_ferrite_grid(joinpath(@__DIR__, "porous_media", "porous_media_0p75.inp"))

    # Create cellsets for each fieldhandler
    addcellset!(grid, "solid3", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS3")))
    addcellset!(grid, "solid4", intersect(getcellset(grid, "solid"), getcellset(grid, "CPS4R")))
    addcellset!(grid, "porous3", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS3")))
    addcellset!(grid, "porous4", intersect(getcellset(grid, "porous"), getcellset(grid, "CPS4R")))

    # Create faceset for the sides and top
    addfaceset!(grid, "sides", x->(first(x) < eps() || first(x) ≈ 5.0))
    addfaceset!(grid, "top", x->(last(x) ≈ 10.0))
    return grid
end

function create_definition(;t_rise=0.1, p_max=100.0)
    grid = get_grid()

    # Define interpolations
    ipu_quad = Lagrange{RefQuadrilateral,2}()^2
    ipu_tri  = Lagrange{RefTriangle,2}()^2
    ipp_quad = Lagrange{RefQuadrilateral,1}()
    ipp_tri  = Lagrange{RefTriangle,1}()

    # Quadrature rules
    qr_quad = QuadratureRule{RefQuadrilateral}(2)
    qr_tri  = QuadratureRule{RefTriangle}(2)

    # CellValues
    cvu_quad = CellValues(qr_quad, ipu_quad)
    cvu_tri = CellValues(qr_tri, ipu_tri)
    cvp_quad = CellValues(qr_quad, ipp_quad)
    cvp_tri = CellValues(qr_tri, ipp_tri)

    # Setup the DofHandler
    dh = DofHandler(grid)
    # Solid quads
    sdh_solid_quad = SubDofHandler(dh, getcellset(grid,"solid4"))
    add!(sdh_solid_quad, :u, ipu_quad)
    # Solid triangles
    sdh_solid_tri = SubDofHandler(dh, getcellset(grid,"solid3"))
    add!(sdh_solid_tri, :u, ipu_tri)
    # Porous quads
    sdh_porous_quad = SubDofHandler(dh, getcellset(grid, "porous4"))
    add!(sdh_porous_quad, :u, ipu_quad)
    add!(sdh_porous_quad, :p, ipp_quad)
    # Porous triangles
    sdh_porous_tri = SubDofHandler(dh, getcellset(grid, "porous3"))
    add!(sdh_porous_tri, :u, ipu_tri)
    add!(sdh_porous_tri, :p, ipp_tri)

    close!(dh)

    # Setup each domain
    domains = Dict{String,DomainSpec}()
    # Solid domain with Triangle elements, quadratic displacement interpolation
    domains["solid3"] = DomainSpec(sdh_solid_tri, elastic_material(), cvu_tri)

    # Solid domain with Quadrilateral elements, quadratic displacement interpolation
    domains["solid4"] = DomainSpec(sdh_solid_quad, elastic_material(), cvu_quad)

    # Porous domain with Triangle elements
    domains["porous3"] = DomainSpec(sdh_porous_tri, poroelastic_material(), (u=cvu_tri, p=cvp_tri))

    # Porous domain with Quadrilateral elements
    domains["porous4"] = DomainSpec(sdh_porous_quad, poroelastic_material(), (u=cvu_quad, p=cvp_quad))

    # Add boundary conditions
    ch = ConstraintHandler(dh);
    # Fix bottom in y and sides in x
    add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), Returns(0.0), [2]))
    add!(ch, Dirichlet(:u, getfaceset(grid, "sides"), Returns(0.0), [1]))
    # Zero pressure on top surface
    add!(ch, Dirichlet(:p, getfaceset(grid, "top"), Returns(0.0)))
    close!(ch)

    # Add Neumann boundary conditions - normal traction on top
    lh = LoadHandler(dh);
    add!(lh, Neumann(:u, 2, getfaceset(grid, "top"), (x,t,n) -> -n*clamp(t/t_rise,0,1)*p_max))

    return FEDefinition(domains; ch, lh)
end;

struct PM_PostProcess{PVD}
    pvd::PVD
    filestem::String
end
function PM_PostProcess(filestem="porous_media")
    pvd = paraview_collection("$filestem.pvd")
    return PM_PostProcess(pvd, filestem)
end

function FESolvers.postprocess!(post::PM_PostProcess, p, step, solver)
    vtk_grid("$(post.filestem)-$step", FP.get_dofhandler(p)) do vtk
        vtk_point_data(vtk, FP.get_dofhandler(p), FP.getunknowns(p))
        vtk_save(vtk)
        post.pvd[step] = vtk
    end
end

FP.close_postprocessing(post::PM_PostProcess, args...) = vtk_save(post.pvd);

problem = FerriteProblem(create_definition(), PM_PostProcess())
solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper(map(x->x^2, range(0, 1, 41))))
solve_problem!(problem, solver)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
