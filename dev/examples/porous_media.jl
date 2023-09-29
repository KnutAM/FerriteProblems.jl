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

    # Setup the interpolation and integration rules
    dim=Ferrite.getdim(grid)
    ip3_lin = Lagrange{dim, RefTetrahedron, 1}()
    ip4_lin = Lagrange{dim, RefCube, 1}()
    ip3_quad = Lagrange{dim, RefTetrahedron, 2}()
    ip4_quad = Lagrange{dim, RefCube, 2}()
    qr3 = QuadratureRule{dim, RefTetrahedron}(1)
    qr4 = QuadratureRule{dim, RefCube}(2)

    # Setup the MixedDofHandler
    dh = MixedDofHandler(grid)
    fh1 = FieldHandler([Field(:u, ip3_quad, dim)], getcellset(grid,"solid3"))
    add!(dh, fh1)
    fh2 = FieldHandler([Field(:u, ip4_quad, dim)], getcellset(grid,"solid4"))
    add!(dh, fh2)
    fh3 = FieldHandler([Field(:u, ip3_quad, dim), Field(:p, ip3_lin, 1)], getcellset(grid,"porous3"))
    add!(dh, fh3)
    fh4 = FieldHandler([Field(:u, ip4_quad, dim), Field(:p, ip4_lin, 1)], getcellset(grid,"porous4"))
    add!(dh, fh4)
    close!(dh)

    # Setup each domain
    domains = Dict{String,DomainSpec}()
    # Solid domain with Triangle elements, quadratic displacement interpolation
    sdh1 = FerriteAssembly.SubDofHandler(dh, fh1)
    cv1 = CellVectorValues(qr3, ip3_quad, ip3_lin)
    domains["solid3"] = DomainSpec(sdh1, elastic_material(), cv1)

    # Solid domain with Quadrilateral elements, quadratic displacement interpolation
    sdh2 = FerriteAssembly.SubDofHandler(dh, fh2)
    cv2 = CellVectorValues(qr4, ip4_quad, ip4_lin)
    domains["solid4"] = DomainSpec(sdh2, elastic_material(), cv2)

    # Porous domain with Triangle elements
    # Taylor hood: (quadratic displacement and linear pressure interpolation)
    sdh3 = FerriteAssembly.SubDofHandler(dh, fh3)
    cv3 = (u=CellVectorValues(qr3, ip3_quad, ip3_lin), p=CellScalarValues(qr3, ip3_lin))
    domains["porous3"] = DomainSpec(sdh3, poroelastic_material(), cv3)

    # Porous domain with Quadrilateral elements
    # Taylor hood: (quadratic displacement and linear pressure interpolation)
    sdh4 = FerriteAssembly.SubDofHandler(dh, fh4)
    cv4 = (u=CellVectorValues(qr4, ip4_quad, ip4_lin), p=CellScalarValues(qr4, ip4_lin))
    domains["porous4"] = DomainSpec(sdh4, poroelastic_material(), cv4)

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
