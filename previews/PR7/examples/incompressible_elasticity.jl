using Ferrite
using FerriteNeumann, FerriteAssembly, FerriteProblems
using FESolvers
import FerriteProblems as FP
import FerriteAssembly as FA

function create_cook_grid(nx, ny)
    corners = [Vec{2}((0.0,   0.0)),
               Vec{2}((48.0, 44.0)),
               Vec{2}((48.0, 60.0)),
               Vec{2}((0.0,  44.0))]
    grid = generate_grid(Triangle, (nx, ny), corners);
    # facesets for boundary conditions
    addfaceset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0);
    addfaceset!(grid, "traction", x -> norm(x[1]) ≈ 48.0);
    return grid
end;

function create_values(interpolation_u, interpolation_p)
    # quadrature rules
    qr      = QuadratureRule{2,RefTetrahedron}(3)
    face_qr = QuadratureRule{1,RefTetrahedron}(3)

    # geometric interpolation
    interpolation_geom = Lagrange{2,RefTetrahedron,1}()

    # cell and facevalues for u
    cellvalues_u = CellVectorValues(qr, interpolation_u, interpolation_geom)
    facevalues_u = FaceVectorValues(face_qr, interpolation_u, interpolation_geom)

    # cellvalues for p
    cellvalues_p = CellScalarValues(qr, interpolation_p, interpolation_geom)

    return cellvalues_u, cellvalues_p, facevalues_u
end;

function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    push!(dh, :u, 2, ipu) # displacement
    push!(dh, :p, 1, ipp) # pressure
    close!(dh)
    return dh
end;

struct LinearElasticity{T}
    G::T
    K::T
end
function LinearElasticity(;Emod, ν)
    Gmod = Emod / 2(1 + ν)
    Kmod = Emod * ν / ((1+ν) * (1-2ν))
    return LinearElasticity(Gmod, Kmod)
end

function FerriteAssembly.allocate_cell_cache(::LinearElasticity, cv::NamedTuple)
    cellvalues_u = cv[:u]
    return collect([symmetric(shape_gradient(cellvalues_u, 1, i)) for i in 1:getnbasefunctions(cellvalues_u)])
end

function create_definition(ν, ip_u, ip_p)
    grid = create_cook_grid(50, 50)
    dh = create_dofhandler(grid, ip_u, ip_p)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfaceset(dh.grid, "clamped"), Returns(zero(Vec{2}))))
    close!(ch)
    update!(ch, 0.0)

    cv_u, cv_p, fv = create_values(ip_u, ip_p)
    cv = (u=cv_u, p=cv_p)   # Create NamedTuple

    nh = NeumannHandler(dh)
    add!(nh, Neumann(:u, fv, getfaceset(grid, "traction"), (x,t,n)->Vec{2}((0.0, 1/16))))

    m = LinearElasticity(;Emod=1.0, ν=ν)

    # Create and return the `FEDefinition`
    return FEDefinition(;dh=dh, ch=ch, nh=nh, cellvalues=cv, material=m)
end;

function FerriteAssembly.element_routine!(
    Ke, re, state, ue, mp::LinearElasticity, cv::NamedTuple, buffer
    )
    cellvalues_u = cv[:u]
    cellvalues_p = cv[:p]

    # Get the local indices for each field
    udofs = dof_range(buffer, :u)
    pdofs = dof_range(buffer, :p)

    # Extract cached gradients
    ∇Nu_sym_dev = FA.get_cache(buffer)

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    for q_point in 1:getnquadpoints(cellvalues_u)
        for i in 1:length(udofs)
            ∇Nu_sym_dev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
        end
        dΩ = getdetJdV(cellvalues_u, q_point)
        for (i_u, i) in enumerate(udofs)
            for (j_u, j) in enumerate(udofs[1:i_u])
                Ke[i,j] += 2 * mp.G * ∇Nu_sym_dev[i_u] ⊡ ∇Nu_sym_dev[j_u] * dΩ
            end
        end

        for (i_p, i) in enumerate(pdofs)
            δNp = shape_value(cellvalues_p, q_point, i_p)
            for (j_u, j) in enumerate(udofs)
                divδNu = shape_divergence(cellvalues_u, q_point, j_u)
                Ke[i,j] += -δNp * divδNu * dΩ
            end
            for (j_p, j) in enumerate(pdofs[1:i_p])
                Np = shape_value(cellvalues_p, q_point, j_p)
                Ke[i,j] += - 1/mp.K * δNp * Np * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

struct PostProcessing
    vtk_file::String
end

function FESolvers.postprocess!(post::PostProcessing, p, step, solver)
    step == 1 && return nothing # We don't want to save the initial conditions.
    dh = FP.getdh(p)
    vtk_grid(post.vtk_file, dh) do vtkfile
        vtk_point_data(vtkfile, dh, FP.getunknowns(p))
    end
end;

function build_problem(ν, ip_u, ip_p)
    def = create_definition(ν, ip_u, ip_p)
    ip_u_string = isa(ip_u, Lagrange{2,RefTetrahedron,1}) ? "linear" : "quadratic"
    post = PostProcessing("cook_$(ip_u_string)_linear")
    return FerriteProblem(def, post)
end

solver = QuasiStaticSolver(;nlsolver=LinearProblemSolver(), timestepper=FixedTimeStepper([0.0,1.0]))
ν = 0.4999999
linear    = Lagrange{2,RefTetrahedron,1}()
quadratic = Lagrange{2,RefTetrahedron,2}()
p1 = build_problem(ν, linear, linear)
p2 = build_problem(ν, quadratic, linear)
solve_problem!(p1, solver)
solve_problem!(p2, solver)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
