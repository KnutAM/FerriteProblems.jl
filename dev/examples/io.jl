using Ferrite, FerriteProblems
import FerriteProblems as FP
include("plasticity.jl");

struct CellStress{T}
    σvm::Vector{T}
end
CellStress(grid::Grid) = CellStress(zeros(getncells(grid)))
von_mises(σ) = (s=dev(σ); sqrt((3/2)*s⊡s))

function FerriteAssembly.integrate_cell!(cs::CellStress, state, ae, m::J2Plasticity, cv, buffer)
    σvm = 0.0
    V = 0.0
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        ϵ = function_symmetric_gradient(cv, q_point, ae)
        σ = m.D ⊡ (ϵ - state[q_point].ϵp)
        σvm += von_mises(σ)*dΩ
        V += dΩ
    end
    cs.σvm[cellid(buffer)] = σvm/V
end

plts = FerriteIO("B/FerriteIO.jld2") do io
    def = FP.getdef(io)
    buf = FP.FEBuffer(def)
    post = FP.getpost(io);

    # Then, we get the time history and the displacement data saved to the `post` struct
    t_history = FP.gettimedata(io)
    u_mag = post.umag

    plt1 = plot()
    plot!(plt1, t_history, u_mag*1e3)
    title!(plt1, "maximum displacement")
    xlabel!(plt1, "time [s]")
    ylabel!(plt1, "umax [mm]")

    # The maximum von Mises stress for each step is calculated next. Note that `step` refers
    # to the count of saved steps, and not the actual simulation steps.
    function get_max_vm_stress(cs, buffer, io, step)
        states = FP.getipdata(io, step, "state") # ::Dict{Int,Vector{J2PlasticityMaterialState}}
        FerriteAssembly.update_states!(FerriteAssembly.get_state(buffer), states)
        work!(Integrator(cs), buffer; a=FP.getdofdata(io, step))
        return maximum(cs.σvm)
    end
    dh = FerriteAssembly.get_dofhandler(def)
    cs = CellStress(dh.grid)
    buffer = FP.getassemblybuffer(buf)
    σ_vm_max = Float64[]
    for step in 1:length(t_history)
        push!(σ_vm_max, get_max_vm_stress(cs, buffer, io, step))
    end

    # Plot the analyzed results
    plt2=plot()
    plot!(plt2, t_history, σ_vm_max*1e-6)
    title!(plt2, "Maximum von Mises stress")
    xlabel!(plt2, "time [s]")
    ylabel!(plt2, "stress [MPa]")

    # Finally, cs now contains the von mises stress
    # of each cell at the last step, which we save to vtk.
    u = FP.getdofdata(io, length(t_history))
    vtk_grid("plasticity", dh) do vtkfile
        vtk_point_data(vtkfile, dh, u) # displacement field
        vtk_cell_data(vtkfile, cs.σvm, "von Mises [Pa]")
    end;
    # Return the plots
    (plt1, plt2)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
