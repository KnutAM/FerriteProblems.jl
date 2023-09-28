using Ferrite, FerriteProblems
import FerriteProblems as FP
include("plasticity.jl");

plts = FerriteIO("B/FerriteIO.jld2") do io
    def = FP.getdef(io)
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
    function get_max_vm_stress(io, step)
        states = FP.getipdata(io, step, "state") # ::Dict{Int,Vector{J2PlasticityMaterialState}}
        σ_vm = maximum(cellstates -> maximum(state -> vonMises(state.σ), cellstates), values(states))
        return σ_vm
    end

    σ_vm = get_max_vm_stress.((io,), 1:length(t_history));

    # Plot the analyzed results
    plt2=plot()
    plot!(plt2, t_history, σ_vm*1e-6)
    title!(plt2, "Maximum von Mises stress")
    xlabel!(plt2, "time [s]")
    ylabel!(plt2, "stress [MPa]")

    # Finally, we extract the degrees of freedom and the stresses in each element
    # at the last step to save to vtk.
    step = length(t_history)
    u = FP.getdofdata(io, step)
    states = FP.getipdata(io, step, "state")
    dh = FP.get_dofhandler(def)
    mises_values = zeros(getncells(dh.grid))
    for (el, cell_states) in pairs(states) # states::Dict{Int}
        for state in cell_states
            mises_values[el] += vonMises(state.σ)
        end
        mises_values[el] /= length(cell_states) # average von Mises stress
    end
    vtk_grid("plasticity", dh) do vtkfile
        vtk_point_data(vtkfile, dh, u) # displacement field
        vtk_cell_data(vtkfile, mises_values, "von Mises [Pa]")
    end;
    # Return the plots
    (plt1, plt2)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
