using Ferrite, FerriteProblems
import FerriteProblems as FP
include("plasticity.jl");

io = FerriteIO("B/FerriteIO.jld2")
def = FP.getdef(io)
post = FP.getpost(io);

t_history = FP.gettimedata(io)
u_mag = post.umag

plt1 = plot()
plot!(plt1, t_history, u_mag*1e3)
title!(plt1, "maximum displacement")
xlabel!(plt1, "time [s]")
ylabel!(plt1, "umax [mm]")

display(plt1) #!src

function get_max_vm_stress(step)
    states = FP.getipdata(io, step, "state")
    σ_vm = maximum(cellstates -> maximum(state -> vonMises(state.σ), cellstates), states)
    return σ_vm
end

σ_vm = get_max_vm_stress.(1:length(t_history));

plt2=plot()
plot!(plt2, t_history, σ_vm*1e-6)
title!(plt2, "Maximum von Mises stress")
xlabel!(plt2, "time [s]")
ylabel!(plt2, "stress [MPa]")

display(plt2) #!src

step = length(t_history)
u = FP.getdofdata(io, step)
states = FP.getipdata(io, step, "state")
dh = FP.getdh(def)
mises_values = zeros(getncells(dh.grid))
for (el, cell_states) in enumerate(states)
    for state in cell_states
        mises_values[el] += vonMises(state.σ)
    end
    mises_values[el] /= length(cell_states) # average von Mises stress
end
vtk_grid("plasticity", dh) do vtkfile
    vtk_point_data(vtkfile, dh, u) # displacement field
    vtk_cell_data(vtkfile, mises_values, "von Mises [Pa]")
end;

close(io)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

