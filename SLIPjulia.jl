using JuMP, Ipopt
using Plots

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 4)

# STepping
N = 800
@variable(model, dt >= 5e-4)  # Timestep
@constraint(model, dt <= 0.05)

# State,control variables
@variable(model, x[1:5, 1:(N+1)]) # [x z dx dz theta]
@variable(model, u[1:N])         
@variable(model, ks[1:N])          
@variable(model, th0, start = pi/6) 
@constraint(model, pi/12 <= th0 <= pi/3)

# Model parameters
m = 1.5   # mass (kg)
l0 = 1.0   # rest length (m)
g = 9.81   # gravity (m/s^2)
v_t = 1.0 # target velocity

# Initial targets
@constraint(model, x[1,1] == 0)
@constraint(model, x[2,1] == l0*cos(th0))
@constraint(model, x[3,1] == v_t)
@constraint(model, x[4,1] == 0)
@constraint(model, x[5,1] == th0)

# Final targets (equal to initial to set periodic motion)
@constraint(model, x[5,N+1] == -th0)
@constraint(model, x[3,N+1] == v_t)
@constraint(model, x[4,N+1] == 0)
@constraint(model, x[2,N+1] == l0*cos(th0))

#Mid-cycle constraint
z = Int(N/2)
@constraint(model, x[5,z] == 0)
px = l0*sin(th0)

for k in 1:N
    # Current state
    xc = x[1,k]
    zc = x[2,k]
    xv = x[3,k]
    zv = x[4,k]
    th1 = x[5,k]
    
    # Leg length and spring forces
    lx = px - xc
    ly = zc
    Fs_lx = -1*(l0*sin(th1) - lx)*ks[k]
    Fs_ly = (l0*cos(th1) - ly)*ks[k]
    
    # Accelerations
    ax = Fs_lx/m
    az = (Fs_ly - m*g)/m
    
    # State update
    @constraint(model, x[1,k+1] == x[1,k] + dt*x[3,k])
    @constraint(model, x[2,k+1] == x[2,k] + dt*x[4,k])
    @constraint(model, x[3,k+1] == x[3,k] + dt*ax)
    @constraint(model,x[4,k+1] == x[4,k] + dt*az)
    @constraint(model,x[5,k+1] == x[5,k] + u[k])
    
    #Constraints
    @constraint(model, -1.5 <= th1 <= 1.5)
    @constraint(model, 0.65 <= x[2,k] <= l0)
    @constraint(model, u[k] <= 0)
    @constraint(model, 5 <= ks[k] <= 200)

end
@objective(model, Min,((x[1,N+1] - x[1,1])/(dt*N) - v_t)^2)
@objective(model, Min, sum((u[k])^2 for k in 1:N))
optimize!(model)

# Plot
if termination_status(model) == MOI.OPTIMAL
    time = [value(dt)*i for i in 0:N]
    th_vals = [value(u[i]) for i in 1:N]
    x_vals = [value(x[1,i]) for i in 1:N+1]
    y_vals = [value(x[2,i]) for i in 1:N+1]

    T = N*value(dt)
    println("Total time: ", T)
    
    p1 = plot(time, x_vals, label="x position")
    p2 = plot(time, y_vals, label="y position")
    plot(p1, p2, layout=(2,1))
end

time = [value(dt)*i for i in 1:N+1]
u_vals = [value(u[i]) for i in 1:N]
x_vals = [value(x[1,i]) for i in 1:N+1]
y_vals = [value(x[2,i]) for i in 1:N+1]
th_vals = [value(x[5,i]) for i in 1:N+1]

T = N.*value(dt)
k = [value(ks[i]) for i in 1:N]
print(value(dt), ' ',T, ' ', value(th0))


b = 0
bf = 5
anim = @animate for i in 1:30:bf*(N+1) 
    global b
    if i > (b+1)*(N+1)  # Check if we've passed into next cycle
        b = b + 1
    end
    
    it = i - b*(N+1)    # Calculate frame within current cycle
    if it > N+1         # Handle edge case
        it = N+1
    end

    plot(xlim=(-0.375, 2.875), ylim=(-0.05, 1.25), 
    aspect_ratio=:equal, 
    legend=false, 
    grid=true)

    dist = x_vals[N+1]
    x_pos = x_vals[it] + dist*b
    y_pos = y_vals[it]
    theta = th_vals[it]
    
    plot!([-0.5, 3], [0, 0], color=:black, linewidth=3)
    
    leg_base_x = l0*sin(value(th0)) + dist*b # ground contact point
    plot!([leg_base_x, x_pos], [0, y_pos], color=:blue, linewidth=5)
    leg2_x = l0*sin(-1*theta) + x_pos
    plot!([leg2_x, x_pos], [0, y_pos], color=:green, linewidth=4, linestyle=:dash)
    scatter!([x_pos], [y_pos], color=:red, markersize=20)

    title!("Time: $(round(time[it], digits=3)) s")
end
gif(anim, "slip_walker.gif", fps=120)

# Save animation

# ramp up and down
# angle of swing leg as control variable
# lots of optimization courses in IOE.
