"""
    Fly dipoles
    ===========

    flydi
        trajectories for a range of starting conditions
    
    flydi_last 
        final position for a range of starting conditions

"""
function k_uniform(n::Int64,m::Int64)
    k_max=n-1-m::Int64
    k_ar = collect(-k_max:2:k_max)

    idx = 0::Int64
    while idx == 0 || idx==n - m + 1
        idx = round(Int64,  rand()*(length(k_ar) + 1.0) )
    end

    k = k_ar[idx]::Int64
    return k
end


function k_laser(n::Int64, m::Int64, k_center::Int64)

    return k
end

function flydi(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, di::DataFrame,
               n::Int64, k::Union{Function, Int64}, m::Int64, dt::Float64, method::Function, killit::Function;
               max_iterations::Int64=Int64(1e6), max_atom::Int64=Int64(0),
			   t_flo::Float64=NaN, t_quench::Float64=NaN, t_anni::Float64=NaN, mass::Float64=m_Ps,
			   ionisation::Function=ion_none)
    """ Dipole trajectory
    """
    # fly dipoles
	
    num, cols = size(di)::Tuple{Int64,Int64}
	if max_atom==0
		result = Array{DataFrame, 1}(undef, num)
		its = num
	else
		result = Array{DataFrame, 1}(undef, max_atom)
		its = max_atom
	end
	
	if max_atom > num
		return "error 1009xF4: max_atom exceeds length of distribution DataFrame"
	end
	
	if ionisation == ion_classic_MaxLFS
		max_field = ion_classic_MaxLFS(n)::Float64
	else
		max_field = NaN
	end

    @showprogress for i in 1:its
	
		#SELECT K VALUE
		
		if k == k_uniform
			k_pick = k_uniform(n,m)::Int64
		# elseif k == k_laser
			# k_pick = k_laser(n,m)
		else
			k_pick = k::Int64
		end
		
		dipole = -(3.0 / 2.0) * n * k_pick * q * a_Ps::Float64
	
        t0, r0, v0 = di[i, :t], Vector(di[i, [:x, :y, :z]]), Vector(di[i, [:vx, :vy, :vz]])
        result[i] = traj(fa, voltages, t0, r0, v0, dipole, mass, dt, method, n, k_pick, m, ionisation, killit;
                         max_field=max_field, max_iterations=max_iterations, df=true, t_flo=t_flo, t_quench=t_quench, t_anni=t_anni)
						 
		if max_atom!=0 && i > max_atom
			break
		end
		
    end
	#result = DataFrame(result)::DataFrame
    return result
end

function flydi_Para(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, di::DataFrame,
               n::Int64, k::Union{Function, Int64}, m::Int64, mass::Float64, dt::Float64, method::Function;
               max_t::Float64=1.0, max_iterations::Int64=Int64(1e6), max_x::Float64=NaN, max_y::Float64=NaN, max_z::Float64=NaN, max_atom::Int64=Int64(0),
			   t_flo::Float64=NaN, t_quench::Float64=NaN, t_anni::Float64=NaN,
			   ionisation::Function=ion_none)
    """ Dipole trajectory
    """
    # fly dipoles
    num, cols = size(di)
	if max_atom==0
		result = Array{DataFrame, 1}(undef, num)
		its = num
	else
		result = Array{DataFrame, 1}(undef, max_atom)
		its = max_atom
	end
	
	if max_atom > num
		return "error 1009xF4: max_atom exceeds length of distribution DataFrame"
	end
	
	if ionisation == ion_classic_MaxLFS
		max_field = ion_classic_MaxLFS(n)
	else
		max_field = NaN
	end

    for i in 1:its
	
		#SELECT K VALUE
		dipole = -(3.0 / 2.0) * n * k * q * a_Ps
	
        t0, r0, v0 = di[i, :t], Vector(di[i, [:x, :y, :z]]), Vector(di[i, [:vx, :vy, :vz]])
        @async result[i] = traj(fa, voltages, t0, r0, v0, dipole, mass, dt, method, n,k,m, ionisation;
                         max_field=max_field, max_t=max_t, max_iterations=max_iterations, df=true, t_flo=t_flo, t_quench=t_quench, t_anni=t_anni,
						 max_x=max_x, max_y=max_y, max_z=max_z)
						 
		if max_atom!=0 && i > max_atom
			break
		end
		
    end
    return result
end

function flydi_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                     di::DataFrame,
                     state::Vector{Int64}, mass::Float64, dt::Float64, method::Function, ionisation::Function;
                     max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ Dipole trajectory (final position)
    """
    # fly dipoles
    t_initial = Array(di[:t])
    r_initial = convert(Matrix, di[[:x, :y, :z]])
    v_initial = convert(Matrix, di[[:vx, :vy, :vz]])
    num, dims = size(r_initial)
    result = Array{Float64, 2}(undef, num, 7)
	dipole = -(3.0 / 2.0) * state[1] * state[2] * q * a_Ps
    for i in 1:num
        t0, r0, v0 = t_initial[i], r_initial[i, :], v_initial[i, :]
        result[i, :] = traj_last(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                                 max_field=max_field, max_t=max_t, max_iterations=max_iterations)
    end
    result = DataFrame(result)
    names!(result, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return result
end
