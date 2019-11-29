"""
    Trajectory of a dipole in an electric field
    ===========================================
"""

function voltages_t(voltages::Vector{Float64}, t::Float64)
    """ electrode voltages at time t
    """
    return voltages
end

function voltages_t(voltages::Function, t::Float64)
    """ electrode voltages at time t
    """
    return voltages(t)
end

function kinetic_energy(v::Vector{Float64}, mass::Float64)
    """ kinetic energy of a moving particle
    """
    return 0.5 * mass * (v[1]^2.0 + v[2]^2.0 + v[3]^2.0)
end

function potential_energy(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, xg::Float64, yg::Float64, zg::Float64, dipole::Float64; get_famp::Bool=false)
    """ potential energy of a static dipole
    """
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)::Float64
    if get_famp
        return dipole * famp, famp
    else
        return dipole * famp
    end
end

function potential_energy(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, xg::Float64, yg::Float64, zg::Float64, dipole::Function; get_famp::Bool=false)
    """ potential energy of a field-dependent dipole
    """
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)::Float64
    dpm = dipole(famp)::Float64
    if get_famp
        return dpm * famp, famp
    else
        return dpm * famp
    end
end

function force(fa::FastAdjust3D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole_f::Function)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (3D)
    """
    xg, yg, zg = grid_r(fa, r)
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)::Float64
    dipole = dipole_f(famp)::Float64
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx::Float64
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy::Float64
    fz = -(potential_energy(fa, voltages, t, xg, yg, zg + 0.5, dipole) - 
           potential_energy(fa, voltages, t, xg, yg, zg - 0.5, dipole)) / fa.dz::Float64
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust2D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole_f::Function)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (2D)
    """
    xg, yg, zg = grid_r(fa, r)
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)::Float64
    dipole = dipole_f(famp)::Float64
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx::Float64
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy::Float64
    fz = 0.0::Float64
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust3D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Float64)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (3D)
    """
    xg, yg, zg = grid_r(fa, r)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx::Float64
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy::Float64
    fz = -(potential_energy(fa, voltages, t, xg, yg, zg + 0.5, dipole) - 
           potential_energy(fa, voltages, t, xg, yg, zg - 0.5, dipole)) / fa.dz::Float64
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust2D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Float64)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (2D)
    """
    xg, yg, zg = grid_r(fa, r)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx::Float64
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy::Float64
    fz = 0.0::Float64
    return Vector([fx, fy, fz])
end

function acceleration(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64)
    """ accelleration of a dipole in an inhomegeous electric field
    """
    return  force(fa, voltages, t, r, dipole) ./ mass
end

function euler(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ Euler method of integration
    """
    a = acceleration(fa, voltages, t, r, dipole, mass)
    rn = r .+ v .* dt
    vn = v .+ a .* dt
    return rn, vn, a
end

function leapfrog(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ leapfrog method of integration
    """
    if any(isnan.(a))
        a = acceleration(fa, voltages, t, r, dipole, mass)
    end
    rn = r .+ v .* dt + a .* (0.5 * dt^2.0)
    an = acceleration(fa, voltages, t + dt, rn, dipole, mass)
    vn = v .+ (a .+ an) .* (0.5 * dt)
    return rn, vn, an
end

function rk4(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ fourth order Runge-Kutta method of integration
    """
    # initial position
    a = acceleration(fa, voltages, t, r, dipole, mass)
    # half-step
    r1 = r .+ (0.5 * dt) .* v
    v1 = v .+ (0.5 * dt) .* a
    a1 = acceleration(fa, voltages, t + 0.5 * dt, r1, dipole, mass)
    # half-step again
    r2 = r .+ (0.5 * dt) .* v1
    v2 = v .+ (0.5 * dt) .* a1
    a2 = acceleration(fa, voltages, t + 0.5 * dt, r2, dipole, mass)
    # full step
    r3 = r .+ dt .* v2
    v3 = v .+ dt .* a2
    a3 = acceleration(fa, voltages, t + dt, r3, dipole, mass)
    # next position
    r4 = r .+ (dt / 6.0) .* (v .+ 2.0 .* v1 .+ 2.0 .* v2 .+ v3)
    v4 = v .+ (dt / 6.0) .* (a .+ 2.0 .* a1 .+ 2.0 .* a2 .+ a3)
    return r4, v4, a3
end

function ion_classic_MaxLFS(n::Int64)
    #((1.286 * 1e9) / (9.0 * (n^4.0))) * 100.0
    return (((2.0 * 1.286 * 1e9) / (9.0 * (n^4.0))) * 100.0)
end

function ion_rates(n::Int64, k::Int64, m::Int64, fldamp::Float64)
    n2 = (-k+n-1.0-m)/2.0 ::Float64
    E_nm = (-1.0*E_hPs / (2.0* n^2.0)) + E_Stark(n,k,m,fldamp) ::Float64
    R = (-2.0*E_nm)^1.5 / (q*a_Ps*sqrt(E_hPs)*fldamp) ::Float64
    A = exp(-(2.0*R/3.0) - (0.25*(n^3.0 *q*a_Ps*fldamp)/(E_hPs))*(34.0*n2^2. + 34.0*n2*m + 46.0*n2 + 7.0*m^2.0
                    + 23.0*m + 53.0/3.0)) ::Float64
    return (A * ((4.0*R).^(2.0*n2+m+1.0)) * E_hPs) / (hbar* n^3.0 * (factorial(convert(Int8, n2))) * (factorial(convert(Int8, n2+m))))
end

function ion_none()
end

function E_Stark(n::Int64, k::Int64, m::Int64, fldamp::Float64) #CONVERT
    first_term = (3.0/2.0) * n * k * q * a_Ps * fldamp ::Float64
    second_term = -(1.0/16.0) * n^4.0  * (17.0 * n^2.0 - 3.0 * k^2.0 -9.0 *m^2.0 +19.0)* q^2.0 * a_Ps^2.0 * fldamp^2.0 / E_hPs  ::Float64 
	third_term = (3.0/32.0) * n^7.0 *k*(23.0* n^2.0 - k^2.0 + 11.0*m^2.0 + 39.0)* q^3.0 * a_Ps^3.0 * fldamp^3.0 / (E_hPs^2.0)  ::Float64
    fourth_term = -(1.0/1024.0) * n^10.0 * (5487.0* n^4.0 + 35182.0* n^2.0 -1134.0*(m^2.0)*k^2.0 + 1806.0*(n^2)*k^2.0 - 3402.0* (m^2.0)*n^2.0 + 147.0*k^4.0 -549.0*m^4.0 + 5754.0*k^2.0 - 8622.0*m^2.0 +16211.0)* q^4.0 * a_Ps^4.0 * fldamp^4.0 / (E_hPs^3.0)::Float64
    
    Stark_Energy = first_term + second_term + third_term + fourth_term ::Float64
    return Stark_Energy
end



function traj(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
              t0::Float64, r0::Vector{Float64}, v0::Vector{Float64},
              dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function,
			  n::Int64, k::Int64, m::Int64, 
			  ionisation::Function, killit::Function;
              max_field::Float64=NaN, max_iterations::Int64=Int64(1e6), df::Bool=false,
			  t_flo::Float64=NaN, t_quench::Float64=NaN, t_anni::Float64=NaN)
    """ full trajectory
    """
    # initialise
    i = 1::Int64
    t, r, v = t0, r0, v0
	quenched = false
    a = Vector([NaN, NaN, NaN])
    xg, yg, zg = grid_r(fa, r)
    KE = kinetic_energy(v, mass)::Float64
	#dipole = -(3.0 / 2.0) * state[1] * state[2] * q * a_Ps
    PE = potential_energy(fa, voltages, t, xg, yg, zg, dipole)::Float64
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)::Float64
    result = [t r[1] r[2] r[3] KE PE famp;]
	
	if ~isnan(t_flo)
		t_flo_kill = (-log(rand()) * t_flo)::Float64
	end
	if ~isnan(t_anni)
		t_anni_kill = (-log(rand()) * t_anni)::Float64
	end

    # step-by-step trajectory
	death = "Max Iteration"::String
    while i < max_iterations
        try
            r, v, a = method(fa, voltages, t, r, v, a, dipole, mass, dt)
            t += dt::Float64
        catch
			death = "Tried and Failed"::String
            break
        end
		
		#Check for invalid value
        if any(isnan, r) | any(isnan, v)
			death = "isnan"::String
            break
        end
        # grid position
        xg, yg, zg = grid_r(fa, r)
        # energy
        KE = kinetic_energy(v, mass)::Float64
        PE, famp = potential_energy(fa, voltages, t, xg, yg, zg, dipole, get_famp=true)::Tuple{Float64,Float64}

        # Deaths
        if electrode_g(fa, xg, yg, zg)
            # hit an electrode or left the pa
			death = "Splat"::String
            break
        elseif lost_g(fa, xg, yg, zg)
            # left the pa?
			death = "Lost"::String
            break
			
        elseif ~isnan(max_field) && famp > max_field
            # ionisation field
			death = "Ionised (Classic)"::String
            break
		elseif ionisation==ion_rates && rand() > (1.0 - exp(ion_rates(n,k,m,famp)*dt*-1.0))
            # ionisation field
			death = "Ionised (Rate)"::String
            break
			
		elseif ~isnan(t_flo) && (t-t0) > t_flo_kill
			death = "Fluoresced"::String
			#Fluoresced
			break
		elseif ~isnan(t_anni) && (t-t0) > t_anni_kill
			death = "Self Annihilated"::String
			#Fluoresced
			break
			
		elseif killit(t,r,v)
			death = "Kill It"::String
			break
		end
		
		# Test for quenching
		if ~isnan(t_quench) && t > t_quench && quenched == false
			t_anni = 142e-9::Float64
			t_anni_kill = (-log(rand()) * t_anni)::Float64
			quenched = true::Bool
		end
		
        # record
        result = vcat(result, [t r[1] r[2] r[3] KE PE famp])
        # next step
        i += 1
    end
	result = vcat(result, [death n k m 0.0 0.0 0.0])
	
    # output
    if df
        result = DataFrame(result)::DataFrame
        names!(result, [:t, :x, :y, :z, :KE, :PE, :famp]);
    end
    return result
end




function traj_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                   t0::Float64, r0::Vector{Float64}, v0::Vector{Float64},
                   dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                   max_field::Float64=NaN, max_t::Float64=NaN, max_iterations::Int64=Int64(1e6))
    """ last point of a trajectory
    """
    # initialise
    i = 1
    t, r, v = t0, r0, v0
    a = Vector([NaN, NaN, NaN])
    result = [t r[1] r[2] r[3] v[1] v[2] v[3];]
    # step-by-step trajectory
    while i < max_iterations
        try
            r, v, a = method(fa, voltages, t, r, v, a, dipole, mass, dt)
            t += dt
        catch
            break
        end
        if any(isnan, r) | any(isnan, v)
            break
        end
        # grid position
        xg, yg, zg = grid_r(fa, r)
        # checks
        if electrode_g(fa, xg, yg, zg)
            # hit an electrode or left the pa
            break
        elseif ~isnan(max_field) && amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg) > max_field
            # ionisation field
            break
        elseif ~isnan(max_t) && t > max_t
            # time's up!
            break
        end
        # record
        result = [t r[1] r[2] r[3] v[1] v[2] v[3];]
        # next step
        i += 1
    end
    # output
    return result
end
