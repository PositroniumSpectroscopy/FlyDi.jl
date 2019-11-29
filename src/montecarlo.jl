"""
    Monte-Carlo distributions
    =========================
"""

function mc_vel(num::Int64; v_s=1.95e3, delta_v=80.0, rho=2.5e-4, z=0.15, epsilon=1.0)
    """ A random sample of velocities with a mean longitudional velocity of v_s
        and a spread of delta_v.  The transverse velocity is found by assuming
        a point source and that the beam passes through a hole at z with a radius of rho
        and with x/y eccentricity of epsilon.        
    """
    # longitudional beam velocity (Gaussian)
    vz = abs.(v_s .+ randn(num) .* (delta_v / (2.0^0.5)))
    # position at the laser
    r = rho .* rand(num).^0.5
    theta = 2.0 * pi * rand(num)
    # transverse velocity (determined by skimmer)
    tof = z ./ vz
    vx = epsilon * r .* cos.(theta) ./ tof
    vy = r .* sin.(theta) ./ tof
    return Array([vx vy vz])
end

function mc_flat(num, xmin=0.0, xmax=5.0e-6)
    """ Random values, evenly sampled from xmin < x < xmax.        
    """
    x = xmin .+ rand(num) .* (xmax - xmin)
    return x
end

function mc_uniform(num::Int64; tmin=0.0, tmax=5e-6, distance=0.3, z0=0.198, reverse_vz=false,
                                v_s=1.86e3, delta_v=65.0, rho=1.5e-4, epsilon=1.0)
    """ 3D He distribution with a flat, oval distribution in xy.
    """
    t =  mc_flat(num, tmin, tmax)
    vel = mc_vel(num, v_s=v_s, delta_v=delta_v, rho=rho, z=distance, epsilon=epsilon)
    tof = distance ./ vel[:, 3]
    x = tof .* vel[:, 1]
    y = tof .* vel[:, 2]
    z = zeros(num) .+ z0
    if reverse_vz
    vel[:, 3] = - vel[:, 3]
    end
    df = DataFrame(hcat(t, x, y, z, vel))
    names!(df, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return df
end

function mc_supersonic(num::Int64; tmin=0.0, tmax=5e-6, distance=0.2, z0=0., reverse_vz=false,
                                   v_s=2.0e3, delta_v=65.0, sigma_x=5e-4, sigma_y=1e-4, rho=1e-4)
    """ A random sample of velocities with a mean longitudional velocity of v_s
        and a spread of delta_v.  The transverse velocity is found by assuming
        a small source (rho) located at z=(z0 - distance) and that the beam has a
        2D Gaussian distribution in xy. 
    """
    # time
    t =  mc_flat(num, tmin, tmax)
    # longitudional velocity
    vz = abs.(v_s .+ randn(num) .* (delta_v / (2.0^0.5)))
    tof = distance ./ vz
    x = randn(num) .* sigma_x
    y = randn(num) .* sigma_y
    # position at the source
    r = rho .* rand(num).^0.5
    theta = 2.0 * pi * rand(num)
    x0 = r .* sin.(theta)
    y0 = r .* cos.(theta)
    # vel
    vx = (x .- x0) ./ tof
    vy = (y .- y0) ./ tof
    z = zeros(num) .+ z0
    if reverse_vz
        vz = - vz
    end
    df = DataFrame(hcat(t, x, y, z, vx, vy, vz))
    names!(df, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return df
end

function Ps_formation(N::Int64, T::Float64;    
                      sigma_t::Float64=2.0e-9, sigma_xy::Float64=2.0e-3, e_Ps_eff::Float64=0.3,
                      x0::Float64=0.0,y0::Float64=0.0,z0::Float64=0.0, t0::Float64=0.0)

    #ADD IN CONVERT POSITRONS TO Ps
    
    #Assign initial time and position data
    t0 = rand(Normal(t0,sigma_t), N)
    x0 = rand(Normal(x0,sigma_xy), N)
    y0 = rand(Normal(y0,sigma_xy), N)
    z0 = fill(z0, N)

    #Efficiency of e+ to o-Ps conversion
#     direct = df[np.random.random(n_positrons) > eff].index
#     df.loc[direct, 'status'] = 'direct'
#     df.loc[direct, 'life'] = 0.0
#     df.loc[direct, ['vx', 'vy', 'vz']] = 0.0
#     #Designate o-Ps
#     oPs = df[df.status != 'direct'].index
#     df.loc[oPs, 'status'] = 'oPs'

    #Velocity distribution
    sigma_v = sqrt((kB * T)/(2.0 * m_e))
    speed_sq = (rand(Normal(0.0,sigma_v), N)).^2 +
    (rand(Normal(0.0,sigma_v), N)).^2 +
    (rand(Normal(0.0,sigma_v), N)).^2
    speed = sqrt.(speed_sq)
    
    # cosine dist. from Greenwood, J. (2002) Vacuum 67 217
    phi = rand(Float64, N) * 2.0 * pi
    
    #theta = np.arcsin(np.sqrt(np.random.random(n_Ps)))
    theta = acos.(rand(Float64, N).^0.5)
    
    # angular velocity distribution assignment
    angles_x = sin.(theta) .* sin.(phi)
    angles_y = sin.(theta) .* cos.(phi)
    angles_z = cos.(theta)
    vel_x = (speed .* angles_x)
    vel_y = (speed .* angles_y)
    vel_z = (speed .* angles_z)
    
    result = Array([t0 x0 y0 z0 vel_x vel_y vel_z])
    df = DataFrame(result)
    names!(df, [:t, :x, :y, :z, :vx, :vy, :vz]);
    
    return df
end


function Ps_Laser_Slice(df::DataFrame;    
                lambda_0::Float64=243.2e-9, linewidth::Float64=1.0/(2*pi*3.2e-9),
                selection_threshold::Float64=0.001, ls_energy::Float64=0.001,
                ls_wavelength::Float64=243.2e-9, ls_bandwidth::Float64=1e11,distance::Float64=0.0005,
                ls_height::Float64=0.01, ls_width::Float64=0.003, trigger::Float64=6e-9,
                ls_sigma_t::Float64=2e-9)
 
    doppler = lambda_0 .* (1.0 .- df.vx ./ c)
    #dopp_vals = doppler.values, dtype='float')
    
    spectral_hits = lineshape(doppler, ls_sigma_wl(ls_bandwidth_wl(ls_bandwidth=ls_bandwidth, ls_wavelength=ls_wavelength))) * ls_delta_lambda(linewidth=linewidth, lambda_0=lambda_0)

    # spatial overlap
    spatial_hits = fluence(df)
    #println(length(spatial_hits)," ", length(spectral_hits))
    #println(spectral_hits)
    hits = df[spectral_hits .* spatial_hits .> selection_threshold,:]
    
    
    return hits
end

function ls_area(;ls_height::Float64=0.01, ls_width::Float64=0.003)
    return ls_height * ls_width
end

function ls_peak_power(;energy::Float64=0.001, signma_t::Float64=2e-9)
    return energy * (2.0 * pi * sigma_t^2.0)^-0.5
end

function ls_peak_intensity(peak_power::Float64, area::Float64)
    return  peak_power / area
end

function ls_bandwidth_wl(;ls_bandwidth::Float64=1e11, ls_wavelength::Float64=243.2e-9)
    return ls_bandwidth * ls_wavelength ^2.0 / c
end

function ls_sigma_wl(ls_bandwidth_wl::Float64)
    return ls_bandwidth_wl/ (2.0 * sqrt(2.0 * log(2)))
end

function ls_power(peak_power::Float64, time::Float64; trigger::Float64=6e-9, ls_sigma_t::Float64=2e-9)
    return peak_power * exp(- (time - trigger)^2.0 / (2.0 * ls_sigma_t^2))
end

function ls_intensity(ls_power::Float64, area::Float64)
    return ls_power / area
end

function ls_fluence(time_1::Float64, time_2::Float64, area::Float64; trigger::Float64=6e-9, ls_sigma_t::Float64=2e-9, ls_energy::Float64=0.001)
    return 0.5*ls_energy/area*(erf((time_2-trigger)/(sqrt(2.0)*ls_sigma_t))-erf((time_1-trigger)/(sqrt(2.0)*ls_sigma_t)))
end

function lineshape(wav::Vector{Float64}, sigma_wl::Float64; ls_wavelength::Float64=243.2e-9)
    return ((2.0 .* pi .* sigma_wl.^2.0).^-0.5) .* exp.(.-(wav.-ls_wavelength).^2.0./(2.0 .* sigma_wl.^2.0))
end

function ls_delta_lambda(;linewidth::Float64=1.0/(2*pi*3.2e-9), lambda_0::Float64=243.2e-9)
""" natural linewidth of the transition"""
    return linewidth * lambda_0^2.0 / c
end

function fluence(df::DataFrame; ls_distance::Float64=0.0005, ls_width::Float64=0.003, ls_height::Float64=0.01)
    """ find the laser fluence experienced by each Ps atom in df
    """
    # which Ps reach the laser in their lifetime and don't go too high or too low?
    #hits = df[(abs.(df.y .+ (ls_distance .+ ls_width ./ 2.0) .* df.vy ./ df.vz) .< ls_height ./ 2.0),:]
    hits = Array{Float64}(undef,size(df,1))
    i = 1::Int64
    @showprogress for rw in eachrow(df)
        if abs(rw.y + (ls_distance + ls_width / 2.0) * rw.vy / rw.vz) < (ls_height / 2.0)
            time_1 = rw.t + (ls_distance / rw.vz) ::Float64
            time_2 = rw.t + ((ls_distance + ls_width) / rw.vz) ::Float64
            hits[i] = ls_fluence(time_1, time_2, ls_area(ls_height=ls_height, ls_width=ls_width))
        else
            hits[i] = 0.0
        end
        i+=1::Int64
    end
    return hits
end
