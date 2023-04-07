# -----------------------------------------------------------
# Global Surface Temperature Change
# -----------------------------------------------------------

@defcomp temperature_withadj begin

    #d   = Parameter(index=[3])      # Thermal response timescales: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (years).
    decay_factor = Parameter(index=[3]) # Thermal response decay factor, calculated as exp(-1/d) where d represents thermal response timescale of jth thermal box.
    q       = Parameter(index=[3])      # Raditive forcing coefficient: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (K W⁻¹m²).
    F       = Parameter(index=[time])   # Total radiative forcing (Wm⁻²).
    Tj_0    = Parameter(index=[3])
    T_0     = Parameter()

    T1_adjustment = Parameter(index=[time], unit="K")

    T   = Variable(index=[time])    # Global mean surface temperature anomaly (K).
    Tj  = Variable(index=[time,3])  # Temperature change for three thermal pools (K).

    function run_timestep(p, v, d, t)

        if is_first(t)

            # Set initial condition for three thermal boxes.
            v.Tj[t,:] = p.Tj_0

            # Set initial temperature.
            v.T[t] = p.T_0

        else

            #Calculate temperature change for the three different thermal response times.
            for j=1:3
                v.Tj[t,j] = v.Tj[t-1,j] * p.decay_factor[j] + p.F[t] * p.q[j] * (1.0 - p.decay_factor[j])
            end

            # Add adjustment to fastest box
            v.Tj[t,1] = v.Tj[t,1] + p.T1_adjustment[t]
            # Don't need to recalculate, since no interaction between boxes

            #Calculate global mean surface temperature anomaly.
            v.T[t] = sum([v.Tj[t-1,:] v.Tj[t,:]]) / 2
        end
    end
end
