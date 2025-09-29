module VisualizeNbodyGradient
        # """
        #     This module allows user to create figures a la rebound from NbodyGradient, which uses masses in solar masses 𝑀⊙, time in days, distance in AU, and angles in radians.
        #     The right-handed coordinate system is set up so that the positive z-axis is pointing away from the observer, 
        #     and the positive x-axis points to the right along the horizontal. 

        #     In NbodyGradient the orientation is such that the x-y plane is the sky plane.
        #     For inclinations close to 90 degrees this means that the x-z plane represents the top-down view.
        #     Therefore we choose to take the z coordinates for plotting.
        # """

    using NbodyGradient,PyPlot,PyCall
    @pyimport matplotlib.animation as anim
    """
    #Example for Kepler-16:
        stara=Elements(m=0.6897, P=0.0  ,    
            t0=0.0   ,
            ecosω=0.0  ,  
            esinω=0.0  , 
            I=0.0   ,
            Ω=0.0)
        starb=Elements(m=0.20255,
        P=41.079220,I=pi/2,e=0.15944)#,t0=212.12316)
        pl=Elements(m=0.333*0.00095,P=229.0,I=pi/2,e=0.0069 ) # circumbinary planet 
        t0=0.0; tmax=1000.0;h=5.0
        ic=ElementsIC(t0,[-1 1 0;-1 -1 1; -1 -1 -1],[stara;starb;pl])
        intr=Integrator(h,t0,tmax)
        nsteps=100 #number of time steps that we want optional
        kep=Keplerians(intr,ic,nsteps) 
        make_plot(kep)
        animate_plot(kep,save=true,filename="Kepler-16.mp4") 
    """
    abstract type AbstractOrbit end
    # Main Keplerians structure must be built with an InitialConditions instance. Currently assumes that the hierarchy in the InitialConditions is correctly setup for the desired system. 

        # TODO: 
        # -figure out how to extract primary/secondary indices of objects based on heirarchy.
        # -use CartesianOutputs structure (which saves states) instead
    struct Keplerians{T <: AbstractFloat} <: AbstractOrbit
        # Cartesian coordinate positions and velocities of each body [dimension, body, step]
        xs::Array{T,3} # in AU
        vs::Array{T,3} # in AU/day
        intr::Integrator # Integration scheme (with t0, tmax, and h defined)
        nsteps::Int64  # The number of steps to integrate
        save_interval::Int64     # How often to save a step
        state::State{T}  # Current state of simulation
        ic::InitialConditions{T}  # InitialConditions for system
        names::Vector{String}
        states::Vector{State{T}}
        function Keplerians(intr::Integrator,ic::InitialConditions{T},nsteps::Int64)  where (T <: AbstractFloat)
            s=State(ic)
            nbody=s.n 
            xs=zeros(3,nbody,nsteps);
            vs=zeros(3,nbody,nsteps);
            names = ["body_$i" for i = 1:nbody]
            states = Array{State,1}(undef,nsteps+1)
            return new{T}(xs,vs,intr,nsteps,1,s,ic,names,states)
        end

    end
    """
        Keplerians{T <: AbstractFloat}(ic,t0,tmax,h)

    ### Arguments            
    - t0::Real : start time for integration
    - tmax::Real : end time for integration 
    - h::Real : integrator step size

    ### Optional
    - nsteps::Int64 : the number of points to save for plotting purposes, tmax = t0 + (h * nsteps).
            If nsteps is not provided, we set the value to an appropriate value (based on the period of outermost object).
    - save_interval::Int64   : if set to N, will save every Nth step.
            Can use large N for when you want to plot really long orbits.
            If using save_intergral, increase h to take less integration steps and make orbit paths more smooth.
    """
    function Keplerians(ic::InitialConditions{T},t0::T,tmax::T,h::T) where (T <: AbstractFloat)
        s=State(ic)
        nbody=s.n
        intr=Integrator(h,t0,tmax)
        nsteps=Integer(round(tmax-t0))
        xs=zeros(3,nbody,nsteps);
        vs=zeros(3,nbody,nsteps);
        names = ["body_$i" for i = 1:nbody]
        # states = []
        return Keplerians(xs,vs,intr,nsteps,1,s,ic,names)
    end
    """
    # If user wants to provide existing Integrator,

        Keplerians{T <: AbstractFloat}(intr,ic)
    """
    function Keplerians(intr::Integrator,ic::InitialConditions{T})  where (T <: AbstractFloat)
        s=State(ic)
        nbody=s.n
        P_max =  round(maximum(ic.elements[:,2]))
        nsteps=Integer(clamp(P_max,100,400))
        xs=zeros(3,nbody,nsteps);
        vs=zeros(3,nbody,nsteps);

        names = ["body_$i" for i = 1:nbody]
        # states = Array{State,1}(undef,nsteps)
        return Keplerians(xs,vs,intr,nsteps,1,s,ic,names)
    end
    # function Keplerians(intr::Integrator,ic::InitialConditions{T})  where (T <: AbstractFloat)
    #     s=State(ic)
    #     nbody=s.n
    # #     println("Number of integrations to plot (nsteps) not provided, setting value to outer planet period/h.")
    #     nsteps=Integer(round(maximum(P)/intr.h) )
    #     xs=zeros(3,nbody,nsteps);
    #     vs=zeros(3,nbody,nsteps);
    #      Keplerians(xs,vs,intr,nsteps,0,s,ic,[""])
    # end
    function Keplerians(intr::Integrator,ic::InitialConditions{T},nsteps::Int64,save_interval::Int64)  where (T <: AbstractFloat)
            @assert (nsteps >= save_interval)
            s=State(ic)
            nbody=s.n 
            nsaves=Integer(round(nsteps/save_interval))
            xs=zeros(3,nbody,nsaves);
            vs=zeros(3,nbody,nsaves);
            names = ["body_$i" for i = 1:nbody]
            # states = Array{State,1}(undef,nsaves)
            return Keplerians(xs,vs,intr,nsteps,save_interval,s,ic,names)
    end

    """ Allow keywords."""
    Keplerians(intr,ic;nsteps,save_interval)=Keplerians(intr,ic,nsteps,save_interval)
    Keplerians(ic;t0=0.5,tmax=100.0,h=0.5)=Keplerians(ic,t0,tmax,h)

    # Get masses (in solar masses) and periods (in days) of all bodies 
    get_masses(o) = o.state.m
    get_periods(o) = o.ic.elements[:,2]
    _get_heirarchy(o) = o.ic.ϵ

    """
        SimOrbits!(o::Keplerians{T})

    Does integration for nsteps
    """
    function SimOrbits!(o::Keplerians{T}) where T <: AbstractFloat
        xs=o.xs
        vs=o.vs
        nsteps=o.nsteps ;   
        save_interval=o.save_interval 
        states=[]
        ic=o.ic ;   
        intr=o.intr 
        s=o.state
        nbody=s.n
        init_state=deepcopy(s)
        o.states[1]=init_state
        t0 = s.t[1] # Initial time
        # Integrate in proper direction
        h = intr.h * NbodyGradient.check_step(t0,intr.tmax)
        tmax = t0 + (h * nsteps) 

        if !isinf(nsteps/save_interval)
            nsaves=Integer(round(nsteps/save_interval))
            for j in 1:nsaves    
                for i=1:nsteps
                    if nsteps%save_interval==0
                    # Take integration step and advance time
                        intr.scheme(s,h)
                        s.t[1] = t0 +  (h * i)
                        # Save State from current step
                        xs[:,:,j] = s.x[:,:]
                        vs[:,:,j] = s.v[:,:]
                        o.states[i+1]=deepcopy(s)
                    end
                end
            end
        else
            for i=1:nsteps
                # Take integration step and advance time
                intr.scheme(s,h)
                s.t[1] = t0 +  (h * i)
                # Save State from current step
                xs[:,:,i] = s.x[:,:]
                vs[:,:,i] = s.v[:,:]
                o.states[i+1] = deepcopy(s)
            end
        end
        return 
    end

    # Set up figure for static plot and animation
    function _setup(show_primary::Bool)
        fig,ax=plt.subplots(figsize=(4,4),dpi=100)
        ax.set_xlabel("[au]")
        ax.set_ylabel("[au]")
        if show_primary
            pc=ax.scatter([],[],s=25,color="black")
        end
        pc=ax.scatter([],[],s=25,color="black")
        return fig,ax
    end

    function _setup_set(show_primary::Bool)
        fig,ax=plt.subplots(2,2,figsize=(4,4),dpi=100)
        ax[2].set_xlabel("x [au]")
        ax[2].set_ylabel("y [au]")
        ax[1].set_ylabel("z [au]")
        ax[4].set_xlabel("z [au]")
        ax[1].xaxis.set_visible(false)
        ax[4].yaxis.set_visible(false)
        ax[3].set_frame_on(false)
        ax[3].set_xticks([]);ax[3].set_yticks([])
        if show_primary
            pc=ax[2].scatter([],[],s=25,color="black")
        end
            pc=ax[2].scatter([],[],s=25,color="black")
        return fig,ax
    end
    julia_colors=["#000000","#389826","#CB3C33","#9558B2","#4063D8"] # starts with black
    """
        make_plot(o::Keplerians{T}) 

    Plot static orbits, with orbital path as faded line and location of bodies at the initial state.
    # Optional Arguments
    - show_primary::Bool        : plots initial position of the primary star. Never plots the path. 
    - lw::Real                  : Linewidth
    - ms::Real                  : Markersize
    - colors::Vector{String}    : user provided list of colors to use for orbits.
    - legend::Bool              : Choose whether to show legend. (Default = false)

    """
    function make_plot(o::Keplerians{Float64},save::Bool=false,filename::String="../test/test_orbits.png",
    show_primary::Bool=true,lw::Real=1.5,ms::Real=25,use_colors::Bool=true;colors::Vector{String}=julia_colors,legend::Bool=false,figsize::Tuple{T,T}=(5,5),projection="xz") where T<:Real

        nbody=o.ic.nbody ; 
        nsteps=o.nsteps ;
        SimOrbits!(o)
        fig,ax=_setup(show_primary);
        fig.set_figheight(figsize[1]);
        fig.set_figwidth(figsize[2])
        proj_axis=1
        projection=="xy" ? proj_axis=2 : proj_axis=3 
        ax.set_title("Top Down",loc="left")
        if show_primary
            ax.scatter(o.state.x0[1,1],o.state.x0[proj_axis,1],marker="*",color="black",s=35*lw)
        end
        # Update colors.
        if use_colors && length(colors)<=nbody
            colors=julia_colors
        elseif !use_colors
            colors = repeat(["black"],nbody)
        else 
            @assert(length(colors) >= nbody)
            colors = colors
        end
        for body in 2:nbody
            ax.scatter(o.state.x[1,body],o.state.x[proj_axis,body],s=ms*lw,color=colors[body],zorder=3)
            ax.plot(o.xs[1,body,:],o.xs[proj_axis,body,:],alpha=0.4,lw=lw,color=colors[body])
        end
        if legend
            fig.legend()
        end
        fig.tight_layout()
        if save
            fig.savefig(filename,dpi)
        end
        return  fig,ax
    end


    """
        animate_plot(o::Keplerians{T}) 

    Animates orbits, with orbital path as transparent line that is same color as the body.
    # Optional Arguments
    - show_primary::Bool        : plots initial position of the primary star. Never plots the path. 
    - lw::Real                  : Linewidth
    - ms::Real                  : Markersize
    - colors::Vector{String}    : user provided list of colors to use for orbits.
    - legend::Bool              : Choose whether to show legend. (Default = false)
    - figsize::Tuple            
    """
    function animate_plot(o::Keplerians{Float64},save::Bool=false,filename::String="../test/test2.mp4",show_primary::Bool=true,lw::Real=1.5,ms::Real=50,use_colors::Bool=true;colors::Vector{String}=julia_colors,legend::Bool=false,figsize::Tuple{T,T}=(5,5)) where T<:Real
        nbody=o.ic.nbody;
        nsteps=o.nsteps;
        h=o.intr.h;
        save_interval=o.save_interval
        projection="xz"
        # if h < 1.0 && nsteps < 100
        # This will create an animation which may not cover a full orbit for long period objects.
        # end
        fig,ax=make_plot(o,false,filename,show_primary,lw,ms,use_colors;colors,legend,figsize)
        limits_x=maximum(abs.(o.xs[1,:,:]))*1.1
        limits_y=maximum(abs.(o.xs[3,:,:]))*1.1
        xlabel=ax.get_xlabel()
        ylabel=ax.get_ylabel()
        proj_axis=1
        projection=="xy" ? proj_axis=2 : proj_axis=3 
        for body in 2:nbody
            ax.plot(o.xs[1,body,:],o.xs[proj_axis,body,:],alpha=0.4,color=colors[body],lw=lw)
        end
        function update_plot(i)
            frame_mult=1 # can change to make faster animation. may have unexpected results
            if i*frame_mult < nsteps
            ax.clear();
            ax.set_ylim(-limits_y,limits_y);
            ax.set_xlim(-limits_x,limits_x)
            ax.set_ylabel(ylabel);
            ax.set_xlabel(xlabel);
            if show_primary
                ax.scatter(o.xs[1,1,i*frame_mult],o.xs[proj_axis,1,i*frame_mult],s=lw*ms*1.4,marker="*",color="black")
            end
            # Update colors.
            if use_colors && length(colors)<=nbody
                colors=julia_colors
            elseif !use_colors
                colors = repeat(["black"],nbody)
            else 
                colors = colors
            end
            for body in 2:nbody
            ax.plot(o.xs[1,body,:],o.xs[proj_axis,body,:],alpha=0.4,color=colors[body],lw=lw)
            ax.scatter(o.xs[1,body,i*frame_mult],o.xs[proj_axis,body,i*frame_mult],color=colors[body],s=ms*lw)
            end
            end
        return ax.get_lines()
        end
        if isinf(nsteps/save_interval) 
            # If save_interval or nsteps is not defined, save 200 steps. 
            nsaves=200
        else
            nsaves=nsteps
        end
        frames=[1:nsaves;]
        movie=anim.FuncAnimation(fig,update_plot,frames=frames,repeat=true,blit=true,interval=40)
        if save
        anim.FuncAnimation.save(movie,filename)
        end
        return 
    end

    animate_plot(o;save::Bool=false,filename::String="test.mp4")=animate_plot(o,save,filename)

    ### To fix:
    function make_plot_set(o::Keplerians{Float64},save::Bool=false,filename::String="../test/test_orbits.png",show_primary=true,lw::Real=1.5,ms::Real=25;colors::Vector{String}=[],legend::Bool=false,figsize::Tuple{T,T}=(4,4)) where T <:Real
        nbody=o.ic.nbody ; 
        nsteps=o.nsteps ;
        SimOrbits!(o)
        fig,ax=_setup_set(show_primary);
        fig.set_figheight(figsize[1]);
        fig.set_figwidth(figsize[2])    
        if show_primary
            ax[1].scatter(o.state.x0[1,1],o.state.x0[3,1],marker="*",color="black",s=35*lw)
            ax[2].scatter(o.state.x0[1,1],o.state.x0[2,1],marker="*",color="black",s=35*lw)
            ax[4].scatter(o.state.x0[3,1],o.state.x0[1,1],marker="*",color="black",s=35*lw)
        end
        # Update colors.
        if isempty(colors)
            colors=julia_colors
        elseif  colors==["black"]
            colors = repeat(["black"],nbody)
        end
        for body in 2:nbody
            ax[1].scatter(o.state.x[1,body],o.state.x[3,body],s=ms*lw,color=colors[body],zorder=3)
            ax[1].plot(o.xs[1,body,:],o.xs[3,body,:],alpha=0.5,lw=lw,color=colors[body])
            ax[2].scatter(o.state.x[1,body],o.state.x[2,body],s=ms*lw,color=colors[body],zorder=3,label=o.names[body])
            ax[2].plot(o.xs[1,body,:],o.xs[2,body,:],alpha=0.5,lw=lw,color=colors[body])
            ax[4].scatter(o.state.x[3,body],o.state.x[2,body],s=ms*lw,color=colors[body],zorder=3)
            ax[4].plot(o.xs[3,body,:],o.xs[2,body,:],alpha=0.5,lw=lw,color=colors[body])
        end
        if legend
            fig.legend()
        end
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0,hspace=0.0,right=0.98,top=0.98)

        if save
        fig.savefig(filename)
        end
        return  fig,ax
    end
    export Keplerians, SimOrbits! 
    export animate_plot,make_plot,make_plot_set


    function _pos_wrt_primary(r0,r2,P,t) # not in 3D, only for simple_plot()
        x0=r0[1];y0=r0[3]
        x=r2[1];y=r2[3]
        d=sqrt(abs((x-x0).^2 .- (y-y0) .^2))
        n=2*π/P
        new_x= x0 + (d *cos(n*t))
        new_y= y0 + (d * sin(n*t))
        return [new_x new_y]
    end

    function _new_pos(o,primary_indx,secondary_indx)
        P=o.ic.elements[secondary_indx,2] # mean motion of satellite about barycenter with host
        pos_about_parent=zeros(o.nsteps,2);
        for i=1:o.nsteps-1
        pos_about_parent[i,:]=_pos_wrt_primary(o.state.x[:,primary_indx] ,o.state.x[:,secondary_indx],P,i)
        end
        return pos_about_parent
    end

    function _transform_coord_system(Ω,ω,I,r)
        x,y,z=r
        X = x * (cos(Ω)* cos(ω) - sin(Ω)*sin(ω)*cos(I)) - y *(cos(Ω)* sin(ω) - sin(Ω)*cos(ω)*cos(I))
        Y = x * (sin(Ω)* cos(ω) - cos(Ω)*sin(ω)*cos(I)) - y *(sin(Ω)* sin(ω) - cos(Ω)*cos(ω)*cos(I))
        Z = x * sin(ω) * sin(I) + y * cos(ω)*sin(I)
        return X,Y,Z
    end

    function _update_3D_cartesian(o)
        intr= o.intr 
        nbody = o.ic.nbody
        elems = o.ic.elements
        ω = atan.(elems[:,4]./elems[:,5])
        I = elems[:,6]
        Ω = elems[:,7]
        new_xs=zeros(3,nbody,o.nsteps);
        # transform all objects to a common coordinate system
        for body=2:nbody
            for i=1:o.nsteps
                new_xs[:,body,i] .= _transform_coord_system(Ω[body],ω[body],I[body],o.xs[1:3,body,i])
            end
        end
        return new_xs
    end

end # module

### <===
    # function convert_to_elements(x::Vector{T}, v::Vector{T}, mstar::T, m_pl::T, t::T) where T <: AbstractFloat
    # G = 39.4845/(365.242 * 365.242) # AU^3 Msol^-1 Day^-2
    # R = norm(x)
    # v2 = dot(v, v)
    # hvec = cross(x, v)
    # h = norm(hvec)
    # Gm = (mstar+m_pl)*G # km3/s2/Msun
    # # (2.134)
    # a = 1. / (2/R - v2/Gm)
    # # (2.135)
    # e = sqrt(1 - h^2 / (Gm * a))
    # # (2.136)
    # hx,hy,hz=hvec
    # I = acos(clamp(hz / h, -1, 1))
    # # (2.137)
    # hz >= 0.0 ? hy *= -1 : hx *= -1
    # sinΩ = hx / (h * sin(I))
    # cosΩ = hy / (h * sin(I))
    # Ω = atan(sinΩ, cosΩ)
    # # (2.138)
    # sin_ω_plus_f = x[3] / (R * sin(I))
    # cos_ω_plus_f = (1 / cos(Ω)) * (x[1] / R + sin(Ω) * sin_ω_plus_f * cos(I))
    # ω_plus_f = atan(sin_ω_plus_f, cos_ω_plus_f)
    # # (2.139) assuming Rdot = (x.v)/R which is projection of velocity in radial direction
    # sinf = a * (1 - e^2) * dot(x, v) / (h * e * R)
    # cosf = (1 / e) * (a * (1 - e^2) / R - 1)
    # f = atan(sinf, cosf) 
    # ω = ω_plus_f - f
    # # (2.140)
    # # Step 1: Solve for E explicitly from Eq. (2.42)
    # cosE = (a - R) / (a * e)
    # E = acos(clamp(cosE, -1, 1))
    # # Adjust E based on radial velocity
    # Rdot = dot(x, v) / R
    # if Rdot < 0
    #     E = 2π - E
    # end

    # #  Mean anomaly M from Eq. (2.51)
    # M = E - e * sin(E)
    # #  Mean motion n from Eq. (2.26)
    # n_mean = sqrt(Gm / a^3)
    # P = 2π/n_mean
    # # True longitude from Eq. (2.118)
    # ϖ = Ω + ω
    # θ = ϖ + f
    # #  τ from M and n (Eq. 2.51 rearranged)
    # τ = t - M / n_mean
    # # Mean longitude, λ from Eq. (2.53)
    # λ =  n_mean * (t - τ)
    # return (
    #     a = a, # in astronomical units
    #     e = e,
    #     I = I,
    #     Ω = mod2pi(Ω),
    #     ω = mod2pi(ω),
    #     f = mod2pi(f),
    #     M = mod2pi(M),
    #     E = mod2pi(E),
    #     τ = τ, # in day
    #     n = n_mean, # per day
    #     P = P, # in days
    #     λ = λ,
    #     θ = θ,
    #     h = h # in 
    # )
    # end
# using Plots
# @recipe function f(ic_kep::ElementsIC)
#     pl_elems=[]
#     N=10
#     t=range(0,stop=2π,length=N)
#     s=State(ic_kep)
#     X,V = NbodyGradient.get_relative_positions(s.x,s.v,ic_kep)
#     μ=NbodyGradient.get_relative_masses(ic_kep)
#     aspect_ratio --> 1
#     xguide --> "x [au]"
#     yguide --> "y [au]"
# #     zguide --> "z [au]"
#     i = 1; b = 0
#     while i < ic_kep.nbody
#         if first(ic_kep.ϵ[i, :]) == zero(Float64)
#             b += 1
#         end
# #         @series begin
#         @series begin 
#             pl_elems=convert_to_elements(unpack.(X)[i+b],unpack.(V)[i+b],μ[i+b],0.0)
#     #         linewidth --> range(0, 10, length = N)
#             seriesalpha --> range(0, 1, length = N) 
#             markersizes --> range(0,5,length=10)
#             seriestype --> :scatter
#             x=pl_elems.a .*sin.(pl_elems.f .* collect(t)/pl_elems.n) 
#             y=pl_elems.a .* cos.(pl_elems.f .* collect(t)/pl_elems.n)
#             primary=true
#             x,y
#         end
#         if b > 0
#             b -= 2
#         elseif b < 0
#             i += 1
#         end
#         i += 1
#     end
# end
# plot(ic_kep)