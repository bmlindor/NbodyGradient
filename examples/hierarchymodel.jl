module AgolModels
using Photodynamics,NbodyGradient
# For fitting, define photodynamic model
struct PhotometryModel{T<:Real} 
    t0::T
    tmax::T
    h::T
    # Structures for computing the various models
    ic::ElementsIC{T}
    tt::TransitTiming{T}
    ts::TransitSeries{T,Photodynamics.ComputedTimes}
    d::NbodyGradient.Derivatives{T}
    lc::Lightcurve{T}
    dlc::dLightcurve{T}
    intr::Integrator{T}
    J::Matrix{T}                    # Jacobian matrix
    jac_inds_elems::Vector{Int64}   # varied parameter indices
    # H::Matrix{T}                
#     jac_inds_q::Vector{Int64}
end

function create_elements_matrix(θ, N)
#     @assert (length(θ[1:end]) % (N-1)) == 1 "Lengths need to agree (N and size of elements matrix)."
    #   Stellar mass at the start; then elements ordered by planet, with planet masses last
    elements_matrix = zeros(N, 7)
    elements_masses = permutedims(reshape(θ[2:(N-1)*7+1], 7, N-1))
    elements_matrix[2:end, 1] .= elements_masses[:, end] # pl_masses
    elements_matrix[2:end, 2:end] .= elements_masses[:, 1:end-1]
    elements_matrix[1,1] = θ[1] 
    return elements_matrix
end

# For now assume we're doing every parameter. We then wrap in a function that fixes the ones we want.
function compute_photometry(model::PhotometryModel, θ; tol=1e-6, maxdepth=6)
    # Adopt hierarchy
    H= model.ic.ϵ 
    N = size(H,1)
    t0 = model.t0
    tmax = model.tmax
    # Replace parameters in lightcurve
    # Start with the lightcurve parameters
    model.lc.rstar .= θ[end]
    model.lc.u_n .= θ[end-2:end-1]
    model.lc.k .= θ[7*(N-1) + 2:end-3]

    # Now get new initial conditions from θ
    # Get the elements matrix 
    elements = create_elements_matrix(θ[1:7*(N-1)+1], N)
    ic = ElementsIC(t0, H, elements)
    @show ic
    # Reset computed transit times
    NbodyGradient.zero_out!(model.tt)
    model.ts.count .= 0
    
    # Setup integrator and run N-body integrator
    s = State(ic)
    
    # Compute the photometry
    try
        model.intr(s, model.ts, model.tt; grad=false)
        compute_lightcurve!(model.lc, model.ts; tol=tol, maxdepth=maxdepth)
    catch e
        model.lc.flux .= 1.0
        return model.lc.flux
    end
    return model.lc.flux #copy(lc.flux)
end

function grad_compute_photometry(model::PhotometryModel, θ; tol=1e-6, maxdepth=6)
    H= model.ic.ϵ 
    N = size(H,1)
    t0 = model.t0
    tmax = model.tmax

    lc = model.lc
    dlc = model.dlc
    # Replace parameters in both lightcurves
    lc.rstar .= θ[end]
    dlc.rstar .= θ[end]
    lc.u_n .= θ[end-2:end-1]
    dlc.u_n .= θ[end-2:end-1]
    lc.k .= θ[7*(N-1)+ 2:end-3]
    dlc.k .= θ[7*(N-1)+ 2:end-3]

    # Now get new initial conditions from θ
    # Get the elements matrix
    elements = create_elements_matrix(θ[1:7*(N-1)+1], N)
    ic = ElementsIC(t0,H, elements)

    # Reset computed transit times
    NbodyGradient.zero_out!(model.tt)
    model.ts.count .= 0
    
    # Setup integrator and run N-body integrator
    s = State(ic)
    
    # Compute the photometry and derivatives
    try
        model.intr(s, model.ts, model.tt; grad=true)
        compute_lightcurve!(model.dlc, model.ts; tol=tol, maxdepth=maxdepth)
        lc.flux .= dlc.flux
        
        # Transform derivatives from wrt Cartesian back to wrt orbital elements
        transform_to_elements!(s,model.dlc)
        # BL: Hurum said there was a bug in this transformation
        # J = dfd{P1,t01,ec1,es1,I1,Ω1,m1...mN,k1,...kN,u1,u2,rs}
        # Collect Jacobian arrays
        # jac_flux = hcat(dlc.dfdr, dlc.dfdu, dlc.dfdelements[:, 7:end], dlc.dfdk)

        model.J .= hcat(dlc.dfdelements[:, 7:end], dlc.dfdk, dlc.dfdu, dlc.dfdr)
    catch e
        @warn "Error in computing the photometry."
        model.lc.flux .= 1.0
        model.dlc.flux .= 1.0
        model.J .= 0.0
    end
    return model.dlc.flux , model.J[:,model.jac_inds_elems] #copy(lc.flux)
end

function build_model_ic(t0,θ,H)
    N=size(H,1)
    elems = create_elements_matrix(θ,N)
    return ElementsIC(t0,H,elems)
end

# Initialize model with data and vector of elements θ  
function PhotometryModel(t0::T,tmax::T,cadence::T,tobs::Vector{T},fobs::Vector{T},ferr::Vector{T},jac_inds_elems::Vector{Int64},H::Matrix{<:Real},θ::Vector{T}) where T <: AbstractFloat
    N=size(H,1)
    @assert length(jac_inds_elems)<=7*N
    ic_model=build_model_ic(t0,θ,T.(H))
    obs_duration=tmax - t0
    rstar = θ[end]
    u_n = θ[end-2:end-1]
    k = θ[7*(N-1)+ 2:end-3]
    # Initialize flux model, creating Transit and Lightcurve structs
    tt_model = TransitTiming(obs_duration+t0, ic_model)
    ts_model = TransitSeries(obs_duration+t0, ic_model)
    lc_model = Lightcurve(cadence, copy(tobs), copy(fobs), ferr,u_n,k,rstar);
    dlc_model = dLightcurve(cadence, copy(tobs), copy(fobs), ferr, u_n,k,rstar);
    d = NbodyGradient.Derivatives(Float64, N);
    intr=Integrator(cadence, t0,tmax)

    return PhotometryModel(t0,tmax,cadence,
    ic_model,  tt_model,ts_model,   d,    lc_model, dlc_model, intr,
    zeros(length(tobs),length(θ)),jac_inds_elems)
end

## Let user pass hierarchy vector
PhotometryModel(t0::T,tmax::T,cadence::T,tobs::Vector{T},fobs::Vector{T},ferr::Vector{T},jac_inds_elems::Vector{Int64},H::Vector{Int64},θ::Vector{T}) where T <: AbstractFloat = PhotometryModel(t0,tmax,cadence,tobs,fobs,ferr,jac_inds_elems,NbodyGradient.hierarchy(H),θ)

## Let user specify only the number of bodies; the hierarchy is filled in as fully-nested; 
PhotometryModel(t0::T,tmax::T,cadence::T,tobs::Vector{T},fobs::Vector{T},ferr::Vector{T},jac_inds_elems::Vector{Int64},H::Int64,θ::Vector{T}) where T <: AbstractFloat = PhotometryModel(t0,tmax,cadence,tobs,fobs,ferr,jac_inds_elems,[H, ones(Int64,H-1)...],θ)

Base.show(io::IO,::MIME"text/plain",model::PhotometryModel{T}) where {T} = begin
# println(io,"PhotometryModel{$T}\n Hierarchy: " Elements: "); show(io,"text/plain",model.ic.elements); end;
)
# struct TimingModel{T<:Real}
#     N::Int
#     t0::T
#     tmax::T
#     h::T
#     # Structures for computing the various models
#     ic::ElementsIC{T}
#     tt::TransitTiming{T}
#     ts::TransitSeries{T,Photodynamics.ComputedTimes}
#     d::NbodyGradient.Derivatives{T}
#     # H::Matrix{Int64} # heirarchy 
# end
# struct LightcurveModel{T<:Real}
#     N::Int
#     t0::T
#     tmax::T
#     h::T
#     lc::Lightcurve{T}
#     dlc::dLightcurve{T}
# end
function compute_timing()
    
end
export PhotometryModel
end # module