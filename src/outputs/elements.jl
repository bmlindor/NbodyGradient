# Converting cartesian coordinates to orbital elements

"""

A 3d point in cartesian space
"""
struct Point{T<:AbstractFloat}
    x::T
    y::T
    z::T
end

# Constructors
Point(x::AbstractVector) = Point(x...)
Point(x::AbstractMatrix) = [Point(x[:,i]) for i in eachindex(x[1,:])]
Point(x::Real) = Point(float(x),float(x),float(x))
unpack(x::Point) = [x.x,x.y,x.z]

# Overloads for products
LinearAlgebra.dot(x::Point,v::Point) = x.x*v.x + x.y*v.y + x.z*v.z
LinearAlgebra.dot(x::Point) = x.x*x.x + x.y*x.y + x.z*x.z
LinearAlgebra.cross(x::Point,v::Point) = Point(x.y*v.z - x.z*v.y,-(x.x*v.z - x.z*v.x),x.x*v.y - x.y*v.x)

""" Calculate the relative positions from the A-Matrix. """
function get_relative_positions(x,v,ic::InitialConditions)
    n = ic.nbody
    X = zeros(3,n)
    V = zeros(3,n)

    X .= permutedims(ic.amat*x')
    V .= permutedims(ic.amat*v')

    return Point(X), Point(V)
end

""" Get relative masses(*G) from initial conditions. """
function get_relative_masses(ic::InitialConditions)
    N = length(ic.m)
    M = zeros(N-1)
    G = 39.4845/(365.242 * 365.242) # AU^3 Msol^-1 Day^-2
    for i in 1:N-1
        for j in 1:N
            M[i] += abs(ic.ϵ[i,j])*ic.m[j]
        end
    end
    return G .* M
end

# Position and velocity magnitudes
mag(x) = sqrt(dot(x))
mag(x,v) = sqrt(dot(x,v))
Rdotmag(R,V,h) = sqrt(V^2 - (h/R)^2)

function hvec(r,rdot)
    hx,hy,hz = unpack(cross(r,rdot))
    hz >= 0.0 ? hy *= -1 : hx *= -1
    return Point(hx,hy,hz)
end

function calc_Ω(hx,hy,h,I)
    sinΩ = hx/(h*sin(I))
    cosΩ = hy/(h*sin(I))
    return atan(sinΩ,cosΩ)
end

function calc_ω(x,R,Rdot,I,Ω,a,e,h)
    # Find ω + f
    wpf = 0.0
    if I != 0.0
        sinwpf = x.z/(R*sin(I))
        coswpf = ((x.x/R) + sin(Ω)*sinwpf*cos(I))/cos(Ω)
        wpf = atan(sinwpf,coswpf)
    end
    # Find f
    sinf = a*Rdot*(1.0 - e^2)/(h*e)
    cosf = (a*(1.0-e^2)/R - 1.0)/e
    f = atan(sinf,cosf)

    return wpf - f,f
end

function convert_to_elements(x::Point,v::Point,Gmm::T,t::T) where T<: AbstractFloat
    R = mag(x)
    V = mag(v)
    h = mag(hvec(x,v))
    hx,hy,hz = unpack(hvec(x,v))

    Rdot = sign(dot(x,v))*Rdotmag(R,V,h)

    a = 1.0/((2.0/R) - (V*V)/(Gmm))
    e = sqrt(1.0 - (h*h/(Gmm*a)))
    I = acos(hz/h)

    # Make sure Ω is defined
    I != 0.0 ? Ω = calc_Ω(hx,hy,h,I) : Ω = 0.0

    ω,f = calc_ω(x,R,Rdot,I,Ω,a,e,h)
    P = (2.0*π)*sqrt(a*a*a/Gmm)

    n = 2π/P
    esinω, ecosω= e.*sincos(ω)

    tp = (-sqrt(1.0-e*e)*ecosω/(n*(1.0-esinω)) - (2.0/n)*atan(sqrt(1.0-e)*(esinω+ecosω+e), sqrt(1.0+e)*(esinω-ecosω-e))) % P

    return [P,0.0,ecosω,esinω,I,Ω,a,e,ω,tp]
end

function convert_to_elements(x::Vector{T}, v::Vector{T}, Gm::T, t::T) where T <: AbstractFloat
    R = norm(x)
    v2 = dot(v, v)
    hvec = cross(x, v)
    h = norm(hvec)

    # (2.134)
    a = 1. / (2/R - v2/Gm)
    # (2.135)
    e = sqrt(1 - h^2 / (Gm * a))
    # (2.136)
    hx,hy,hz=hvec
    I = acos(clamp(hz / h, -1, 1))
    # (2.137)
    hz >= 0.0 ? hy *= -1 : hx *= -1
    sinΩ = hx / (h * sin(I))
    cosΩ = hy / (h * sin(I))
    Ω = atan(sinΩ, cosΩ)
    # (2.138)
    sin_ω_plus_f = x[3] / (R * sin(I))
    cos_ω_plus_f = (1 / cos(Ω)) * (x[1] / R + sin(Ω) * sin_ω_plus_f * cos(I))
    ω_plus_f = atan(sin_ω_plus_f, cos_ω_plus_f)
    # (2.139) assuming Rdot = (x.v)/R which is projection of velocity in radial direction
    sinf = a * (1 - e^2) * dot(x, v) / (h * e * R)
    cosf = (1 / e) * (a * (1 - e^2) / R - 1)
    f = atan(sinf, cosf)
    ω = ω_plus_f - f
    # (2.140)
    # Step 1: Solve for E explicitly from Eq. (2.42)
    cosE = (a - R) / (a * e)
    E = acos(clamp(cosE, -1, 1))
    # Adjust E based on radial velocity
    Rdot = dot(x, v) / R
    if Rdot < 0
        E = 2π - E
    end
    #  Mean anomaly M from Eq. (2.51)
    M = E - e * sin(E)
    #  Mean motion n from Eq. (2.26)
    n_mean = sqrt(Gm / a^3)
    #  τ from M and n (Eq. 2.51 rearranged)
    τ = t - M / n_mean
    return (
        a = a,
        e = e,
        I = I,
        Ω = mod2pi(Ω),
        ω = mod2pi(ω),
        f = mod2pi(f),
        M = mod2pi(M),
        E = mod2pi(E),
        τ = τ,
        n = n_mean,
        h = h
    )
end

function get_orbital_elements(s::State{T},ic::InitialConditions{T}) where T<:AbstractFloat
    # For comparison purposes only.
    elems = Elements{T}[]
    μ = get_relative_masses(ic)
    X,V = get_relative_positions(s.x,s.v,ic)
    N = ic.nbody

    push!(elems,Elements(m=ic.m[1]))

    i = 1; b = 0
    while i < N

        # Check if new binary
        if first(ic.ϵ[i,:]) == zero(T)
            b+=1
        end
        t0 = ic.elements[i+b,3]
        new_elems = convert_to_elements(X[i+b],V[i+b],μ[i+b],0.0)
        push!(elems,Elements(ic.m[i+1],new_elems...))

        # Compensate for new binary in step
        if b > 0
            b -= 2
        elseif b < 0
            i += 1
        end

    i+=1
    end
    return elems
end
""" Get orbital elements from initial conditions (and provided t0) if x and v are arrays. """
function get_orbital_elements_actual(s::State{T}, ic::InitialConditions{T}) where T <: AbstractFloat
    elems = Elements{T}[]
    μs = get_relative_masses(ic)
    X,V = get_relative_positions(s.x,s.v,ic)

    push!(elems, Elements(m=ic.m[1]))  # Central body
    i = 1; b = 0
    while i < ic.nbody
        if first(ic.ϵ[i, :]) == zero(T)
            b += 1
        end
        t0 = ic.elements[i+b,3]
        a, e, I, Ω, ω, f, M, E, τ, n, h = convert_to_elements(unpack(X[i+b]), unpack(V[i+b]), μs[i+b], s.t[1])
        # println("Additional elements. \n f: $f \n M: $M \n E: $E ")
        push!(elems, Elements(ic.m[i+1], 2π / n, 0.0, e*cos(ω), e*sin(ω), I, Ω, a, e, ω, τ))
        if b > 0
            b -= 2
        elseif b < 0
            i += 1
        end
        i += 1
    end
    return elems
end