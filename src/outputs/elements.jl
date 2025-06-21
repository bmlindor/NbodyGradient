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
unpack(x::Point) = (x.x,x.y,x.z)

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

    return wpf - f
end

# function convert_to_elements(x,v,M)
#     R = mag(x)
#     V = mag(v)
#     h = mag(hvec(x,v))
#     hx,hy,hz = unpack(hvec(x,v))
#     Rdot = sign(dot(x,v))*Rdotmag(R,V,h)

#     Gmm = M

#     a = 1.0/((2.0/R) - (V*V)/(Gmm))
#     e = sqrt(1.0 - (h*h/(Gmm*a)))
#     I = acos(hz/h)

#     # Make sure Ω is defined
#     I != 0.0 ? Ω = calc_Ω(hx,hy,h,I) : Ω = 0.0

#     ω = calc_ω(x,R,Rdot,I,Ω,a,e,h)
#     P = (2.0*π)*sqrt(a*a*a/Gmm)

#     n = 2π/P
#     ecosω, esinω = e.*sincos(ω)
#     tp = (-sqrt(1.0-e*e)*ecosω/(n*(1.0-esinω)) - (2.0/n)*atan(sqrt(1.0-e)*(esinω+ecosω+e), sqrt(1.0+e)*(esinω-ecosω-e))) % P

#     return [P,0.0,ecosω,esinω,I,Ω,a,e,ω,tp]
# end
function convert_to_elements(x::Point, v::Point, Gmm::Float64, t::Float64)
    # (2.126) R^2 = x^2 + y^2 + z^2
    r = mag(x)
    # (2.127) V^2 = vx^2 + vy^2 + vz^2
    v2 = dot(v)
    # (2.129) h = r × v
    hvec = cross(x, v)
    h = mag(hvec)
    # Specific orbital energy (2.28): ε = v^2/2 - Gmm/r
    ε = 0.5 * v2 - Gmm / r
    # (2.134) a = [2/r - v^2 / Gmm]^-1
    if abs(ε) < 1e-10
        # Parabolic
        a = Inf
        e = 1.0
    else
        a = -Gmm / (2ε)
        # (2.135) e = sqrt(1 - h^2 / (Gmm * a))
        e = sqrt(1.0 + 2ε * h^2 / Gmm^2)
    end
    # (2.131) inclination I = arccos(hz / h)
    hx, hy, hz = unpack(hvec)
    I = acos(clamp(hz / h, -1, 1))
    # (2.132-2.133) Ω from hvec
    Ω = I > 1e-8 ? atan2(hx, -hy) : 0.0
    evec = (cross(v, hvec) .- Point(x.x * Gmm / r, x.y * Gmm / r, x.z * Gmm / r)) .* (1 / Gmm)
    ω = 0.0
    if I > 1e-8
        nvec = Point(-hy, hx, 0.0)
        n = mag(nvec)
        if n > 0
            ω = acos(clamp(dot(nvec, evec) / (n * e), -1, 1))
            if evec.z < 0
                ω = 2π - ω
            end
        end
    end
    # (2.20) true anomaly f = angle between evec and r
    cosf = clamp(dot(evec, x) / (e * r), -1, 1)
    sinf = dot(x, v) / (e * sqrt(Gmm * a))
    f = atan(sinf, cosf)
    # (2.42) r = a(1 - e cosE)
    # (2.50, 2.51) Eccentric anomaly and mean anomaly
    E = 0.0
    M = 0.0
    if e < 1.0
        E = 2 * atan( sqrt((1 - e) / (1 + e)) * tan(f / 2) )
        M = E - e * sin(E)
    elseif e ≈ 1.0
        M = 0.0
    else
        # Hyperbolic orbit
        F = asinh( sqrt((e - 1) / (1 + e)) * tanh(f / 2) )
        M = e * sinh(F) - F
    end
    # (2.39) and (2.140)  mean motion , time of periapsis
    n = sqrt(Gmm / abs(a^3))
    τ = t - M / n
    return (
        a=a,
        e=e,
        I=I,
        Ω=Ω % (2π),
        ω=ω % (2π),
        f=f % (2π),
        M=M % (2π),
        E=E,
        τ=τ,
        n=n,
        h=h
    )
end
# function get_orbital_elements(s::State{T},ic::InitialConditions{T}) where T<:AbstractFloat
#     elems = Elements{T}[]
#     μ = get_relative_masses(ic)
#     X,V = get_relative_positions(s.x,s.v,ic)
#     N = ic.nbody

#     push!(elems,Elements(m=ic.m[1]))

#     i = 1; b = 0
#     while i < N

#         # Check if new binary
#         if first(ic.ϵ[i,:]) == zero(T)
#             b+=1
#         end

#         new_elems = convert_to_elements(X[i+b],V[i+b],μ[i+b])
#         push!(elems,Elements(ic.m[i+1],new_elems...))

#         # Compensate for new binary in step
#         if b > 0
#             b -= 2
#         elseif b < 0
#             i += 1
#         end

#     i+=1
#     end
#     return elems
# end
function get_orbital_elements(s::State{T}, ic::InitialConditions{T}) where T <: AbstractFloat
    elems = Elements{T}[]
    μs = get_relative_masses(ic)
    X, V = get_relative_positions(s.x, s.v, ic)
    push!(elems, Elements(m=ic.m[1]))  # Central body
    i = 1; b = 0
    while i < ic.nbody
        if first(ic.ϵ[i, :]) == zero(T)
            b += 1
        end
        a, e, I, Ω, ω, f, M, E, τ, n, h = convert_to_elements(X[i+b], V[i+b], μs[i+b], s.t[1])
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
