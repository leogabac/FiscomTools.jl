function HermitePoly(n::Integer,x)
    upper = floor(n/2) |> Int
    Hnt = sum( (-1)^k/factorial(k)/factorial(n-2k)*(2*x)^(n-2k) for k in 0:upper )
    return factorial(n) * Hnt
end

function LaguerrePoly(n::Integer,alpha::Integer,x)
    sum( binomial(n+alpha,n-k) * (-x)^k / factorial(k) for k in 0:n)
end


# =============================================================================
# Beams by benjamimoso
# =============================================================================

"""
    LaguerreGaussBeam(x, y, z, w0, phi, lambda, l, p)
    Laguerre-Gaussian beam
"""
function LaguerreGaussBeam(x::Real, y::Real, z::Real, w0::Real, phi::Real, lambda::Real, l, p)
    LG::ComplexF64=0.0 + im*0.0

    zr = pi*(w0^2)/lambda
    wz2 = (w0^2) * (1 + (z/zr)^2)
    wz = sqrt(wz2)
    rr2=(x^2 + y^2)/(wz2)
    C = sqrt(2 * factorial(p) / (pi*factorial(p + abs(l))))
    k = 2*pi/lambda

    if p==0
        LG = (C/wz) * ((sqrt(2*rr2))^abs(l)) * exp(-rr2) *
            1.0 * exp(-im*rr2*z/zr) * exp(im*(l*atan(y,x)+phi)) *
            exp(-im*(2*p+abs(l)+1) * atan(z,zr))
    else
        LG = (C/sqrt(wz)) * ((sqrt(2*rr2))^abs(l)) * exp(-rr2) *
            LaguerrePoly.(p, abs(l), 2*rr2) * exp(-im*rr2*z/zr) * exp(im*(l*atan(y,x)+phi)) *
            exp(-im*(2*p+abs(l)+1) * atan(z,zr))
    end

    return LG::ComplexF64
end

"""
    HermiteGaussBeam(x, y, z, w0, phi, lambda, m, n)
    Hermite-Gaussian beam
"""
function HermiteGaussBeam(x::Real, y::Real, z::Real, w0::Real, phi::Real, lambda::Real, m::Integer, n::Integer)
    HG::ComplexF64=0.0 + im*0.0

    zr = pi*(w0^2)/lambda
    wz2 = (w0^2) * (1 + (z/zr)^2)
    wz = sqrt(wz2)
    rr2=(x^2 + y^2)/(wz2)
    C = sqrt(2/pi) * (2.0)^(-(m+n)/2) * (1/sqrt(factorial(n)*factorial(m)*(w0^2)))
    k = 2*pi/lambda

    HG= C * HermitePoly(m,sqrt(2)*x/wz) * HermitePoly(n,sqrt(2)*y/wz) * exp(-im*rr2*z/zr) *
        exp(-rr2) * exp(-im*(m+n+1)*atan(z,zr)) * exp(im*(phi))   # Im not sure about adding the dot in LaguerrePoly...
    return HG::Complex{Float64}
end
