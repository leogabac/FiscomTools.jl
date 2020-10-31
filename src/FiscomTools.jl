module FiscomTools
include("test.jl")
include("incremental.jl")
include("bisection.jl")
include("numdiff.jl")
include("newton.jl")
include("gauss.jl")
include("intsimp.jl")
include("lucrout.jl")
include("trapz.jl")

export incremental, bisection, newtonR # For root approximation
export numdiff # Numerical differentiation
export redGauss, lgsolve, lucrout # Matrices and systems of equations
export trapz,intsimp # Integration techniques
end
