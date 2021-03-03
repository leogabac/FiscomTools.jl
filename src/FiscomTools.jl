module FiscomTools
include("test.jl")
include("incremental.jl")
include("bisection.jl")
include("ndiff.jl")
include("newton.jl")
include("gauss.jl")
include("intsimp.jl")
include("lucrout.jl")
include("trapz.jl")
include("diff_operators.jl")

export incremental, bisection_it, bisection_err, newtonR # For root approximation
export ndiff # Numerical differentiation
export redGauss, lgsolve, lucrout # Matrices and systems of equations
export trapz,intsimp # Integration techniques
export d1hat, d2hat # Differential operators of the first and second numerical derivative
end
