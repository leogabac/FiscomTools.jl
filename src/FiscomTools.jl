module FiscomTools
include("test.jl")
include("bisection.jl")
include("numdiff.jl")
include("newton.jl")
include("gauss.jl")
include("intsimp.jl")
include("lucrout.jl")
include("trapz.jl")
include("diff_operators.jl")

export incremental, bisection, newtonroot # For root approximation
export dv, pdv, gradient # Numerical differentiation
export ref, lgsolve, lucrout # Matrices and systems of equations
export trapz,intsimp # Integration techniques
export d1hat, d2hat # Differential operators of the first and second numerical derivative
end
