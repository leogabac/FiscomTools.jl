
# ========== INFO ========== #
# Creator: GitHub user leogabac
# Contributors: none
# Last modification: 02/04/2022
# ========================== #

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
    include("linear_regression.jl")

    export incremental, bisection, newtonroot # For root approximation
    export dv, pdv, gradient # Numerical differentiation
    export ref, lgsolve # Matrices and systems of equations
    export trapz,intsimp # Integration techniques
    export Lhat, D1hat, D2hat # Differential operators of the first and second numerical derivative
    export linreg, predict, correlation, slope, intercept

end
