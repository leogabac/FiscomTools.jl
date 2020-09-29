module FiscomTools
include("test.jl")
include("incremental.jl")
include("bisection.jl")
include("numdiff.jl")
include("newton.jl")
export bisection,incremental,numdiff,newtonR
end
