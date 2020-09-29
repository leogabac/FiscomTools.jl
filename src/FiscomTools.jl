module FiscomTools
include("test.jl")
include("incremental.jl")
include("bisection.jl")
include("diff.jl")
include("newton.jl")
export bisection,incremental,diff,newtonR
end
