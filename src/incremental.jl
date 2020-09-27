function incremental(f,min::Int64,max::Int64,N::Int64)
    t = range(min,stop = max,length = N) #makes the domain
    intervalo = Array{Float64}(undef,N-1,2) #prepares all intervals
    intervalo[:,1], intervalo[:,2] = t[1:end-1], t[2:end]
    return intervalo[(f.(t[1:end-1])) .* (f.(t[2:end])) .< 0 ,:] #chooses all intervals that the product of the edges is negative
end

function incremental(f,min::Float64,max::Float64,N::Int64)
    t = range(min,stop = max,length = N) #makes the domain
    intervalo = Array{Float64}(undef,N-1,2) #prepares all intervals
    intervalo[:,1], intervalo[:,2] = t[1:end-1], t[2:end]
    return intervalo[(f.(t[1:end-1])) .* (f.(t[2:end])) .< 0 ,:]
end