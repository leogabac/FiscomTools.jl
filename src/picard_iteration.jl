function picard(F, x_0::Float64; t_i = 0, t_f = 10, Δt = 0.1, k = 1000)
    time = collect(range(t_i, stop = t_f, step = Δt))
    len = length(time);
    x = [x_0 for _ in time]
    for _ in 2:k+1
        x_prime = F.(x)
        picard_integral = 0;
        for t in 2:len
            picard_integral += (x_prime[t] + x_prime[t-1]) * (Δt / 2)
            x[t] = x_0 + picard_integral
        end
    end
    return time, x
end

function picard(F, X_0::Tuple; t_i = 0, t_f = 10, Δt = 0.1, k = 1000)
    time = collect(range(t_i, stop = t_f, step = Δt))
    len = length(time)
    X = [float.(X_0) for _ in time]
    for _ in 2:k+1
        X_prime = [F(i...) for i in X]
        picard_int = X_0 .* 0
        for t in 2:len
            picard_int = picard_int .+ (X_prime[t] .+ X_prime[t-1]) .* (Δt / 2)
            X[t] = X_0 .+ picard_int
        end
    end
    return time, X
end

function picard(F, X_0::Tuple, t::Vector{Float64}, X; k = 100)
    Δt = t[2] - t[1];
    len = length(t)
    for _ in 2:k+1
        X_prime = [F(i...) for i in X]
        picard_int = X_0 .* 0
        for t in 2:len
            picard_int = picard_int .+ (X_prime[t] .+ X_prime[t-1]) .* (Δt / 2)
            X[t] = X_0 .+ picard_int
        end
    end
    return X
end
