#= ==========================================================================================
n-th dimensional interpolation
========================================================================================== =#


# Method for a variable number of inputs, which are the variables themselves
# it returns a vector of constants

function interpolate(z, variables::Array{Float64,1}...)
    m = Array{Int64}(undef, length(variables));
    for i in 1:length(variables)
        m[i] = length(unique(variables[i]));
    end
    Vt = 1;
    for i in 1:length(variables)
        x = unique(variables[i]);
        V = Array{Float64}(undef, m[i], m[i]);
        for j in 0:m[i] - 1
            for k in 1:m[i]
                V[k, j + 1] = x[k]^j; # Si truena fue esto 
            end
        end
        Vt = kron(Vt, inv(V))
    end
    return Vt * z, m
end

# Method that takes a data array as an input
# it returns a vector of constants

function interpotale(data::Array{Float64,2})
    m = Array{Int64}(undef, size(data, 2) - 1);
    for i in 1:size(data, 2) - 1
        m[i] = length(unique(data[:,i]));
    end
    Vt = 1;
    for i in 1:size(data, 2) - 1
        x = unique(data[:,i]);
        V = Array{Float64}(undef, m[i], m[i]);
        for j in 0:m[i] - 1
            for k in 1:m[i]
                V[k, j + 1] = x[k]^j; # Si truena fue esto 
            end
        end
        Vt = kron(Vt, inv(V))
    end
    return Vt * data[:,end], m
end

# Function that writes the polynomial in LaTeX
function pol2tex(c::Array{Float64,1}, m::Array{Int64,1})
    M = 1;
    for i in 1:length(m)
        M *= m[i];
    end
    c = round.(c, sigdigits=5);
    for i in 1:length(c)
        if abs(c[i]) < 1e-10
            c[i] = 0;
        end
    end

    function foranid(n::Int64, i::Array{Int64,1})
        if n == length(m) + 1
            if c[k] != 0
                if abs(c[k]) != 1
                    model[k] = string(abs(c[k]));
                else
                    model[k] = "";
                end
                for j in 1:n - 1
                    if i[j] != 0                
                        if i[j] == 1
                            model[k] *= "x_" * "{" * string(j) * "}"; 
                        else
                            model[k] *= "x_" * "{" * string(j) * "}^{" * string(i[j]) * "}"; 
                        end
                    end
                end
            end
            k += 1;
        else
            for q in 0:m[n] - 1
                model = foranid(n + 1, [i;q]);
            end
        end
        return model
    end
    k = 1;
    model = Array{String}(undef, M, 1);
    for q in 0:m[1] - 1
        model = foranid(2, [q]);
    end

    first = false;
    if c[1] != 0
        modelo = string(c[1]);
        first = true;
    end

    for i in 2:M
        if first == false
            if c[i] > 0
                modelo = model[i];
            elseif c[i] <  0
                modelo = " -" + model[i];
            end
            first = true;
        else
            if c[i] > 0             
                modelo *= " + " * model[i];
            elseif c[i] <  0
                modelo *= " - " * model[i];
            end
        end
    end
    for i = 10:-1:1
        modelo = replace(modelo, "e-" * string(i) => "\times 10^{-" * string(i) * "}");
    end
    return modelo
end

# Function that substitutes the variable names of the LaTeX string into sth desired by the user
function sustex(modelo::String, vars::Array{String,1})
    for i = length(vars):-1:1
        modelo = replace(modelo, "x_{" * string(i) * "}" => vars[i]);
    end
    return modelo
end

#= # Función del polinomio
function polin(c::Array{Float64,1}, m::Array{Int64,1}, x::Array{Float64,1})

    function foranid(n::Int64, i::Array{Int64,1})
        if n == length(m) + 1
            term = c[k];
            for j in 1:n - 1
                term *= x[j]^i[j]
            end
            z += term;
            k += 1;
        else
            for q in 0:m[n] - 1
                z = foranid(n + 1, [i;q]);
            end
        end
        return z
    end
    k, z = 1, 0; # initialize
    for q in 0:m[1] - 1
        z = foranid(2, [q]);
    end
    return z;
end =#