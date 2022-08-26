struct LinearModel
    slope
    intercept
end

slope(model::LinearModel) = model.slope
intercept(model::LinearModel) = model.intercept

"""
    linreg(x,y)

Computes a linear regression of the data vectors `x` and `y`.

"""
function linreg(x::Vector,y::Vector)
    @assert length(x) == length(y) "Vectors must have same length"

    n = length(x)
    A = [n sum(x); sum(x) sum(x.^2)]
    b = [sum(y), sum(x .* y)]
    (b,m) = tuple(A\b...)
    return LinearModel(m,b)
end

function predict(model::LinearModel, x::Number)
    m = slope(model)
    b = intercept(model)
    return m*x + b
end

function correlation(model::LinearModel, x::Vector,y::Vector)
    y_bar = sum(y) / length(y)
    St = sum( (y .- y_bar).^2 )
    Sr = sum( (y .-  map(xi->predict(model,xi),x) ).^2 )
    return (St - Sr)/St
end
