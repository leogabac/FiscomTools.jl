function redGauss(A)
    # First, we make sure to convert all to a Float64 array
    A = convert(Array{Float64,2},A)
    numRow, numCol = size(A); # Then retrieve some important variable names
    #
    for k in 1:numCol-1
        # Is the pivot zero?
        if A[k,k] == 0.0 # Yes
            nonzero = findall(x->x!==0.0, A[:,k]) # Find all non-zero elements in the current column
            selected = nonzero[1] # We then select the first non-zero element in the column
            A[k,:], A[selected,:] = A[selected,:],A[k,:] # We swap the rows
            pivot = A[k,k]; # And make a new pivot
            #
        else # Nope, everything was fine :D
            pivot = A[k,k]; #The pivot is the k,k since we want an identity matrix
        end
        #
        A[k,:] = A[k,:]./pivot #We scale the whole current row (Pivot row)
        for j in (k+1):numRow # Last, we go from the current+1 row, up until the bottom
            A[j,:] += -A[j,k].*A[k,:] # We multiply the pivot row by the negative of the current row, and add it to the current row
        end
    end
    #
    return A
end

function lgsolve(A)
    # First, we make sure to convert all to a Float64 array
    A = convert(Array{Float64,2},A)
    numRow, numCol = size(A); # Then retrieve some important variable names
    #x=zeros(numRow,1); 
    x = Array{Float64}(undef,numRow)
    #
    for k in 1:numCol-1
        # Is the pivot zero?
        if A[k,k] == 0.0 # Yes
            nonzero = findall(x->x!==0.0, A[:,k]) # Find all non-zero elements in the current column
            selected = nonzero[1] # We then select the first non-zero element in the column
            A[k,:], A[selected,:] = A[selected,:],A[k,:] # We swap the rows
            pivot = A[k,k]; # And make a new pivot
            #
        else # Nope, everything was fine :D
            pivot = A[k,k];
        end
        #
        A[k,:] = A[k,:]./pivot #We scale the whole current row
        for j in (k+1):numRow # Last, we go from the current+1 row, up until the bottom
            A[j,:] += -A[j,k].*A[k,:] # We multiply the pivot row by the negative of the current row, and add it to the current row
        end
    end


# Backwards solution
    x[end] = A[end,end];
    for k in reverse(1:length(x)-1)
        S=0; #initialize
        for m = 1:length(x)
            S += A[k,m]*x[m];
        end
        x[k] = A[k,end] - S;
    end
    return x
end

function redGauss2(A)
    # First, we make sure to convert all to a Float64 array
    A = convert(Array{Float64,2},A)
    numRow, numCol = size(A); # Then retrieve some important variable names
    #
    for k in 1:numCol-1
        # Finding the maximum value to take as a pivot
        ab = abs.(A[:,k])
        pivotIndex = findall(x->x == maximum(ab),ab) #Looking for the index of the pivot

        # Taking the first maximum, just in case it repeats
        if length(pivotIndex) > 1
            pivotInd = pivotIndex[1];
        else #No type warn
            pivotInd = pivotIndex;
        end

        pivot = A[pivotInd,k];

        # Rearranging the rows
        A[k,:], A[pivotInd,:] = A[pivotInd,:],A[k,:]

        A[k,:] = A[k,:]./pivot #We scale the whole current row
        for j in (k+1):numRow # Last, we go from the current+1 row, up until the bottom
            A[j,:] += -A[j,k].*A[k,:] # We multiply the pivot row by the negative of the current row, and add it to the current row
        end
    end
    return A
end


function lgsolve2(A)
    # First, we make sure to convert all to a Float64 array
    A = convert(Array{Float64,2},A);
    numRow, numCol = size(A); # Then retrieve some important variable names
    x=zeros(numRow,1); 
    #
    for k in 1:numCol-1
        # Finding the maximum value to take as a pivot
        ab = abs.(A[:,k])
        pivotIndex = findall(x->x == maximum(ab),ab) #Looking for the index of the pivot

        # Taking the first maximum, just in case it repeats
        if length(pivotIndex) > 1
            pivotInd = pivotIndex[1];
        else #No type warn
            pivotInd = pivotIndex;
        end

        pivot = A[pivotInd,k];

        # Rearranging the rows
        A[k,:], A[pivotInd,:] = A[pivotInd,:],A[k,:]

        A[k,:] = A[k,:]./pivot #We scale the whole current row
        for j in (k+1):numRow # Last, we go from the current+1 row, up until the bottom
            A[j,:] += -A[j,k].*A[k,:] # We multiply the pivot row by the negative of the current row, and add it to the current row
        end
    end

    # Backwards solution
    x[end] = A[end,end];
    for k in reverse(1:length(x)-1)
        S=0; #initialize
        for m = 1:length(x)
            S += A[k,m]*x[m];
        end
        x[k] = A[k,end] - S;
    end
    return x
end
