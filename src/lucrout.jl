function lucrout(A)
    A = convert(Array{Float64,2},A);
    numRow, numCol = size(A);
    if numRow != numCol
        return println("Non-square matrix")
    end

    L, U = zeros(numRow,numCol), zeros(numRow,numCol); # initialize

    # We first put the first column of L and the ones in U
    for j in 1:numRow
        L[j,1] = A[j,1];
        U[j,j] = 1.0;
    end

    # We then obtain the fist row of U by dividing
    for j = 2:numRow
        U[1,j] = A[1,j]/L[1,1];
    end

    # Now we try to calculate the rest

    for i = 2:numCol
        # We calculate all the i-th column of L
        for j in i:numRow
            S = 0; #initialize
            for k in 1:i-1
                S += L[j,k]*U[k,i];
            end
            L[j,i] = A[j,i] - S;
        end
        
        #We know calculate the U-th row of U 
        for j = (i+1):numCol
            S = 0; # Initialize
            for k= 1:(i-1)
                S += L[i,k]*U[k,j];
            end
            U[i,j] = 1/L[i,i]*(A[i,j]-S);
        end

    end
        return L,U
end

