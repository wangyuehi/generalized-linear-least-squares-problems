function [ U] = partialChol( A, stopIndex )
    % partial reverse Cholesky fectorization
    % A = U * U'
    n = size(A,1);
    U = zeros(size(A));
    % compute from the last element..
    
    for i = n:-1:stopIndex
        if i > n-1
           U(n,n) = sqrt(A(n,n));
        else
            U(i,i) = sqrt(A(i,i) - U(i,i+1:n) * U(i,i+1:n)');
        end
        for j = 1:i-1
            U(j,i) = (A(i,j) - U(i,i+1:n)*U(j,i+1:n)')/U(i,i);
        end
    end
    


end

