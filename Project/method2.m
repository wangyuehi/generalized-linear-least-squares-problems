function [xhat] = method2(y,A,B)
    %%%%%
    % y: m x 1
    % A: m x n
    % B: m x m, lower triangular form
    %%%%%
    % reference C.C. Paige. Fast numerically stable computations for generalized linear least squares problems. SIAM J. Num. Anal., 16:165-171, 1979
    % Perform the first stage described in reference
    % Annihilate superdiagonal of [y,A] while preserve lower triangular structure of B
    %%%%%
    
    [m,n] = size(A);
    yA = [y,A];
    [yAm, yAn] = size(yA);
    
    %yAB = [y,A,B];
    %IP = eye(m+n+1, m+n+1);
    Q = eye(m,m);
    P = eye(m,m);
    %{
    givens` Givens rotation matrix.
    G = givens(x,y) returns the complex Givens rotation matrix
 
        | c       s |                  | x |     | r |
    G = |           |   such that  G * |   |  =  |   |
        |-conj(s) c |                  | y |     | 0 |
 
    where c is real, s is complex, and c^2 + |s|^2 = 1.
 
    
    %}
    % apply givens rotation, start from the upper right corner
    
    % a matrix to translate [r, 0] to [0, r]
    trans = [0,1;1,0];
    for i = 1 : n
        GivensColIndex = n+1;
        GivensRowIndex = i;
        while GivensRowIndex > 0 && GivensColIndex > 0
            g = givens(yA(GivensRowIndex, GivensColIndex),yA(GivensRowIndex+1, GivensColIndex) );
            g = trans*g;
            
            G = eye(m,m);
            G(GivensRowIndex,GivensRowIndex) = g(1,1);
            G(GivensRowIndex,GivensRowIndex+1) = g(1,2);
            G(GivensRowIndex+1,GivensRowIndex) = g(2,1);
            G(GivensRowIndex+1,GivensRowIndex+1) = g(2,2);
            
            Q = G * Q;
            
            % apply G to [y,A,B];
            yA = G * yA;
            B = G * B;
            
            % new non-zero introduced at current row, super diagonal
            % position
            
            % apply givens from right
            
             g = givens(B(GivensRowIndex, GivensRowIndex),B(GivensRowIndex, GivensRowIndex+1) );
             g = g';
             
             G2 = eye(m,m);
             G2(GivensRowIndex,GivensRowIndex) = g(1,1);
             G2(GivensRowIndex,GivensRowIndex+1) = g(1,2);
             G2(GivensRowIndex+1,GivensRowIndex) = g(2,1);
             G2(GivensRowIndex+1,GivensRowIndex+1) = g(2,2);
            
            
             P = P * G2;
            B = B * G2;
            
            GivensRowIndex = GivensRowIndex-1;
            GivensColIndex = GivensColIndex-1;
            
        end
         
    end
    
    
    
    for i = n+1 : m-1
        GivensColIndex = n+1;
        GivensRowIndex = i;
        while GivensRowIndex > 0 && GivensColIndex > 0
            g = givens(yA(GivensRowIndex, GivensColIndex),yA(GivensRowIndex+1, GivensColIndex) );
            g = trans*g;
            
            G = eye(m,m);
            G(GivensRowIndex,GivensRowIndex) = g(1,1);
            G(GivensRowIndex,GivensRowIndex+1) = g(1,2);
            G(GivensRowIndex+1,GivensRowIndex) = g(2,1);
            G(GivensRowIndex+1,GivensRowIndex+1) = g(2,2);
            
            Q = G * Q;
            
            % apply G to [y,A,B];
            yA = G * yA;
            B = G * B;
            
            % new non-zero introduced at current row, super diagonal
            % position
            
            % apply givens from right
            
             g = givens(B(GivensRowIndex, GivensRowIndex),B(GivensRowIndex, GivensRowIndex+1) );
             g = g';
             
             G2 = eye(m,m);
             G2(GivensRowIndex,GivensRowIndex) = g(1,1);
             G2(GivensRowIndex,GivensRowIndex+1) = g(1,2);
             G2(GivensRowIndex+1,GivensRowIndex) = g(2,1);
             G2(GivensRowIndex+1,GivensRowIndex+1) = g(2,2);
            
            
             P = P * G2;
            B = B * G2;
            
            GivensRowIndex = GivensRowIndex-1;
            GivensColIndex = GivensColIndex-1;
            
            
        end
        
%         eta = yA(m-n,1);
%         w = yA(m-n+1 : m,1);
%         lambda = B(m-n,m-n);
%         
%         L = yA(m-n+1:m,2:end);
%         
%         z = zeros(m,1);
%         z(1:m-n-1) = 0;
%         z(m-n) = eta / lambda;
%         z(m-n+1:m) = 0;
        eta = yA(m-n, 1);
        z = yA((m-n)+1:m,1);
        Rt = yA(m-n+1:m, 2:n+1);
        
        rho = B(m-n,m-n);
        r = B(m-n+1:m,m-n);
        
        miuhat = eta/rho;
        Rtxhat = z - miuhat * r;
        
        xhat = (Rt) \ Rtxhat;
        
        
         
    end
end