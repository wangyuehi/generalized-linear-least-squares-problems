function [xhat ] = method1(y, A, B)
    [m, n] = size(A);
    
    
     [Q,R] = qr(A);

    
    R = R(1:n,1:n);
    Q1 = Q(1:m, 1:n);
    Q2 = Q(1:m,n+1:end);
    temp = Q' * B;
    [U, P] = householder (temp);
    U2 = U(n+1:end, n+1:end);
    U12 = U(1:n, n+1:end);
    z2hat = U2\(Q2'*y);
    xhat = R\(Q1' * y  - U12 * z2hat );
end

function [U,P] = householder(A)
% apply householder from right hand side
% to transform A into upper triangular
[m, n] = size(A);
U = A;
P = eye(n,n);
for i = m: -1: 2
    x = (U(i, 1: i))';
    alpha = norm(x);
    s = - sign(A(i,i));
    e1 = zeros(i,1); e1(i) = 1;
    u = x - alpha * s * e1;
    v = u / norm(u);
    I = eye(i);
    H = I - 2 * (v * v'); 
    Hb = eye(n,n);
    Hb(1:i, 1:i) = H;
    U = U * Hb;
    P = P * Hb;
    
end

end
