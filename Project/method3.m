function [xhat] = method3(y,A,B)
    
% compute QR of A by householder 
% [Q,R] = qr(A);
Sigma = B * B';
[m, n] = size(A);


   [Q,R] = qr(A);
    R = R(1:n,1:n);
    Q1 = Q(1:m, 1:n);
    Q2 = Q(1:m,n+1:end);
% apply Q to sigma
   Sbar = Q' * Sigma * Q;

% partial Cholesky
U = partialChol(Sbar, n+1);
% [U] = chol(Sbar);

% B = Q*U;
  U2 = U(n+1:end, n+1:end);

   U12 = U(1:n, n+1:end);
    z2hat = U2\(Q2'*y);
    xhat = R\(Q1' * y  - U12 * z2hat );

end