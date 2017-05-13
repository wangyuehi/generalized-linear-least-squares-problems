% B = [1,0,0,0,0,0,0,0;
%      2,3,0,0,0,0,0,0;
%      4,5,6,0,0,0,0,0;
%      7,8,9,10,0,0,0,0;
%      11,12,13,14,15,0,0,0;
%      12,13,14,15,16,17,0,0;
%      14,15,16,17,18,19,20,0;
%      10,11,12,11,12,15,12,16];

 n = 10;
 relaError = zeros(25,3);
 time  = zeros(25,3);
 
 for m = 20:20:250
 B = rand(m,m);
 ind = m / 20;
 
 cor = B * B';
 A = 10 * rand(m,n);
%  v = mvnrnd(zeros(1,m),cor);
 x = ones(n,1);
 y = A * x;
B = chol(cor);

t0 = cputime;
x1 = method1(y,A,B);
t1 = cputime - t0;

time(ind,1) = t1;
% display(['cputime of method 1 ', num2str(t1)]);

t = cputime;
x2 = method2(y,A,B);
t2 = cputime - t;
% display(['cputime of method 1 ', num2str(t2)]);
time(ind,2) = t2;


t = cputime;
x3 = method3(y,A,B);
t3 = cputime - t;
% display(['cputime of method 1 ', num2str(t3)]);
time(ind,3) = t3;

% e1 = relativeError(x, x1);
% e2 = relativeError(x, x2);
% e3 = relativeError(x, x3);
relaError(ind,1) = relativeError(x, x1);
relaError(ind,2) = relativeError(x, x2);
relaError(ind,3) = relativeError(x, x3);
% display(['relative error of method 1 ', num2str(e1)]);
% display(['relative error of method 2 ', num2str(e2)]);
% display(['relative error of method 3 ', num2str(e3)]);
 end