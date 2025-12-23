function [Ad, Bd] = c2d_exact(A, B, Ts)
    n = size(A,1);
    M = [A, B; zeros(size(B,2), size(A,1)+size(B,2))];
    Md = expm(M*Ts);
    Ad = Md(1:n,1:n);
    Bd = Md(1:n,n+1:end);
end
