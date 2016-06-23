%This function computes the maximum residual value. If this is positive,
%the matrix violates the constraints.
function f=maxresidual(x)
    global g
    n = int16(1/2*(sqrt(length(x)*4+1)-1));
    H = xtoH(n,x);            %get H from x
    f=max([0,-min(eig(g-H)),-min(eig(eye(2*n)+H))]); %conditions not met if smallest eigenvalues are negative.
end