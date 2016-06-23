%optimization function: construct matrices A+iB from vector x given and
%calculate the product of the cayley transform of its entries. Note that
%this might be infinite quite often.
function f = objective(x)
    global convfail
    n = int16(1/2*(sqrt(length(x)*4+1)-1));
    H=xtoH(n,x);
    D = sort(eig(H));
    f=max(-(D(1,1)+1)*convfail,1/2*log(prod((1+D(n+1:2*n,1))./(1-D(n+1:2*n,1)))));
end