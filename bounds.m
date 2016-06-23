%This function calculates a number of bounds. The first bound is the best
%spectral lower bound, the sdp-lower bound, the third being the lower bound
%achieved from assuming superadditivity and the fourth is the best spectral
%upper bound.
function [v]=bounds(n,gamma)
    %Define the symplectic form and make sure that gamma is positive def.
    J=zeros(2*n,2*n);
    J(1:n,n+1:2*n)=eye(n,n);
    J(n+1:2*n,1:n)=-eye(n,n);
    g=logm(sqrtm(transpose(gamma)*gamma));
    
    %Call cvx: constraints are that log(gamma)>=H and H in pi(n)
    cvx_begin sdp
        variable H(2*n,2*n) symmetric
        minimize(norm(H,1))
        H<=g;
        1i*(H*J+J*H)<=0;
        1i*(H*J+J*H)>=0;
    cvx_end
    sdp=1/4*cvx_optval;
    
    %compute spectral bounds.
    v=sort(eig(gamma));
    lower=1;
    for i=1:n
        if v(i)<1
            lower=lower*1/v(i);
        end
    end
    l=1/2*log(lower);
    u=1/2*log(prod(v(n+1:2*n)));
    
    %compute lower bound from assuming superadditivity using the analytic
    %calculations from the 2x2-case.
    sm=zeros(n,1);
    for i=1:n
        sm(i,1)=1/2*log(1/min(svd(gamma([i n+i],[i n+i]))));
    end
    small=sum(sm);
    v=[l,sdp,small,u];
end