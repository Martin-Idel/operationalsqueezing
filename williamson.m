%Given a positive definite matrix g, computes the williamson normal form.
%Output: S^{-T}S^{-1} if S^{-T}LS^{-1}=g. One can easily change this to
%output S,L.
function [s0,s,l] = williamson(g)
    n=int16(length(g)/2);
    J=[zeros(n,n) eye(n); -eye(n) zeros(n,n)]; 
    %Create the antisymmetric matrix g^{-1/2}*J*g^{-1/2} and 
    %block-diagonalize it (using Schur).
    middle=eye(2*n)/sqrtm(g)*J/sqrtm(g);
    [Q,T]=schur(middle);
    
    %Switch basis to block basis of J.
    P=zeros(2*n,2*n);
    for i=1:n
       P(2*i-1,i)=1;
       P(2*i,n+i)=1;
    end
    nQ=Q*P;
    nT=transpose(P)*T*P;
    for i=1:n
        if nT(i,n+i)<0
            nT(i,n+i)=-nT(i,n+i);
            nT(n+i,i)=-nT(n+i,i);
            a=nQ(:,i);
            nQ(:,i)=nQ(:,n+i);
            nQ(:,n+i)=a;
        end
    end
    
    %Create the Williamson normal form (D is inverse).
    D=[nT(1:n,n+1:2*n) zeros(n,n); zeros(n,n) nT(1:n,n+1:2*n)];

    %save 1/D as normal form, s as symplectic matrix.
    l=eye(2*n)/D;
    s=eye(2*n)/sqrtm(g)*nQ/sqrtm(D);
    s0=eye(2*n)/transpose(s)*eye(2*n)/s;
end