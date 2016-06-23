%This function computes the gradient at the maximum residual value.
%If this is positive, the matrix violates the constraints.
function f=maxresidualgrad(x)
    global g F
    n = int16(1/2*(sqrt(length(x)*4+1)-1));
    H = xtoH(n,x);            %get H from x
    [V1,D1]=eig(g-H);
    [V2,D2]=eig(eye(2*n)+H);
    if (min(diag(D1))<0) && (min(diag(D1))<=min(diag(D2)))
        [~,I1] = sort(diag(D1));
        V1 = V1(:,I1);
        %element i of gradient given by v_1^T*F(i)*v_1 (via reshape)
        fint1=reshape(V1(:,1)'*F,[1 2*n n*(n+1)]);
        fint2=reshape(permute(fint1,[1 3 2]),[n*(n+1) 2*n]);
        f=fint2*V1(:,1);
    elseif (min(diag(D2))<0) && (min(diag(D2))<min(diag(D1)))
        [~,I2] = sort(diag(D1));
        V2 = V2(:,I2);
        %element i of gradient given by v_1^T*F(i)*v_1 (via reshape)
        fint1=reshape(V2(:,1)'*F,[1 2*n n*(n+1)]);
        fint2=reshape(permute(fint1,[1 3 2]),[n*(n+1) 2*n]);
        f=fint2*V2(:,1);
    else
        f=zeros(n*(n+1));
    end
end