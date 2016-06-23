%This function computes a subgradient of the objective function. The
%formula is given in the main text.
function grad = objectivegrad(x)
    global F
    n = int16(1/2*(sqrt(length(x)*4+1)-1));
    H=xtoH(n,x);
    [V,D] = eig(H);
    [D,I] = sort(diag(D));
    if D(1,1)<-1 
        grad=zeros(n*(n+1),1); %outside the set of H, we just set the gradient to zero (this case should not happen).
    else
        V = V(:,I); %eigenfunctions are normalized.
        %idea: compute first V^T*F, then reshape matrix to only compute the
        %diagonal elements of V^T*F(i)*V corresponding to the elements
        %v^T_j*F(i)*v_j.
        part1=reshape(V(:,n+1:2*n)'*F,[n 2*n n*(n+1)]); %make it three dimensional to work with it
        part2=reshape(permute(part1,[2 1 3]),[1 2*n*n n*(n+1)]);
        part3=reshape(permute(part2,[1 3 2]),[n*(n+1) 2*n*n]);
        %we will compute sum_j v^T_j*F(i)*v_j directly, therefore we need
        %to first include the scaling 1/((lambda_j+1)(lambda_j-1))
        Dneu=ones(n,1)./((1+D(n+1:2*n,1)).*(1-D(n+1:2*n,1)));
        Vmult=bsxfun(@times,V(:,n+1:2*n),Dneu'); %apparently faster than reshape, multiplies columns by scalar in vector.
        Vnew=reshape(Vmult,[2*n*n 1]); 
        grad=part3*Vnew;
    end
end