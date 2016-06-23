%This function constructs an (approximate) solution to the minimal
%squeezing problem. It also outputs all bounds known.
%options: options(1)==1 --> calculates bounds. options(1)==2 no bounds.
%options(2)==1 numerical gradients, options(2)==2 numerical objective,
%options(2)==3 analytical gradients
function [mi,cs,sympap,bds,preperror]=minimum(dim,gamma,solvoptions,options)
    global n g F convfail
    %Define n and predefine convfail (later to be defined)
    n=dim;
    convfail=1;
    
    if(options(2)==2 || options(2)==3)
        %Produce matrices to calculate the gradient.
        Fbegin(:,:,:)=zeros(2*n,2*n,n*(n+1));
        for i=1:n*(n+1)
            xi=zeros(n*(n+1));
            xi(i)=1;
            Fbegin(:,:,i)=xtoH(n,xi);
        end
        F=reshape(Fbegin,2*n,2*n*(n*(n+1))); 
    end

    %Check positivity and produce Cayley-transformed g:
    diff=eig(gamma-sqrtm(transpose(gamma)*gamma));
    for i=1:2*n
        if diff(i)>1.e-8
            error('Error: The provided matrix is not positive definite');
        end
    end
    I=eye(2*n);
    g=(gamma-I)/(gamma+I);

    %produce starting point vector:
    [symp,~,l]=williamson(gamma);    
    s0=(symp-I)/(symp+I);
    x0=Htox(n,s0);

    %check matrix: it must be a valid covariance matrix and should not lie
    %on the boundary.
    warn=0;
    for i=1:n
        if l(i,i)<1-1.e-10
            error('Starting point no valid covariance matrix');
        end
        if l(i,i)>1-1.e-10 && l(i,i)<1+1.e-10
            warn=1;
        end
    end
    
    %Set the value such that the function appears convex.
    convfail=10000000*objective(x0);
    
    %Solve
    if options(2)==1
        [y,z]=solvopt(x0,'objective',[],solvoptions,'maxresidual');
    elseif options(2)==2
        [y,z]=solvopt(x0,'objective','objectivegrad',solvoptions,'maxresidual');
    else
        [y,z]=solvopt(x0,'objective','objectivegrad',solvoptions,'maxresidual','maxresidualgrad');
    end

    %rewrite the output vector as matrix.
    H=xtoH(n,y(:,1));
    G=(eye(2*n)+H)\(eye(2*n)-H);
    [~,S,~]=williamson(G);

    %this is the amount of necessary squeezing and the matrix.
    cs=maxresidual(y(:,1)); %check that matrix satisfies constraints
    mi=z; %matrix result
    sympap=S; 
    
    %Calculate the preparation error.
    correction=min(eig(-eye(2*n)+eye(2*n)/S*gamma/transpose(S)));
    if correction>0
        correction=0;
    end
    preperror=max(eig(eye(2*n)/transpose(sympap(:,:,1))/sympap(:,:,1))*abs(correction));
    
    if options(1)==1
        bds=bounds(n,gamma); %bounds
    else
        bds=[0 0 0 0];
    end
    
    %Write output to console. Maximum constraint violation: 10^{-7}
    fprintf('------------------------------------------------------------------------\n');
    if warn==1
        fprintf('Warning: starting point lies on the boundary!\n');
    end
    if preperror>solvoptions(6)*10^8
        fprintf('Error: Matrix cannot be prepared in the way SolvOpt implies.\n');
        fprintf('Lower bounds: %6.5f, %6.5f, %6.5f\n',bds(1),bds(2),bds(3));
        fprintf('Upper bounds: %6.5f, %6.5f\n',z,bds(4));
        fprintf('Preparation error: %6.5f\n',preperror);
    elseif preperror<=solvoptions(6)*10^8 && preperror>solvoptions(6)*10^3
        fprintf('Warning: Minimum preparation not very accurate.\n');
        fprintf('Minimum: %6.5f \n',z);
        fprintf('Lower bounds: %6.5f, %6.5f, %6.5f \n',bds(1),bds(2),bds(3));
        fprintf('Upper bounds: %6.5f \n',bds(4));
        fprintf('Preparation error: %6.5f\n',preperror);
    else
        fprintf('Minimum: %6.5f \n',z);
        fprintf('Lower bounds: %6.5f, %6.5f, %6.5f \n',bds(1),bds(2),bds(3));
        fprintf('Upper bounds: %6.5f \n',bds(4));        
    end
    fprintf('------------------------------------------------------------------------\n');
end