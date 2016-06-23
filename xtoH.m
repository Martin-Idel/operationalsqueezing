%Function takes size and parameters for input and outputs our H.
function H=xtoH(n,x)
    C=triu(ones(n,n));              %Init. upper triangular matrix.
    A=C;
    B=C;
    A(C==1)=x(1:n*(n+1)/2);         %Fill upper part with vector.
    B(C==1)=x(n*(n+1)/2+1:n*(n+1));  
    A=A+transpose(triu(A,1));       %Fill lower part with transpose.
    B=B+transpose(triu(B,1));
    H=[A B; B -A];                  %Fill H.
end