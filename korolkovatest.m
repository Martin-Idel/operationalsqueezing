close all;
clear all;

%three mode states.
n=3; 

p=30;
m=30;

minoptions(2)=3;
minoptions(1)=1;

A=zeros(m*p,7); %Array to store: d,r,minimum,lower,upper,prep,kor,otherupper

for i=1:p
    for j=1:m
        %Create the state of Mista and Korolkova
        r=0.1+j*.05;
        d=r+i*.03;
        
        a=cosh(2*r);
        c=sinh(2*r);
        f=atan(exp(-2*r)*sinh(2*d)+sqrt(1+exp(-4*r)*sinh(2*d)*sinh(2*d)));
        x=2*sinh(2*r)/(exp(2*d)*sin(f)*sin(f)+exp(-2*d)*cos(f)*cos(f));

        gamma=zeros(2*n,2*n);
        gamma(1,1)=exp(2*d)*a+cos(f)*cos(f)*x;
        gamma(1,3)=-exp(2*d)*c+cos(f)*cos(f)*x;
        gamma(1,5)=x*sqrt(2)*cos(f);
        gamma(1,6)=gamma(1,5);
        gamma(2,2)=exp(-2*d)*a+sin(f)*sin(f)*x;
        gamma(2,4)=exp(-2*d)*c-sin(f)*sin(f)*x;
        gamma(2,5)=x*sqrt(2)*sin(f);
        gamma(2,6)=gamma(2,5);
        gamma(3,3)=gamma(1,1);
        gamma(3,5)=gamma(1,5);
        gamma(3,6)=gamma(1,5);
        gamma(4,4)=gamma(2,2);
        gamma(4,5)=-gamma(2,5);
        gamma(4,6)=-gamma(2,5);
        gamma(5,5)=4*x+1;
        gamma(5,6)=4*x;
        gamma(6,6)=gamma(5,5);

        for j1=2:2*n
            for i1=1:(j1-1)
                gamma(j1,i1)=gamma(i1,j1);
            end
        end

        gamma=sigmaToJ(gamma); %This is the final state of Mista and Korolkova.
        
        %Devine options for solvopt
        options(2)=1.e-6;
        options(3)=1.e-8;
        options(6)=1.e-8;
        %options(5)=10;

        %Call function to solve.
        [mi,cs,sympap,bds,preperror]=minimum(n,gamma,options,minoptions);
                
        %Fill the array via: d,r,minimum,lower,upper,prep,kor
        A(m*(i-1)+j,1)=d;
        A(m*(i-1)+j,2)=r;
        A(m*(i-1)+j,3)=mi;
        A(m*(i-1)+j,4)=max(bds(1:3)); %best lower bound.
        A(m*(i-1)+j,5)=bds(4);
        A(m*(i-1)+j,6)=preperror;
        A(m*(i-1)+j,7)=2*d; %Cost of Korolkova procedure.
    end
end

%Create graphics
x=[A(:,1); A(:,1); A(:,1)];
y=[A(:,2); A(:,2); A(:,2)];
z=[A(:,3); A(:,4); A(:,7)];
s=5;
C=repmat([1,2,3],numel(A(:,1)),1);
c=C(:);

plot1=scatter3(x,y,z,s,c); %Plot with lower bounds and upper bounds.
xlabel('d');
ylabel('r');

plot2=scatter3(A(:,1), A(:,2),A(:,6)); %Plot error in preparation.
xlabel('d');
ylabel('r');

%xlswrite('korolkova.xls',A); %Save array