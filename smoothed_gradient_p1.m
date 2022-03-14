clear;
Y1=[-0.41258837, 0.02793465, -0.30337507, 0.01861933, 0.09304836, -0.17540167, 0.72967253, 0.19329906, -0.21152674, -0.01773728];
Y2=[0.82963183, 0.84673504, 0.92476774, 0.88551538, 0.85099446, 0.90818556, 1.33082285, 1.02993049, 1.56206885, 0.74440564];
m1=10; m2=10; M=10; lambda=5; N=2; max_ite=50; eps=1e-2; theta=3/4;
b1=ones(m1,1)/10;
b2=ones(m2,1)/10;
X=(Y1+Y2)/2;
a=ones(M,1)/M;
changea=99;
changex=99;
count=0;
tic;
while changex>eps
    t=0;
    ahat=ones(M,1)/M;
    atilde=ones(M,1)/M;
    M1=makeM(X,Y1);
    K1=exp(-lambda*M1);
    U1=K1.*M1;
    M2=makeM(X,Y2);
    K2=exp(-lambda*M2);
    U2=K2.*M2;
    for ite=1:100
        beta=(t+1)/2;
        anew=(1-1/beta)*ahat+1/beta*atilde;
        if (norm(anew-a)/norm(a)<eps && ite>2)
            break;
        end
        a=anew;
        [~,~,u1,v1]=sinkhornTransport(a,b1,K1,U1,lambda);
        [~,~,u2,v2]=sinkhornTransport(a,b2,K2,U2,lambda);
        alpha1=-1/lambda*log(u1)+sum(log(u1))/lambda/m1*ones(M,1);
        alpha2=-1/lambda*log(u2)+sum(log(u2))/lambda/m2*ones(M,1);
        alpha=(alpha1+alpha2)/N;
        atilde=softmax(beta*alpha);
        ahat=(1-1/beta)*ahat+1/beta*atilde;
        t=t+1;
    end
    
    T1=diag(u1)*K1*diag(v1);
    T2=diag(u2)*K2*diag(v2);
    newX=(1-theta)*X+theta/2*(Y1*T1'+Y2*T2')*diag(1./a);
    changex=norm(newX-X)/norm(X);
    X=newX;
    count=count+1;
    disp(changex);
    %lambda=lambda*1.2;
end
figure();
bar(X,a,10);
xlim([-0.7,1.7]);
disp(sum(T1.*M1,'all')+sum(T2.*M2,'all'))
disp(toc)
    

