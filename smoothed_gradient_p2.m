Y1=[ 2.13647389, -2.64470809,  5.84631894, -6.98020947, -2.94254153, 0.04958081, -6.58950368,...
     0.95205625, -5.73286034, -0.09772471, -1.41441751,  4.08656681, -3.81226385,  0.27255582,...
     -5.02069786, -3.17954019,  7.19494216, -8.38522149, -9.47602305,  4.71295821];
Y2=[23.77506235, 25.90793152, 14.44796339, 13.38625086, 30.95782642,...
    27.63501945, 21.09838795, 20.29635571, 22.97877629,  7.55646666,...
    10.80359335, 26.43684217, 21.24387921, 36.68363111,  8.70643644,...
    37.63070613,  6.77734705, 17.01520291, 10.23293397, 13.87057592,...
    28.4110526 ,  3.84849208, 36.09820298, 31.27980933, 34.43577525,...
    21.61291983,  8.28579134, 17.27605099, 13.98069211, 18.94498462];
Y3=[51.01611139, 61.37384657, 54.58824101, 54.24744202, 57.46986682,...
       51.05774993, 44.22655963, 47.22199274, 47.96861929, 48.65183471];
m1=20; m2=30; m3=10; M=60; lambda=1e-1; N=3; max_ite=100; eps=1e-2; theta=1/2;
b1=ones(m1,1)/m1;
b2=ones(m2,1)/m2;
b3=ones(m3,1)/m3;
X=[Y1,Y2,Y3];
a=ones(M,1)/M;
changea=99;
changex=99;
count=0;
tic;
while changex>eps && count<50
    t=0;
    ahat=ones(M,1)/M;
    atilde=ones(M,1)/M;
    M1=makeM(X,Y1);
    K1=exp(-lambda*M1);
    U1=K1.*M1;
    M2=makeM(X,Y2);
    K2=exp(-lambda*M2);
    U2=K2.*M2;
    M3=makeM(X,Y3);
    K3=exp(-lambda*M3);
    U3=K3.*M3;
    for ite=1:max_ite
        beta=(t+1)/2;
        anew=(1-1/beta)*ahat+1/beta*atilde;
        if (norm(anew-a)/norm(a)<eps && ite>2)
            break;
        end
        a=anew;
        [~,~,u1,v1]=sinkhornTransport(a,b1,K1,U1,lambda);
        [~,~,u2,v2]=sinkhornTransport(a,b2,K2,U2,lambda);
        [~,~,u3,v3]=sinkhornTransport(a,b3,K3,U3,lambda);
        alpha1=-1/lambda*log(u1+1e-7)+sum(log(u1+1e-7))/lambda/m1*ones(M,1);
        alpha2=-1/lambda*log(u2+1e-7)+sum(log(u2+1e-7))/lambda/m2*ones(M,1);
        alpha3=-1/lambda*log(u3+1e-7)+sum(log(u3+1e-7))/lambda/m3*ones(M,1);
        alpha=(alpha1+alpha2+alpha3)/N;
        atilde=softmax(beta*alpha);
        ahat=(1-1/beta)*ahat+1/beta*atilde;
        t=t+1;
    end
    T1=diag(u1)*K1*diag(v1);
    T2=diag(u2)*K2*diag(v2);
    T3=diag(u3)*K3*diag(v3);
    newX=(1-theta)*X+theta/3*(Y1*T1'+Y2*T2'+Y3*T3')*diag(1./a);
    changex=norm(newX-X)/norm(X);
    X=newX;
    count=count+1;
    disp(changex);
    lambda=lambda*1.2;
end
figure();
bar(X,a,10);
%xlim([-0.7,1.7]);
xlim([10,40]);
disp(sum(T1.*M1,'all')+sum(T2.*M2,'all')+sum(T3.*M3,'all'))
disp(toc)
    

