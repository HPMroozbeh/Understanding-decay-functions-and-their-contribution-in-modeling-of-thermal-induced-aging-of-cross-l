clear;
clc;
close all
format long
% Network peobability fitted parameters
% reletive chain length
Rbar_s = 4;
Rbar_b = 3.8;
% an indicator for crosslink density
q_s = 1.0;
q_b = 0.97;
% we won't need it in this section but can be an indicator of lost chain
beta_s = 1;
beta_b = 1;
nue = 1.002;
n_max = 100;
N_s00= 19.04*10^6;
% calculating the area under probability functions for network s & b
area_s = ProbIntegral(nue*Rbar_s,n_max,Rbar_s,q_s);
area_b = ProbIntegral(nue*Rbar_b,n_max,Rbar_b,q_b);
% making sure that there is nosegment generation
C1=((FirstMoment(nue*Rbar_s,n_max,Rbar_s,q_s)*area_b)/(FirstMoment(nue*Rbar_b,n_max,Rbar_b,q_b)*area_s));
N_b0inf = C1*N_s00;
% scale parameter between chi and lambda
C=0.2;
% decay function parameters
A1 = 0.4;
A2 = 0.6;
T1 = 1*10^2;
T2 = 2*10^2;
%Constants
E1 = 35000;
E2 = 28500;
R = 8.314;
%reference temperature
temp = 333;
% desired temperature
Des_temp = 368;
%Time-temperature superposition
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
af2 = exp((E2/R)*((1/temp)-(1/Des_temp)));
% aging time
t=30*24*3600;
% micro_sphere
xc=[1,0,0,0.7071067812,0.7071067812,0.7071067812,0.7071067812,0,0'...
    ,0.3879073041,0.3879073041,0.3879073041,0.3879073041,0.3879073041'...
    ,0.3879073041,0.3879073041,0.3879073041,0.8360955967,0.8360955967'...
    ,0.8360955967,0.8360955967];
yc=[0,1,0,0.7071067812,-0.7071067812,0,0,0.7071067812,0.7071067812'...
    ,0.3879073041,0.3879073041,-0.3879073041,-0.3879073041,0.8360955967'...
    ,0.8360955967,-0.8360955967,-0.8360955967,0.3879073041'...
    ,0.3879073041,-0.3879073041,-0.3879073041];
zc=[0,0,1,0,0,0.7071067812,-0.7071067812,0.7071067812,-0.7071067812'...
    ,0.8360955967,-0.8360955967,0.8360955967,-0.8360955967,0.3879073041'...
    ,-0.3879073041,0.3879073041,-0.3879073041,0.3879073041,-0.3879073041'...
    ,0.3879073041,-0.3879073041];
wc=[0.0265214244,0.0265214244,0.0265214244,0.0199301476,0.0199301476'...
    ,0.0199301476,0.0199301476,0.0199301476,0.0199301476,0.0250712367'...
    ,0.0250712367,0.0250712367,0.0250712367,0.0250712367,0.0250712367'...
    ,0.0250712367,0.0250712367,0.0250712367,0.0250712367,0.0250712367'...
    0.0250712367];
syms Lambda_s Lambda_max 
F_s=[Lambda_s,0,0;0,1/sqrt(Lambda_s),0;0,0,1/sqrt(Lambda_s)];
Lambda_D=sym(zeros(1,21));
Lambda_Dmax=zeros(1,21);
Lambda_maxdir=sym(zeros(1,21));
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_D(i)=sqrt(D*transpose(F_s)*F_s*transpose(D));
end
F_smax=[Lambda_max,0,0;0,1/sqrt(Lambda_max),0;0,0,1/sqrt(Lambda_max)];
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_maxdir(i)=sqrt(D*transpose(F_smax)*F_smax*transpose(D));
end
Lambda_rel=1;
Lambda_s=1;
q=1;
P_s=zeros([3 3 1]);
T=zeros([3 3 1]);
P_pp_s=zeros(1,21);
P_b=zeros([3 3 1]);
P_pp_b=zeros(1,21);
cunt = 1;
% first loading
% 2.35 is about the same amount for chi = 1.56
if Lambda_s>=1 && Lambda_s<=2.35
    while Lambda_s<2.35
    for i=1:21
        Lambda_Dmax(i)=max(1,eval(Lambda_D(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
     q=q+1;
     cunt = cunt+1;
     Lambda_s=Lambda_s+0.01;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
    end
end

%start unloading
Lambda_max= 2.35;
Lambda_s=Lambda_s-0.01;

cuntu1 =1;

if Lambda_s>=1 && Lambda_s<=2.35
    while Lambda_s>=1
    for i=1:21  
        Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
    if T(1,1,q)<=0
        break
    end
    Lambda_s=Lambda_s-0.01;
    q=q+1;
    cuntu1 = cuntu1+1;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
    end
end

% Second loading
% 3.7 is about the same amount for chi = 2.12

Lambda_s=Lambda_s+0.01;

cuntl2 = 1;
if Lambda_s>=1 && Lambda_s<=3.7
    while Lambda_s<3.7
    for i=1:21
        Lambda_Dmax(i)=max([1,eval(Lambda_D(i)),eval(Lambda_maxdir(i))]);
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
     q=q+1;
     cuntl2 = cuntl2+1;
     Lambda_s=Lambda_s+0.01;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
    end
end

%start unloading
Lambda_max= 3.7;
Lambda_s=Lambda_s-0.01;
cuntu2 =1;
if Lambda_s>=1 && Lambda_s<=3.7
    while Lambda_s>=1
    for i=1:21  
        Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
        nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
        nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
     P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
            *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
     P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
            *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
        P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
        P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
    end
    T(:,:,q)=(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-exp(-E1/(R*temp))*af1*(t/T1))+A2*exp(-exp(-E2/(R*temp))*af2*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
    if T(1,1,q)<=0
        break
    end
    Lambda_s=Lambda_s-0.01;
    q=q+1;
    cuntu2 = cuntu2+1;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
    end
end

% start of third loading

% 4.61 is about the same amount for chi = 2.5

% Lambda_s=Lambda_s+0.01;
% cuntl3 = 1;

% if Lambda_s>=1 && Lambda_s<=4.61
%     while Lambda_s<4.61
%     for i=1:21
%         Lambda_Dmax(i)=max([1,eval(Lambda_D(i)),eval(Lambda_maxdir(i))]);
%         nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
%         nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
%      P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
%             *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
%      P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
%             *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
%         P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%         P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%     end
%     T(:,:,q)=(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
%      q=q+1;
%      cuntl3 = cuntl3+1;
%      Lambda_s=Lambda_s+0.01;
%      P_s(:,:,q)=0;
%      P_b(:,:,q)=0;
%     end
% end
% 
% %start unloading
% Lambda_max= 4.61;
% Lambda_s=Lambda_s-0.01;
% cuntu3 = 1;
% if Lambda_s>=1 && Lambda_s<=4.61
%     while Lambda_s>=1
%     for i=1:21  
%         Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
%         nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
%         nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
%      P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
%             *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
%      P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
%             *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
%         P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%         P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%     end
%     T(:,:,q)=(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
%     if T(1,1,q)<=0
%         break
%     end
%     Lambda_s=Lambda_s-0.01;
%     q=q+1;
%     cuntu3 = cuntu3+1;
%      P_s(:,:,q)=0;
%      P_b(:,:,q)=0;
%     end
% end

% start of forth loading

% 5.82 is about the same amount for chi = 3

% Lambda_s=Lambda_s+0.01;
% cuntl4 = 1;
% 
% if Lambda_s>=1 && Lambda_s<=5.82
%     while Lambda_s<5.82
%     for i=1:21
%         Lambda_Dmax(i)=max([1,eval(Lambda_D(i)),eval(Lambda_maxdir(i))]);
%         nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
%         nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
%      P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
%             *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
%      P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
%             *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
%         P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%         P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%     end
%     T(:,:,q)=(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
%      q=q+1;
%      cuntl4 = cuntl4+1;
%      Lambda_s=Lambda_s+0.01;
%      P_s(:,:,q)=0;
%      P_b(:,:,q)=0;
%     end
% end
% 
% %start unloading
% Lambda_max= 5.82;
% Lambda_s=Lambda_s-0.01;
% cuntu4 = 1;
% if Lambda_s>=1 && Lambda_s<=5.82
%     while Lambda_s>=1
%     for i=1:21  
%         Lambda_Dmax(i)=max(1,eval(Lambda_maxdir(i)));
%         nmin_s=(nue*(Lambda_Dmax(i))*Rbar_s);
%         nmin_b=(nue*(Lambda_Dmax(i))*Rbar_b);
%      P_pp_s(i)=(Rbar_s)*(N_s00/(1-C^(1/3)))*Phi(nmin_s,n_max,Rbar_s,q_s)'...
%             *NIntegral(nmin_s,n_max,Rbar_s,q_s,eval(Lambda_D(i)));
%      P_pp_b(i)=(Rbar_b)*(N_b0inf/(1-C^(1/3)))*Phi(nmin_b,n_max,Rbar_b,q_b)'...
%             *NIntegral(nmin_b,n_max,Rbar_b,q_b,eval(Lambda_D(i)));
%         P_s(:,:,q)=P_s(:,:,q)+((P_pp_s(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%         P_b(:,:,q)=P_b(:,:,q)+((P_pp_b(i))*(wc(i)/((1-C^(1/3))*eval(Lambda_D(i))+C^(1/3)))*eval(F_s)*(transpose([xc(i),yc(i),zc(i)])*[xc(i),yc(i),zc(i)]));
%     end
%     T(:,:,q)=(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))*(P_s(:,:,q)-(P_s(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1))))+((1-(((1/(A1+A2))*(A1*exp(-af*(t/T1))+A2*exp(-af*(t/T2))))))*((P_b(:,:,q)-(P_b(3,3,q)/sqrt(Lambda_s))*eval(transpose(F_s^-1)))));
%     if T(1,1,q)<=0
%         break
%     end
%     Lambda_s=Lambda_s-0.01;
%     q=q+1;
%     cuntu4 = cuntu4+1;
%      P_s(:,:,q)=0;
%      P_b(:,:,q)=0;
%     end
% end

% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end

Mic_strain(1) = 0;
Mic_strain(2) = 0.01;
for i = 3:(cunt-1)
Mic_strain(i) = Mic_strain(i-1)+0.01;
end
for i = cunt:(cunt+cuntu1-2)
Mic_strain(i) = Mic_strain(i-1)-0.01;
end
for i = (cunt+cuntu1-1):(cunt+cuntu1+cuntl2-3)
Mic_strain(i) = Mic_strain(i-1)+0.01;
end
for i = (cunt+cuntu1+cuntl2-2):(cunt+cuntu1+cuntl2+cuntu2-4)
Mic_strain(i) = Mic_strain(i-1)-0.01;
end
% for i = (cunt+cuntu1+cuntl2+cuntu2-3):(cunt+cuntu1+cuntl2+cuntu2+cuntl3-5)
% Mic_strain(i) = Mic_strain(i-1)+0.01;
% end
% for i = (cunt+cuntu1+cuntl2+cuntu2+cuntl3-4):(cunt+cuntu1+cuntl2+cuntu2+cuntl3+cuntu3-6)
% Mic_strain(i) = Mic_strain(i-1)-0.01;
% end
% for i = (cunt+cuntu1+cuntl2+cuntu2+cuntl3+cuntu3-5):(cunt+cuntu1+cuntl2+cuntu2+cuntl3+cuntu3+cuntl4-7)
% Mic_strain(i) = Mic_strain(i-1)+0.01;
% end
% for i = (cunt+cuntu1+cuntl2+cuntu2+cuntl3+cuntu3+cuntl4-6):(cunt+cuntu1+cuntl2+cuntu2+cuntl3+cuntu3+cuntl4+cuntu4-8)
% Mic_strain(i) = Mic_strain(i-1)-0.01;
% end

for i=1:(q-1)
Macro_strain(i) = (1-(0.2^(1/3))).*(Mic_strain(i)+1)+(0.2^(1/3)-1);
end
% Experimental data point
num = xlsread('BPU_95C.xlsx');
strain_exp=num(:,1)./100;
stress_exp=num(:,2);
plot(Macro_strain,stress*0.72,'-k','LineWidth',3)
hold on
plot(strain_exp,stress_exp,'d','Color','r','LineWidth',3,'MarkerSize',10)
xlim([0 1.5])
ylim([0 6])

% error calculation
% square_er = 0;
% max_error = 0;
% 
% i =2;
% 
% for j= 1:(q-1)
% if  strain_exp(i) <= Macro_strain(j) && i <= length(strain_exp)-1
%     
%     square_er = square_er + ((stress(j)-stress_exp(i))/stress_exp(i))^2;
%     max_error = max(max_error , sqrt(square_er));
%     i = i+1;
% end
% end
% 
% Error = sqrt(square_er/length(stress_exp))*100
% max_error*100