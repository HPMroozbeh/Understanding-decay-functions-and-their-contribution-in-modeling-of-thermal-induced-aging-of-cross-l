clear;
clc;
format long
%% Network peobability fitted parameters
%%% reletive chain length
%R_s is the end to end distance of chain length of unaged material
Rbar_s = 4;

%R_s is the end to end distance of chain length of totally aged material
Rbar_b = 3.8;

% an indicator for crosslink density
% In my paper I have  p_s which is (1-q_s). 
% q boundary is 0< q_s or q_b <1
%since totally aged is usually harder than unaged, q_b usually is smalle
%than q_s. In general higher q means softer material.
q_s = 1.0;
q_b = 0.97;

% we won't need it in this section but can be an indicator of lost chain
% I did not use these, but they can be used like paper of Vahid as an
% indicator of lost chain
beta_s = 1;
beta_b = 1;

% nue is the exact same parameter as the one in 2009 paper of Dr. dargazany
nue = 1.02;

% n_max is the highest number of chains. usually I put it a large number
% and it is not going to affect our result. Becareful large number is
% reletive to R values. if you put R_s as 60, then n_max should be at much
% higher than 100. It should be at least higher than lambda*R* nue.
n_max = 100;

%this parameter only shift verything up or down without changing their
%shape.
N_s00= 19.04*10^6;

% calculating the area under probability functions for network s & b. I
% need them for calculating the total number of chain segments in the
% network.
area_s = ProbIntegral(nue*Rbar_s,n_max,Rbar_s,q_s);
area_b = ProbIntegral(nue*Rbar_b,n_max,Rbar_b,q_b);

% making sure that there is no segment generation. Total number of segments
% of aged network is equal to total number of segments of unaged network
% using these two equations:
C1=((FirstMoment(nue*Rbar_s,n_max,Rbar_s,q_s)*area_b)/(FirstMoment(nue*Rbar_b,n_max,Rbar_b,q_b)*area_s));
N_b0inf = C1*N_s00;

% scale parameter between chi and lambda. It is the same parameter as in
% 2009 paper of Dr. dargazany.
C = 0.2;
phi= 1/3;

% change Chi based on the order of macro stretch that you whant to present.
Chi = 6;
Lambda = (Chi - C^(phi))/(1 - C^(phi));

%%
%% decay function parameters
% A_1 and A_2 has the relation of A_2 = 1-A_1. But, I wrote it relatively
% in the code, so you don't have to worry about this relation. Therefore,
% if you are using generic algorithm you can let them be any number in
% general. If you need a boundary keep them between 0 and 1 as their
% relative value is what is important.
A1 = 0.4;
A2 = 0.6;
% In my paper I used tau instead of these which is 1/T. But it does the
% same thing. It only helps us scale our time frame.
T1 = 1*10^2;
T2 = 2*10^2;
%Constants
% E is the activation energy parameter. If you read my paper in polymer
% degradation journal you can know all about it. Their boundary is between
% 20000 - 130000 usually.
E1 = 35000;
E2 = 28500;
% R is constant. You don't need to worry about it.
R = 8.314;
%reference temperature. It is the temperature you use as reference. I used
%60 but in general it is better to use highest temperature as reference.
%You choose it once and keep it constant for the whole process.
temp = 333;
% desired temperature. You change this temperature to the temperature that
% you are using the data for.
Des_temp = 353;
%Time-temperature superposition. To better understand these two read either
%my 2020 paper or celina's paper
af1 = exp((E1/R)*((1/temp)-(1/Des_temp)));
af2 = exp((E2/R)*((1/temp)-(1/Des_temp)));
% aging time. It is the experimental test time in second.
t=0*24*3600;
%% Calculating Stress
% In this section you don't need to change anything in general , this the
% core of code that compute the stress based on what parameters you gave
% it.
% These are all the 21 constants for directions of microsphere.
% 0.5 is the step size. If you don't need accuracy, you can increase it for
% faster runs. 
Step_size = 0.5;

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
syms Lambda_s Lambda_max Lambda_y
% Deformation gradient tensor
F_s=[Lambda_s,0,0;0,1/sqrt(Lambda_s),0;0,0,1/sqrt(Lambda_s)];
Lambda_D=sym(zeros(1,21));
Lambda_Dmax=zeros(1,21);
for i=1:21
    D=[xc(i),yc(i),zc(i)];
    Lambda_D(i)=sqrt(D*transpose(F_s)*F_s*transpose(D));
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

if Lambda_s>=1 && Lambda_s<=Lambda
    while Lambda_s<Lambda
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
     
     Lambda_s=Lambda_s+Step_size;
     P_s(:,:,q)=0;
     P_b(:,:,q)=0;
    end
end

% Changing stress to MPa
for i=1:(q-1)
stress(i)=T(1,1,i)/10^6;
end
Micro_strain(1) = 0;
Micro_strain(2) = Step_size;
for i = 3:(cunt-1)
Micro_strain(i) = Micro_strain(i-1)+Step_size;
end
%Micro_strain = 0:Step_size:11;
Macro_strain = (1-(0.2^(1/3))).*(Micro_strain+1)+(0.2^(1/3)-1);
%% Plot Section
% Experimental data point. Becareful, the temperature you import should be
% the same as the desired temperature. The strucrure of data is as follow:
% First column is strain in %, second column is unaged stress, 3rd column
% is 1 day aging. 4th is 10 days and 5th is 30 days. Therefore, just change
% the "?" in num(:,?), to calculate the respecting stress with respect
% to time. For instance if t is t = 10*24*3600 the stress is num(:,4).
num = xlsread('BPU_80.xlsx');
strain_exp=num(:,1)./100;
stress_exp=num(:,2);
% color code for plotting
orange = [0.93 0.69 0.13];
purple = [0.71 0.27 1];
green = [ 0.39 0.83 0.07];
plot(Macro_strain,stress,'-.','Color',green,'LineWidth',3)
hold on
plot(strain_exp,stress_exp,'o','Color',green,'LineWidth',3,'MarkerSize',9)
xlim([0 5])
ylim([0 10])
title('\theta = 80')
xlabel('Strain') 
ylabel('Stress [MPa]')
%legend({'BPU Aged for 0 days'},'Location','northwest')

% error calculation
square_er = 0;
max_error = 0;

i =2;


for j= 1:(q-1)
if  strain_exp(i) <= Macro_strain(j) && i <= length(strain_exp)-2
    
    square_er = square_er + (((stress(j-1)+(((strain_exp(i)-Macro_strain(j-1))/(Macro_strain(j)-Macro_strain(j-1)))*(stress(j)-stress(j-1))))-stress_exp(i))/stress_exp(i))^2;
    max_error = max(max_error , sqrt(square_er));
    i = i+1;
end
end

Error = sqrt(square_er/length(stress_exp))*100
max_error*100