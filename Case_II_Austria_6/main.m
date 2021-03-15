%%%%%%%% the main file for the 2D TM (scalar wave) highly nonlinear inverse scattering problem %%%%%%%%%
% (1) This is a complete program, and one can run it to obtain the recosntruction result. 
% (2) In order to save time so that readers can get the results quickly, 
% this program only contains inversion algorithms, the synthetic data, 
% which can be generated in forward problem, has been prepared in advance.

% Copyright ? 2021, University of Electronic Science and Technology of China, Xiaoniu Zhang

clear all;close all;clc 
 
%To start off, we firstly need to initialize some problem parameters.
%Secondly, we need to prepare some materials for simulation.

% Problem definition:
% The imaging domain D, is assumed to be of size  2m¡Á2m  
% and the object of interest is located within this domain. 
% The incident field is 400 MHz, and hence the free space wavelength is lamda = 0.75 m.
% The measurements are taken on a circle of radius R_obs = 3.75 m with the centre of 
% the circle coinciding with the center of D and Ns = 40 receivers are placed along 
% the measurement domain in an equiangular manner for each illuminations.
% Total number of different illuminations are Ni = 20.

% Some materials for simulation:
% 1.The Synthetic data: <E_s_6_30db>

% 2.The matrix form of the internal radiation operator
% mapping the induced current within D to scattered fields within D: <Z>

% 3.The true profile, which is used to calculate the error: <chai0>

% 4.The scattering operator that maps the contrast sources 
% within D to the measured scattered fields: <Gs>

% 5.Incident Field: <E_inc>

% 6.The L2-norm of the incident Field and the scattered field:
% <E_p_sca_norm_2> and <E_p_inc_norm_2>

 
%%%%%%%%%%%%% 1.¡¾<E_s_6_30db>¡¿%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Synthetic data prepared in advance;
%The Austria profile with the relative permittivity being 6,SNR=30db
%is used to generate the synthetic data;
load('E_s_6_30db.mat')  
E_s = E_s_6_30db;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%% 2.¡¾<Z>¡¿ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_DOI = 2;
M = 32;
d = size_DOI/M;
A = size_DOI^2;
[X_dif,Y_dif] = meshgrid(((1-M):1:(M-1))*d,((1-M):1:(M-1))*d); 
R_2d = sqrt(X_dif.^2+Y_dif.^2); % (2M-1) x (2M-1)


imp = 120*pi;
lambda = 0.75;
k0 = 2*pi/lambda;
celldia = sqrt(d^2/pi)*2; 
cellrad = celldia/2;

Gs = -imp*pi*cellrad/2*besselj(1,k0*cellrad)*besselh(0,1,k0*R_2d); % (2M-1) x (2M-1)
Gs(M,M) = -imp*pi*cellrad/2*besselh(1,1,k0*cellrad)-1i*imp/k0; 


Z = zeros(2*M-1,2*M-1);
Z(1:M,1:M) = Gs((M):(2*M-1),(M):(2*M-1));
Z((M+1):(2*M-1),(M+1):(2*M-1)) = Gs(1:(M-1),1:(M-1));
Z(1:M,(M+1):(2*M-1)) = Gs((M):(2*M-1),1:(M-1));
Z((M+1):(2*M-1),1:M) = Gs(1:(M-1),(M):(2*M-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% 3.¡¾<chai0>¡¿%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = eye(M^2);
size_DOI = 2;
d = size_DOI/M; 
tx = ((-(M-1)/2):1:((M-1)/2))*d; 
ty = (((M-1)/2):(-1):(-(M-1)/2))*d;
[x,y] = meshgrid(tx,ty);

epr_0 = 6;
epr = ones(M,M);
epr((x-0.3).^2+(y-0.6).^2<=0.2^2) = epr_0;
epr((x+0.3).^2+(y-0.6).^2<=0.2^2) = epr_0;
epr((x).^2+(y+0.2).^2>=0.3^2 & (x).^2+(y+0.2).^2<=0.6^2) = epr_0;
chai0 = epr-ones(M,M);
chai0 = chai0(:);

WasK = zeros(M,M);
WasK((x-0.3).^2+(y-0.6).^2<=0.2^2) = 1;
WasK((x+0.3).^2+(y-0.6).^2<=0.2^2) = 1;
WasK((x).^2+(y+0.2).^2>=0.3^2 & (x).^2+(y+0.2).^2<=0.6^2) = 1;
WasK = WasK(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% 4.¡¾Gs¡¿%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns = 40;
phi_temp = linspace(0,2*pi,41);phi_temp(end) = [];
phi =phi_temp;
R_obs = 3.75; % radius of the circle formed by receiving antennas
X = R_obs*cos(phi); % 1 x Ns % x coordinates of receiving antennas
Y = R_obs*sin(phi); %
X_obs = repmat(X.',[1,M*M]); % Ns x M^2
Y_obs = repmat(Y.',[1,M*M]); % Ns x M^2
R_2d = sqrt((X_obs-repmat(reshape(x,M*M,1).',[Ns,1])).^2+ ...
    (Y_obs-repmat(reshape(y,M*M,1).',[Ns,1])).^2); % Ns x M^2
Gs = -imp*pi*cellrad/2*besselj(1,k0*cellrad)*besselh(0,1,k0*R_2d); % Ns x M^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% 5.¡¾<E_inc>¡¿%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ni = 20;
theta = pi/10:pi/10:2*pi; % angle of incidence
lambda = 0.75;
k0 = 2*pi/lambda;
E_inc = exp(1i*k0*x(:)*cos(theta(:).')+1i*k0*y(:)*sin(theta(:).')); % M^2 x Ni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% 6.¡¾E_p_sca_norm_2¡¿¡¾E_p_inc_norm_2¡¿%%%%%%%%%%%%%%%%%%%%%%%%
E_p_sca_norm_2 = sum(abs(E_s).^2,1);%[1,Ni]
E_p_inc_norm_2 = sum(abs(E_inc).^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except ...
     M Ni k0 imp  tx ty WasK d ...
     E_s chai0 ...
     Z Gs E_inc ...
     I E_p_sca_norm_2 E_p_inc_norm_2

 
%Next, we will carry out the iteration of the inversion based on CIE 
%within the framework of FBE/CSI.

%%%%%%%%%%%%%%%%%%%%%% Initialize the iteration parameters  %%%%%%%%%%%%%%% 
%Before the iteration starts, we need to initialize some parameters, 
%which will be updated in each iteration.

%1.<chai_n>-----the initialization of the contrast of the dielectric profile, which is set
%to be the background air;
chai_n = zeros(M^2,1);
%2.<a_n>-----the initialization of the induced current for each illuminations;
a_n=zeros(M^2,Ni);
%3.<y_n>-----the initialization of the Fourier transform of <a_n>;
y_n = zeros(M^2,Ni);
%4.<n>-----the number of total iterations;
n = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



%%%% Set some important parameters to be used in the 3-rounds of the nested scheme %%%%%%%%%%%%%%%%%%%%%%

iter = [300,300,2600]; %%% This is the parameter defined in Ref. [21] and [22], which controls the 
                      %%% number of the iterations in each round. The example here shows that,
                      %%% in the first round, the number of iterations is 
                      %%% 300; In the second round, the number of iterations
                      %%% is 300; In the third round, the number of
                      %%% iterations is 2600;
                      

                      
beta = [50,15,0.3];  %%% This is the parameter defined in Ref. [20], which is the core parameter
                     %%% being used in CIE. The example here shows that,
                     %%% in the first round, the value of \beta is 50;
                     %%% (Different from Ref. [20], the value of \beta used in the first round
                     %%% can be larger than 10; This is very crucial!)
                     %%% In the second round, the value of \beta is 15;
                     %%% In the third round, the value of \beta is 0.3;

MF = [5,7,10];       %%% This is the parameter defined in Ref. [15], which controls the 
                     %%% number of Fourier bases being used during the inversion. The
                     %%% example here shows that, in the first round, one
                     %%% uses MF = 5, then in the second round, using the
                     %%% results obtained in the first round as the initial guess, one
                     %%% uses M_F = 7. And so on.


alpha = [1,1,6];     %%% This is the parameter defined in the author's
                     %%% paper, which controls the weights of the residue of state
                     %%% equation. The example here shows that, in the first round, 
                     %%% the value of \alpha is 1; In the second round, the
                     %%% value of \alpha is 1; In the third round, the value of \alpha is 6;
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                               
%%%%%%%%%%%%%%   Main loop using FBE/CSI based on ICIE   %%%%%%%%%%%%%%%%%%
%%% Since our scheme is heavily using Contrast Source Inversion (CSI), 
%%% for a detailed theory of FBE-CIE, please refer to the literature [11] and [21]. 
for it_Num=1:size(MF,2)
    
A_0 = alpha(it_Num);

beta0 = beta(it_Num);
b_v = beta0*ones(M^2,1);



M_F =  MF(it_Num); 
MUSK = [ones(M_F,M_F),zeros(M_F,M-2*M_F),ones(M_F,M_F);...
    zeros(M-2*M_F,M);...
    ones(M_F,M_F),zeros(M_F,M-2*M_F),ones(M_F,M_F)];



R_n = (b_v.*chai_n)./(b_v.*chai_n+1);


for ii=1:iter(it_Num)

n = n+1;   

%%%%%%%%%%%%%  the residue calculations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% the residue of data equation: <p_n>
p_n=Gs*I*a_n-E_s;






%%%% the residue of state equation: <r_n>   
R_K = repmat(R_n,1,Ni);

PQK = repmat(b_v,1,Ni);

r_n= -R_K.*E_inc  ...
             +(PQK - R_K.*PQK)./(-1i*k0/imp).*a_n ...
             -R_K.*GGd(a_n,Z,M);


%%%  Calculate the Polak-Ribi`ere gradient of the induced current: <g_n>  %%%
fai_s = E_p_sca_norm_2;

fai_d = E_p_inc_norm_2./A_0; 
            

for rr=1:Ni
    
AAA = Gs'*p_n(:,rr);
  
BBB = ((b_v - R_n(:).*b_v).*(I./(-1i*k0/imp)))'*r_n(:,rr) ...
      -GGd(conj(R_n(:)).*(r_n(:,rr)),conj(Z),M);

temp1 = fft2(reshape(     AAA       ,M,M)).*MUSK;

temp2 = fft2(reshape(     BBB       ,M,M)).*MUSK;

temp3 = temp1./fai_s(rr) + temp2./fai_d(rr);
temp4 = reshape(temp3,M^2,1);
g_n(:,rr) = temp4;
end

%%%%%%%%%%  Update direction vector: <v_n>   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n==1
v_n=g_n;
end
if n>1
v_n = g_n+repmat(real(sum(conj(g_n-g_n_1).*g_n,1))./sum(abs(g_n_1).^2,1),M^2,1).*v_n; 
end


for hh=1:Ni
    temp5 = ifft2(reshape(v_n(:,hh),M,M));
    d_n(:,hh) = reshape(temp5,M^2,1);
end



%%%%  Update the step size parameter in CG optimization: <Num>\<Den>   %%%%

CCC =   Gs*d_n;
DDD =   ((PQK - R_K.*PQK)./(-1i*k0/imp)).*d_n ...
        -R_K .* GGd(d_n,Z,M);
      



Num = -sum(   conj(CCC)  .*  p_n     ,1)./fai_s  ...
      -sum(   conj(DDD)  .*  r_n     ,1)./fai_d;
  
Den = sum(abs(  CCC   ).^2,1)./fai_s ...
      +sum(abs(  DDD   ).^2,1)./fai_d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%%%%  Update the induced current: <a_n>    %%%%%%%%%%%%%%%%%%%%
y_n = y_n+repmat(Num./Den,M^2,1).*v_n; 

for yy=1:Ni
a_n(:,yy) = reshape(ifft2(reshape(y_n(:,yy),M,M)),M^2,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%% Update the estimated modified contrast function : <R_n> %%%%%%%%%%
J_tot_n = ( (PQK).*(I*a_n))./(-1i*k0/imp); 

E_tot_n=E_inc+GGd(I*a_n,Z,M)+J_tot_n;


chai_num = sum(conj(E_tot_n).*J_tot_n,2); 
chai_den = sum(conj(E_tot_n).*E_tot_n,2);

R_n = chai_num./chai_den; 

R_n(real(R_n)<0) = 0;

R_n = real(R_n);   % In this case, we restrict the scatterers to be lossless as a prior information                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




g_n_1=g_n;


%%%%%%%%%%%%%%%%       Calculate two kinds of errors   %%%%%%%%%%%%%%%%%%%%
p_tot = -R_n./(b_v.*R_n-b_v)+ones(M^2,1);q_tot = chai0+ones(M^2,1);
p_sca = p_tot.*WasK;q_sca = q_tot.*WasK;
p_sca(find(p_sca==0))=[];q_sca(find(q_sca==0))=[];

%%% Calculate the total error: <Er_tot>
Er_tot(n) = sqrt(norm((p_tot-q_tot)./q_tot)^2/M^2);

%%% Calculate the internal error: <Er_sct>
Er_sct(n) = sqrt(norm((p_sca-q_sca)./q_sca)^2/M^2);

%%% Display the two kinds of errors at each iteration
disp(['The number of iterations n=' num2str(n),',','the internal error <Er_sct> =',num2str(Er_sct(n)),',','the total error <Er_tot> =',num2str(Er_tot(n))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%%%%%%%%%%  Display the estimated profile at each round  %%%%%%%%%%%%%%%

%%% Tranform the modified contrast <R_n> into real contrast <chai_n>
chai_n = -R_n./(b_v.*R_n-b_v);
chai = reshape(chai_n,M,M);

%%% Display the real part of the estimated profile at each round
figure; 
imagesc(tx,ty,real(chai)+ones(M,M));
set(gca,'YDir','normal');
set(gcf,'color','w');
colormap(parula(256));
colorbar;

%%% Display the imaginary part of the estimated profile at each round
figure; 
imagesc(tx,ty,imag(chai));
set(gca,'YDir','normal');
set(gcf,'color','w');
colormap(parula(256));
colorbar;
end

%%% Display the two kinds of error after 600 iterations in one figure
figure;      
plot(Er_tot(600:end),'r-','LineWidth',3.5);
hold on
plot(Er_sct(600:end),'b-','LineWidth',3.5);
grid on
xlabel('Number of Iterations');
ylabel('Two kinds of error');
legend('the total error Er_{tot}', 'the internal error Er_{sct}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%