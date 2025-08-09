% Linear stability analysis of free bars in dimensionless form
% Plots Growth rate as a function of (beta, theta0, ds0) 


clear; clc;
close all;

imma = 1i;
Delta = 1.65;        % relative sediment density
por = 0.4;           % porosity
g = 0.981e1;         % gravitational acceleration 
r = 0.4e0;           % Ikeda slope parameter

ds=0.03;             % base uniform flow 
theta__o=0.2;        % base uniform flow Shields number


% uniform transport capacity according to MPM (1947) Phi__o = M * (mu * So + theta__o - vartheta__o) ^ m;

vartheta__o = 0.047; % critical Shields number on horizontal bed
M=8;                 % MPM coefficient           
m = 3/2;             % MPM exponent
mu = 0.1;            % longitudinal slope correction coefficient


Ch = 0.6e1 + 0.25e1 * log(0.4000000000e0 / ds);  % Chezy dimensionless
Cd = 0.1e1 / Ch ^ 2;                             % Drag coefficient

%% Computation of uniform flow variables
So = theta__o * Delta * ds;   % reference dimensionless slopes
F = Ch * sqrt(So);
iF = 0.1e1 / F ^ 2;
Qo = -sqrt(Delta) * ds ^ (0.3e1 / 0.2e1) * So ^ (-0.1e1 / 0.2e1) / (-0.1e1 + por) / Ch;

Phi__o = M * (mu * So + theta__o - vartheta__o) ^ m;
Phi__m = m / (mu * So + theta__o - vartheta__o);
r1 = r * theta__o ^ (-0.1e1 / 0.2e1);



% Lateral structure of the perturbation
n__pi = 1/2;
P = pi * n__pi;

% coefficients of the linear matrix L
c__U = (0.2e1 * So + 0.2e1) * Phi__m * theta__o * Phi__o * Qo;
c__V = -P * Phi__o * Qo;
c__0 = Phi__o * Qo * r1 * P ^ 2;
c__k = (mu * Phi__m ) * Phi__o * Qo;



%% Stability analysis

cBK(1,2) = -P ^ 4 * iF ^ 2 * c__0;
cBK(1,3) = -iF * P ^ 2 * ((P ^ 2 * c__k + 0.2e1 * c__0) * iF - 0.2e1 * c__0);
cBK(1,4) = -0.2e1 * c__k * (-0.1e1 + iF) * iF * P ^ 2 - (-0.1e1 + iF) ^ 2 * c__0;
cBK(1,5) = -c__k * (-0.1e1 + iF) ^ 2;
cBK(2,1) = -0.4e1 * P ^ 4 * Cd ^ 2 * iF ^ 2 * c__0;
cBK(2,2) = -0.4e1 * Cd * (((P ^ 2 * c__k + c__0) * Cd - P * (P * c__U + c__V) / 0.4e1) * iF - 0.5e1 / 0.2e1 * Cd * c__0) * iF * P ^ 2;
cBK(2,3) = -0.2e1 * (0.2e1 * c__k * (iF - 0.5e1 / 0.2e1) * Cd ^ 2 - (iF - 0.1e1 / 0.2e1) * c__U * Cd) * iF * P ^ 2 + (0.2e1 * iF + 0.1e1) * c__V * Cd * iF * P - 0.9e1 * Cd ^ 2 * c__0 - Cd ^ 2 * (-0.1e1 + iF) ^ 2 * c__0;
cBK(2,4) = -Cd * (c__k * (iF ^ 2 - 0.2e1 * iF + 0.10e2) * Cd + 0.3e1 * iF * c__U);
cBK(3,2) = 0.3e1 * Cd ^ 3 * (-0.3e1 * Cd * c__0 + P * iF * (P * c__U + 0.3e1 * c__V));
cBK(3,3) = -0.3e1 * Cd ^ 3 * (0.3e1 * Cd * c__k + iF * c__U);
cBKr(2) = cBK(1,2) * cBK(2,1) ^ 2 * cBK(3,2) ^ 2;
cBKr(3) = cBK(2,1) ^ 2 * (0.3e1 * cBK(1,2) * cBK(3,2) * cBK(3,3) + 0.2e1 * cBK(1,3) * cBK(3,2) ^ 2);
cBKr(4) = cBK(2,1) * (0.2e1 * cBK(1,2) * cBK(2,1) * cBK(3,3) ^ 2 + cBK(1,2) * cBK(2,2) * cBK(3,2) * cBK(3,3) - 0.2e1 * cBK(1,2) * cBK(2,3) * cBK(3,2) ^ 2 + 0.6e1 * cBK(1,3) * cBK(2,1) * cBK(3,2) * cBK(3,3) + cBK(1,3) * cBK(2,2) * cBK(3,2) ^ 2 + 0.3e1 * cBK(1,4) * cBK(2,1) * cBK(3,2) ^ 2);
cBKr(5) = cBK(2,1) * ((-0.4e1 * cBK(1,2) * cBK(2,4) - 0.2e1 * cBK(1,3) * cBK(2,3) + 0.2e1 * cBK(1,4) * cBK(2,2) + 0.4e1 * cBK(1,5) * cBK(2,1)) * cBK(3,2) ^ 2 - 0.4e1 * cBK(3,3) * (cBK(1,2) * cBK(2,3) - 0.5e1 / 0.4e1 * cBK(1,3) * cBK(2,2) - 0.9e1 / 0.4e1 * cBK(1,4) * cBK(2,1)) * cBK(3,2) + cBK(3,3) ^ 2 * (cBK(1,2) * cBK(2,2) + 0.4e1 * cBK(1,3) * cBK(2,1)));
cBKr(6) = cBK(1,3) ^ 2 * cBK(3,2) ^ 3 + (-0.2e1 * cBK(1,2) * cBK(3,3) * cBK(1,3) + (-0.5e1 * cBK(1,3) * cBK(2,4) - 0.2e1 * cBK(1,4) * cBK(2,3) + 0.3e1 * cBK(1,5) * cBK(2,2)) * cBK(2,1) + cBK(1,2) * cBK(2,3) ^ 2 - cBK(1,3) * cBK(2,2) * cBK(2,3)) * cBK(3,2) ^ 2 + cBK(3,3) * (cBK(1,2) ^ 2 * cBK(3,3) + 0.12e2 * cBK(1,5) * cBK(2,1) ^ 2 + (-0.9e1 * cBK(1,2) * cBK(2,4) - 0.2e1 * cBK(1,3) * cBK(2,3) + 0.9e1 * cBK(1,4) * cBK(2,2)) * cBK(2,1) - cBK(1,2) * cBK(2,2) * cBK(2,3) + cBK(1,3) * cBK(2,2) ^ 2) * cBK(3,2) - 0.2e1 * cBK(3,3) ^ 2 * cBK(2,1) * (cBK(1,2) * cBK(2,3) - 0.2e1 * cBK(1,3) * cBK(2,2) - 0.3e1 * cBK(1,4) * cBK(2,1));
cBKr(7) = 0.4e1 * cBK(1,3) * cBK(1,4) * cBK(3,2) ^ 3 + ((-0.4e1 * cBK(1,2) * cBK(1,4) + cBK(1,3) ^ 2) * cBK(3,3) + (-0.2e1 * cBK(1,3) * cBK(2,4) - 0.2e1 * cBK(1,4) * cBK(2,3)) * cBK(2,2) + 0.4e1 * cBK(1,2) * cBK(2,3) * cBK(2,4) + (-0.6e1 * cBK(1,4) * cBK(2,4) - 0.2e1 * cBK(1,5) * cBK(2,3)) * cBK(2,1)) * cBK(3,2) ^ 2 - 0.2e1 * cBK(3,3) * (cBK(1,2) * cBK(3,3) * cBK(1,3) - cBK(1,4) * cBK(2,2) ^ 2 + (cBK(1,2) * cBK(2,4) + cBK(1,3) * cBK(2,3) / 0.2e1 - 0.13e2 / 0.2e1 * cBK(1,5) * cBK(2,1)) * cBK(2,2) - cBK(1,2) * cBK(2,3) ^ 2 / 0.2e1 + 0.9e1 / 0.2e1 * cBK(1,3) * cBK(2,1) * cBK(2,4)) * cBK(3,2) + cBK(3,3) ^ 2 * (cBK(1,2) ^ 2 * cBK(3,3) + cBK(1,3) * cBK(2,2) ^ 2 + (-cBK(1,2) * cBK(2,3) + 0.7e1 * cBK(1,4) * cBK(2,1)) * cBK(2,2) - 0.5e1 * cBK(1,2) * cBK(2,1) * cBK(2,4) + 0.8e1 * cBK(1,5) * cBK(2,1) ^ 2);
cBKr(8) = (0.6e1 * cBK(1,3) * cBK(1,5) + 0.4e1 * cBK(1,4) ^ 2) * cBK(3,2) ^ 3 + ((-0.6e1 * cBK(1,2) * cBK(1,5) + 0.6e1 * cBK(1,3) * cBK(1,4)) * cBK(3,3) + (-0.4e1 * cBK(2,2) * cBK(2,4) - cBK(2,3) ^ 2) * cBK(1,4) + 0.4e1 * cBK(1,2) * cBK(2,4) ^ 2 + (cBK(1,3) * cBK(2,3) - 0.7e1 * cBK(1,5) * cBK(2,1)) * cBK(2,4) - 0.3e1 * cBK(1,5) * cBK(2,2) * cBK(2,3)) * cBK(3,2) ^ 2 - 0.6e1 * cBK(3,3) * (cBK(1,2) * cBK(3,3) * cBK(1,4) + (0.3e1 / 0.2e1 * cBK(2,1) * cBK(2,4) + cBK(2,2) * cBK(2,3) / 0.6e1) * cBK(1,4) + (-0.5e1 / 0.6e1 * cBK(1,2) * cBK(2,3) + 0.2e1 / 0.3e1 * cBK(1,3) * cBK(2,2)) * cBK(2,4) - (cBK(2,1) * cBK(2,3) + 0.3e1 / 0.2e1 * cBK(2,2) ^ 2) * cBK(1,5) / 0.3e1) * cBK(3,2) - 0.2e1 * cBK(3,3) ^ 2 * ((-cBK(2,1) * cBK(2,3) - cBK(2,2) ^ 2) * cBK(1,4) + (cBK(1,2) * cBK(2,2) + 0.2e1 * cBK(1,3) * cBK(2,1)) * cBK(2,4) - 0.5e1 * cBK(2,2) * cBK(1,5) * cBK(2,1));
cBKr(9) = -0.2e1 * cBK(1,2) * cBK(1,4) * cBK(3,3) ^ 3 + (0.2e1 * (-0.5e1 * cBK(1,2) * cBK(1,5) + cBK(1,3) * cBK(1,4)) * cBK(3,2) + (cBK(1,2) * cBK(2,3) - 0.2e1 * cBK(1,3) * cBK(2,2) - 0.3e1 * cBK(1,4) * cBK(2,1)) * cBK(2,4) + cBK(1,4) * cBK(2,2) * cBK(2,3) + cBK(1,5) * (0.4e1 * cBK(2,1) * cBK(2,3) + 0.3e1 * cBK(2,2) ^ 2)) * cBK(3,3) ^ 2 + 0.6e1 * ((0.5e1 * cBK(1,3) * cBK(1,5) + 0.4e1 * cBK(1,4) ^ 2) * cBK(3,2) / 0.3e1 + cBK(1,2) * cBK(2,4) ^ 2 + (cBK(1,3) * cBK(2,3) / 0.6e1 - cBK(1,4) * cBK(2,2) - 0.3e1 / 0.2e1 * cBK(1,5) * cBK(2,1)) * cBK(2,4) - cBK(2,3) * (cBK(1,4) * cBK(2,3) + cBK(1,5) * cBK(2,2)) / 0.6e1) * cBK(3,2) * cBK(3,3) + 0.2e1 * (0.6e1 * cBK(1,4) * cBK(1,5) * cBK(3,2) + cBK(1,3) * cBK(2,4) ^ 2 + (-cBK(1,4) * cBK(2,3) - 0.3e1 * cBK(1,5) * cBK(2,2)) * cBK(2,4) - cBK(1,5) * cBK(2,3) ^ 2) * cBK(3,2) ^ 2;
cBKr(10) = -0.4e1 * cBK(1,2) * cBK(1,5) * cBK(3,3) ^ 3 + (0.2e1 * (0.2e1 * cBK(1,3) * cBK(3,2) - cBK(2,1) * cBK(2,4) + cBK(2,2) * cBK(2,3)) * cBK(1,5) + 0.2e1 * cBK(1,2) * cBK(2,4) ^ 2 + 0.5e1 * cBK(1,4) ^ 2 * cBK(3,2) - 0.2e1 * cBK(2,4) * cBK(1,4) * cBK(2,2)) * cBK(3,3) ^ 2 + 0.3e1 * (0.2e1 / 0.3e1 * (0.13e2 * cBK(1,4) * cBK(3,2) - 0.4e1 * cBK(2,2) * cBK(2,4) - cBK(2,3) ^ 2) * cBK(1,5) + cBK(1,3) * cBK(2,4) ^ 2 - cBK(1,4) * cBK(2,3) * cBK(2,4)) * cBK(3,2) * cBK(3,3) + 0.9e1 * cBK(1,5) ^ 2 * cBK(3,2) ^ 3 - 0.5e1 * cBK(1,5) * cBK(2,3) * cBK(2,4) * cBK(3,2) ^ 2;
cBKr(11) = cBK(1,4) ^ 2 * cBK(3,3) ^ 3 + ((0.18e2 * cBK(1,4) * cBK(3,2) - 0.2e1 * cBK(2,2) * cBK(2,4)) * cBK(1,5) + cBK(1,3) * cBK(2,4) ^ 2 - cBK(1,4) * cBK(2,3) * cBK(2,4)) * cBK(3,3) ^ 2 + (0.21e2 * cBK(1,5) ^ 2 * cBK(3,2) ^ 2 - 0.7e1 * cBK(1,5) * cBK(2,3) * cBK(2,4) * cBK(3,2)) * cBK(3,3) - 0.2e1 * cBK(1,5) * cBK(2,4) ^ 2 * cBK(3,2) ^ 2;
cBKr(12) = cBK(1,5) * cBK(3,3) * (0.4e1 * cBK(1,4) * cBK(3,3) ^ 2 + 0.16e2 * cBK(1,5) * cBK(3,2) * cBK(3,3) - 0.2e1 * cBK(2,3) * cBK(2,4) * cBK(3,3) - 0.3e1 * cBK(2,4) ^ 2 * cBK(3,2));
cBKr(13) = 0.4e1 * cBK(1,5) ^ 2 * cBK(3,3) ^ 3 - cBK(1,5) * cBK(2,4) ^ 2 * cBK(3,3) ^ 2;

cOn(1,3) = -4 * P ^ 4 * Cd ^ 2 * iF ^ 2 * c__0;
cOn(3,1) = -P ^ 4 * iF ^ 2 * c__0;
cOn(3,3) = -(0.4e1 * Cd * ((((P ^ 2 * c__k + c__0) * Cd) - (P * (P * c__U + c__V)) / 0.4e1) * iF - 0.5e1 / 0.2e1 * Cd * c__0) * iF * (P ^ 2));
cOn(3,5) = 3 * Cd ^ 3 * (-3 * Cd * c__0 + P * iF * (P * c__U + 3 * c__V));
cOn(5,1) = -iF * P ^ 2 * ((P ^ 2 * c__k + 2 * c__0) * iF - 2 * c__0);
cOn(5,3) = -Cd * (4 * Cd * P ^ 2 * c__k * iF ^ 2 - 10 * Cd * P ^ 2 * c__k * iF - 2 * P ^ 2 * c__U * iF ^ 2 + Cd * c__0 * iF ^ 2 + P ^ 2 * c__U * iF - 2 * P * c__V * iF ^ 2 - 2 * Cd * c__0 * iF - P * c__V * iF + 10 * Cd * c__0);
cOn(5,5) = -9 * Cd ^ 4 * c__k - 3 * Cd ^ 3 * c__U * iF;
cOn(7,1) = -(2 * P ^ 2 * c__k * iF + c__0 * iF - c__0) * (-1 + iF);
cOn(7,3) = -Cd * (c__k * (iF ^ 2 - 2 * iF + 10) * Cd + 3 * iF * c__U);
cOn(9,1) = -c__k * (-1 + iF) ^ 2;
cOd(1,4) = 4 * P ^ 4 * Cd ^ 2 * iF ^ 2;
cOd(3,2) = P ^ 4 * iF ^ 2;
cOd(3,4) = 2 * Cd ^ 2 * P ^ 2 * iF * (2 * iF - 5);
cOd(3,6) = 9 * Cd ^ 4;
cOd(5,2) = 2 * iF * P ^ 2 * (-1 + iF);
cOd(5,4) = Cd ^ 2 * (iF ^ 2 - 2 * iF + 10);
cOd(7,2) = (-1 + iF) ^ 2;
con(2,3) = 2 * iF ^ 2 * Cd ^ 2 * P ^ 3 * (P * c__U + 3 * c__V);
con(4,1) = iF ^ 2 * P ^ 3 * c__V;
con(4,3) = -Cd ^ 2 * iF * ((P * c__U - 3 * c__V) * iF + 4 * P * c__U + 9 * c__V) * P;
con(6,1) = -iF * P * ((P * c__U - c__V) * iF + c__V);
con(6,3) = -iF * Cd ^ 2 * c__U * (-1 + iF);
con(8,1) = -iF * c__U * (-1 + iF);
cod(1,3) = 4 * P ^ 4 * Cd ^ 2 * iF ^ 2;
cod(3,1) = P ^ 4 * iF ^ 2;
cod(3,3) = 2 * Cd ^ 2 * P ^ 2 * iF * (2 * iF - 5);
cod(3,5) = 9 * Cd ^ 4;
cod(5,1) = 2 * iF * P ^ 2 * (-1 + iF);
cod(5,3) = Cd ^ 2 * (iF ^ 2 - 2 * iF + 10);
cod(7,1) = (-1 + iF) ^ 2;


OmegaLin = @(k,b) ...
    ( b.^4 .* k.^4 .* cOn(5,5) + b.^2 .* k.^6 .* cOn(7,3) + k.^8 .* cOn(9,1) + ...
      b.^4 .* k.^2 .* cOn(3,5) + b.^2 .* k.^4 .* cOn(5,3) + k.^6 .* cOn(7,1) + ...
      b.^2 .* k.^2 .* cOn(3,3) + k.^4 .* cOn(5,1) + b.^2 .* cOn(1,3) + k.^2 .* cOn(3,1) ) ...
  ./ ( b.^5 .* k.^2 .* cOd(3,6) + b.^3 .* k.^4 .* cOd(5,4) + b    .* k.^6 .* cOd(7,2) + ...
       b.^3 .* k.^2 .* cOd(3,4) + b    .* k.^4 .* cOd(5,2) + b.^3 .* cOd(1,4) + ...
       b    .* k.^2 .* cOd(3,2) );

% find the roots of the polynomial
n_poly = 12; % max degree of the polynomial [note that cBKr(0) = 0]
coefs_poly = zeros(n_poly+1,1);
for i=1:n_poly
    coefs_poly(n_poly-i+1) = cBKr(i);
end
K_roots = roots(coefs_poly);
K_roots(imag(K_roots)~=0) = []; % remove imaginary roots
K_roots(K_roots<=0) = [];       % remove negative or zero roots
K = min(K_roots);               % select the minimum root

% compute Beta, beta and k
BETA = -K ^ 2 * (2 * K ^ 3 * cBK(1,5) * cBK(3,3) + K ^ 2 * cBK(1,4) * cBK(3,3) + 3 * K ^ 2 * cBK(1,5) * cBK(3,2) + 2 * K * cBK(1,4) * cBK(3,2) - cBK(1,2) * cBK(3,3) + cBK(1,3) * cBK(3,2)) / (K ^ 4 * cBK(2,4) * cBK(3,3) + 2 * K ^ 3 * cBK(2,4) * cBK(3,2) - K ^ 2 * cBK(2,2) * cBK(3,3) + K ^ 2 * cBK(2,3) * cBK(3,2) - 2 * K * cBK(2,1) * cBK(3,3) - cBK(2,1) * cBK(3,2));
beta = sqrt(BETA);
k = sqrt(K);

% critical beta
beta0 = beta;

%% ----------------------------- Plotting -----------------------------------
range_k=linspace(0.001,1,1000);
range_beta=linspace(1,20,1000);

Om=zeros(100,100);

for i=1:length(range_k)
    for j=1:length(range_beta)
        Om(i,j)=OmegaLin(range_k(1,i),range_beta(1,j));
    end
end


figure('Color','w');
% heatmap-style surface
surf(range_k, range_beta, Om.', 'EdgeColor','none'); % transpose to match axes
view(0,90);
colormap parula;
colorbar;
xlabel('$k$',Interpreter='latex',FontSize=18,Rotation=0); ylabel('$\beta$',Interpreter='latex',FontSize=18,Rotation=0); title('Growth rate $\Omega(k,\beta)$',Interpreter='latex',FontSize=18);
zlim([0 1e-3])
clim([-1e-3 1e-3])
hold on;
% neutral stability curve (Omega=0)
contour(range_k, range_beta, Om.', [0 0], 'k-', 'LineWidth', 1.2);
scatter(k,beta,Marker="o")

grid on;


