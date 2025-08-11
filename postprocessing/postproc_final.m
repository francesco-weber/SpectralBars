% read the results of Colombini's spectral model
clear; clc;
close all;

% folder with results
fold = 'inputfiles/';

%% read input data
T = load([fold,'input.txt']);
beta = T(3);
theta=T(1);
lambda = T(4);
nx = T(6);
ny = T(7);
nxp = 2*nx;
nyp = ny+1;

textcase = ['th',num2str(theta),'_b',num2str(beta),'_l',num2str(lambda),'(',int2str(nx),'x',int2str(ny),')'];
%textcase = ['th0.18_beta14'];

%% map
% A = load([fold,fname_map]);
% if(nxp*nyp~=length(A))
%     disp('error on nx*ny')
%     return
% end

% Open the file
maps_eta  = load("res_map_e.txt");
%%

numberofmaps=length(maps_eta)/(nxp*nyp);

evol_eta_fixed=zeros(numberofmaps-1);
evol_wavenumber_sel=zeros(numberofmaps-1);


m=length(maps_eta)/numberofmaps;






for k=1:numberofmaps-3

    if k==1
        x=maps_eta(1:m,1);
        y=maps_eta(1:m,2);
        z=maps_eta(1:m,3);

        x = reshape(x,nyp,nxp);
        y = reshape(y,nyp,nxp);
        z = reshape(z,nyp,nxp);

        z_sec= z(1,:);
        x_sec= x(1,:);
    



        evol_eta_fixed(k)=z(3,1);
        
        

 

    end

    x=maps_eta(((k-1)*m+1):k*m,1);
    y=maps_eta(((k-1)*m+1):k*m,2);
    z=maps_eta(((k-1)*m+1):k*m,3);

    x = reshape(x,nyp,nxp);
    y = reshape(y,nyp,nxp);
    z = reshape(z,nyp,nxp);

    z_sec= z(1,:);
    x_sec= x(1,:);

    %% data
    signal = z_sec; % resampled data
    N = length(signal); %total number of samples
    X = 2*pi/lambda; % total duration of the signal in hours
    f0 = 1/X; % fundamental frequency
    % frequency
    f = (1:N/2)*f0;
    
    
    %% Using the Fourier transform
    FT= fft(signal,N); % determining the a_n coefficients
    c_n = FT./N;
    A0 = real(c_n(1));
    a_n = c_n(2:floor(N/2)+1);
    
    % Computing amplitude and phase of each harmonic
    % A = 2*abs(a_n); %amplitude of the first n_harmonics harmonics
    % phase = angle(a_n);
    % figure;
    % bar([A0 A']);
    
    % Computing the Power Spectral Density
    PSD = abs(a_n).^2/f0;
    power_spectrum = abs(a_n).^2/f0;
    
    %% Plotting the PSD
    %peaks
    mean_frequency = sum(f .* power_spectrum) / sum(power_spectrum);
    evol_wavenumber_sel(k)=mean_frequency;

    



   
    evol_eta_fixed(k)=z(3,1);
   
end

figure
surfc(x,y,z);
%colormap gray;
shading interp;
pbaspect([10, 2, 1]);
%zlim([-1,1]);
title('Evolution of bed elevation','Interpreter','latex');
xlabel('x [-]','Interpreter','latex');
ylabel('y [-]','Interpreter','latex');
zlabel('$\eta$ [-]','Interpreter','latex');
frame = getframe(gcf);
savefig(sprintf('surf_final_%s.fig', textcase));



figure
contourf(x,y,z);
pbaspect([10, 1, 1]);
%zlim([-1,0.5]);
colorbar
xlabel('x')
ylabel('y')
title('map of bottom elevation')
frame = getframe(gcf);
savefig(sprintf('contour_final_%s.fig', textcase));





figure
plot(evol_eta_fixed(:,1));
title('Evolution of pointwise bed elevation $\eta$','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');
frame = getframe(gcf);
savefig(sprintf('eta_t_%s.fig', textcase));



figure
plot(evol_wavenumber_sel(:,1));
title('Evolution of selected $\lambda$','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$\lambda$','Interpreter','latex');
savefig(sprintf('lambda_t_%s.fig', textcase));

