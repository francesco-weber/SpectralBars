% read the results of Colombini's spectral model
clear; clc;
close all;

% folder with results
fold = 'inputfiles/';

%% read input data
T = load([fold,'input.txt']);
beta = T(3);
lambda = T(4);
nx = T(6);
ny = T(7);
nxp = 2*nx;
nyp = ny+1;

textcase = ['b',num2str(beta),'_l',num2str(lambda),'(',int2str(nx),'x',int2str(ny),')'];

%% map
% A = load([fold,fname_map]);
% if(nxp*nyp~=length(A))
%     disp('error on nx*ny')
%     return
% end

% Open the file
maps_eta  = load("res_map_e.txt");
maps_u  = load("res_map_u.txt");
maps_v  = load("res_map_v.txt");
maps_d  = load("res_map_d.txt");

%%

numberofmaps=length(maps_eta)/(nxp*nyp);

evol_eta_fixed=zeros(numberofmaps-1);
evol_Kp=zeros(numberofmaps-1);
evol_K=zeros(numberofmaps-1);
evol_diss=zeros(numberofmaps-1);

m=length(maps_eta)/numberofmaps;


% Define parameters
numFrames = numberofmaps; % Number of frames in the movie
frameRate = 5; % Frames per second
duration = numFrames / frameRate; % Duration of the movie

% v_eta = VideoWriter('video_eta_b16_from_rand.mp4', 'MPEG-4');
% v_eta_section = VideoWriter('video_eta_b16_section_from_rand.mp4', 'MPEG-4');
% v_eta_plan = VideoWriter('video_eta_b16_plan_from_rand.mp4', 'MPEG-4');
% % v_u = VideoWriter('video_u_b18.mp4', 'MPEG-4');
% v_eta.FrameRate = frameRate;
% v_eta_section.FrameRate = frameRate;
% v_eta_plan.FrameRate = frameRate;
% % v_u.FrameRate = frameRate;
% open(v_eta);
% open(v_eta_section);
% open(v_eta_plan);
% % open(v_u);




for k=1:numberofmaps-1

    if k==1
        x=maps_eta(1:m,1);
        y=maps_eta(1:m,2);
        z=maps_eta(1:m,3);
        u=maps_u(1:m,3);
        v=maps_v(1:m,3);
        d=maps_d(1:m,3);

        x = reshape(x,nyp,nxp);
        y = reshape(y,nyp,nxp);
        z = reshape(z,nyp,nxp);
        u = reshape(u,nyp,nxp);
        v = reshape(v,nyp,nxp);
        d = reshape(d,nyp,nxp);

        z_sec= z(1,:);
        x_sec= x(1,:);
    
    
    
        % figure
        % % pcolor(x,y,z); shading flat
        % plot(x_sec,z_sec);
        % ylim([-0.3,0.3]);
        % pbaspect([10, 2, 1]);
        % % surf(x,y,z);
        % frame = getframe(gcf);
        % writeVideo(v_eta_section, frame);
        % close


        U=mean(u(:));
        V=mean(v(:));

        Ekin=u.^2+v.^2;
        Ekin_pert=(u-U).^2+(v-V).^2;
        Ediss=0.064*1000*(sqrt(u.^2+v.^2)).^3;
% 
        Ekin_dxdy=1/2*Ekin.*(1/nxp)*(1/nyp)*1000.*d;
        Ekin_dxdy_pert=1/2*Ekin_pert.*(1/nxp)*(1/nyp)*1000.*d;
% 
        Ekin_tot=sum(Ekin_dxdy(:));
        Ekin_tot_pert=sum(Ekin_dxdy_pert(:));

        evol_eta_fixed(k)=z(1,1);
        evol_K(k)=Ekin_tot;
        evol_Kp(k)=Ekin_tot_pert;
        evol_diss(k)=sum(Ediss(:));
        
        

        % 
        % figure
        % surfc(x,y,z);
        % %colormap gray;
        % shading interp;
        % pbaspect([10, 2, 1]);
        % zlim([-1,1]);
        % title('Evolution of bed elevation','Interpreter','latex');
        % xlabel('x [-]','Interpreter','latex');
        % ylabel('y/$\beta$ [-]','Interpreter','latex');
        % zlabel('$\eta$ [-]','Interpreter','latex');
        % frame = getframe(gcf);
        % writeVideo(v_eta, frame);    
        % close

        % figure
        % % pcolor(x,y,z); shading flat
        % contourf(x,y,z);
        % pbaspect([10, 1, 1]);
        % zlim([-1,0.5]);
        % colorbar
        % xlabel('x')
        % ylabel('y')
        % title('map of bottom elevation')
        % frame = getframe(gcf);
        % writeVideo(v_eta_plan, frame);    
        % close


        % figure
        % surfc(x,y,u);
        % %colormap gray;
        % pbaspect([10, 2, 1]);
        % zlim([-1,1]);
        % frame = getframe(gcf);
        % writeVideo(v_u, frame);
        % close
    end

    x=maps_eta(((k-1)*m+1):k*m,1);
    y=maps_eta(((k-1)*m+1):k*m,2);
    z=maps_eta(((k-1)*m+1):k*m,3);
    u=maps_u(((k-1)*m+1):k*m,3);
    v=maps_v(((k-1)*m+1):k*m,3);
    d=maps_d(((k-1)*m+1):k*m,3);

    x = reshape(x,nyp,nxp);
    y = reshape(y,nyp,nxp);
    z = reshape(z,nyp,nxp);
    u = reshape(u,nyp,nxp);
    v = reshape(v,nyp,nxp);
    d = reshape(d,nyp,nxp);

    z_sec= z(1,:);
    x_sec= x(1,:);



    % figure
    % % pcolor(x,y,z); shading flat
    % plot(x_sec,z_sec);
    % ylim([-0.3,0.3]);
    % pbaspect([10, 2, 1]);
    % % surf(x,y,z);
    % frame = getframe(gcf);
    % writeVideo(v_eta_section, frame);
    % close

    U=mean(u(:));
    V=mean(v(:));

    Ekin=u.^2+v.^2;
    Ekin_pert=(u-U).^2+(v-V).^2;
    Ediss=0.064*1000*(sqrt(u.^2+v.^2)).^3;
% 
    Ekin_dxdy=1/2*Ekin.*(1/nxp)*(1/nyp)*1000.*d;
    Ekin_dxdy_pert=1/2*Ekin_pert.*(1/nxp)*(1/nyp)*1000.*d;
% 
    Ekin_tot=sum(Ekin_dxdy(:));
    Ekin_tot_pert=sum(Ekin_dxdy_pert(:));

    evol_eta_fixed(k)=z(1,1);
    evol_K(k)=Ekin_tot;
    evol_Kp(k)=Ekin_tot_pert;
    evol_diss(k)=sum(Ediss(:));

    % figure
    % surfc(x,y,z);
    % %colormap gray;
    % shading interp;
    % pbaspect([10, 2, 1]);
    % zlim([-1,1]);
    % title('Evolution of bed elevation','Interpreter','latex');
    % xlabel('x [-]','Interpreter','latex');
    % ylabel('y/$\beta$ [-]','Interpreter','latex');
    % zlabel('$\eta$ [-]','Interpreter','latex');
    % frame = getframe(gcf);
    % writeVideo(v_eta, frame);
    % close
    % 
    % figure
    % contourf(x,y,z);
    % pbaspect([10, 1, 1]);
    % zlim([-1,0.5]);
    % colorbar
    % xlabel('x')
    % ylabel('y')
    % title('map of bottom elevation')
    % frame = getframe(gcf);
    % writeVideo(v_eta_plan, frame);    
    % close

    % figure
    % surfc(x,y,u);
    % %colormap gray;
    % pbaspect([10, 2, 1]);
    % zlim([-1,1]);
    % frame = getframe(gcf);
    % writeVideo(v_u, frame);
    % close


% shading interp;
% % axis equal;



end



figure
plot(evol_K(:,1));
title('evolution of Kinetic Energy','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$K$','Interpreter','latex');
savefig(['K_t_', textcase,'.fig'])


figure
plot(evol_Kp(:,1));
title('evolution of Perurbation Kinetic Energy','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$K_p$','Interpreter','latex');
savefig(['Kp_t_', textcase, '.fig'])


figure
plot(evol_diss(:,1));
title('evolution of dissipation','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$\epsilon$','Interpreter','latex');
savefig(['eps_t_', textcase,'.fig'])


figure
plot(evol_eta_fixed(:,1));
title('Evolution of pointwise bed elevation $\eta$','Interpreter','latex');
xlabel('Time','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');
savefig(['eta_t_', textcase, '.fig'])

