%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%        Heat Diffusion in 2D Plate Using Finite Difference Method          %%%%%%%%%%%%%
%%%%%%%%%%%%%                   Numerical Anylysis cource project                       %%%%%%%%%%%%%
%%%%%%%%%%%%%                          Dr. Mohamed Tawfik                               %%%%%%%%%%%%%
%%%%%%%%%%%%%                           Amr Mousa Mohamed                               %%%%%%%%%%%%%
%%%%%%%%%%%%% Space and Communication Engineering - Zewail City of Science & Technology %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is designed to solve the heat equation in a 2D plate.
% Using fixed boundry conditions "Dirichlet Conditions" and initial
% temperature in all nodes, It can solve until reach steady state with
% tolerence value selected below.
% After solution, graphical simulation appears to show you how the heat
% diffuses throughout the plate within time interval selected below

clear; close all; clc;


 %% 1-Inputs section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Please select your material, enter your parameters and your initial conditions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%                              -- Aluminum --            % ** Uncomment if needed **
% name=('Aluminium');                                        % Material name
% conductivity = 204.3;                                      % thermal conductivity (j/m.C.sec)
% spacific_heat = 910;                                       % specific heat (j/kg.C)
% denisty = 2700.0;                                          % density (kg/m^3)

% %%%                               -- Copper --             % ** Uncomment if needed **
% name=('Copper');                                           % Material name
% conductivity = 401;                                        % thermal conductivity (W/m.K)
% spacific_heat = 390;                                       % specific heat (J/kg K)
% denisty = 8940;                                            % density (kg/m^3)

%%%                               -- Silver --
name=('Silver');                                             % Material name
conductivity = 629;                                          % thermal conductivity (W/m.K)
spacific_heat = 233;                                         % specific heat (J/kg K)
denisty = 10490;                                             % density (kg/m^3)

% %%%                          -- Custom Material --         % ** Uncomment and enter your values if needed **
% name=('Custom Material');                                  % Material name
% conductivity =    ;                                        % thermal conductivity (W/m.K)
% spacific_heat =    ;                                       % specific heat (J/kg K)
% denisty =      ;                                           % density (kg/m^3)

%%
Lx= 1;                             % plate width (m)
Ly= 1;                             % plate length (m)
Nx=40;                             % nodes in x direction
Ny=40;                             % nodes in y direction

T_initial= 0 ;                   % Initial temperature in all nodes ( the whole plate )
T_east   = 150 ;                   % temperature on the upper side ( at y=0  "Dirichlet Conditions" )
T_west   = 300 ;                   % temperature on the lower side ( at y=Ly "Dirichlet Conditions" )
T_north  = 50 ;                   % temperature on the left  side ( at x=0  "Dirichlet Conditions" )
T_south  = 100 ;                   % temperature on the right side ( at x=Lx "Dirichlet Conditions" ) 

t_end=100 ;                        % final time for visual simulation (sec)
dt=0.6 ;                           % time step (1 sec)

tolerence = 0.5;                   % tolerence for numerical simulation (0.5 deg Celesius)
tolerence_ss=0.001;                % tolerence for steady state section (0.1 deg Celesius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% From here, You don't need to modify %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2- Constants Section
 
k=1;                                                                       % iteration counter
err_SS_max(k)=1;                                                                  % initial error
err_SS_min(k)=1;                                                                  % initial error
dx=Lx/Nx;                                                                  % delta x
dy=Ly/Ny;                                                                  % delta y
n_time=round(t_end/dt);                                                           % number of iterrations for time
alpha = conductivity/(spacific_heat*denisty);                              % alpha (1/sec)
T_max=max([T_east T_west T_north T_south T_initial]);                      % Max T to set axes limits in plotting
T_min=min([T_east T_west T_north T_south T_initial]);                      % Min T to set axes limits in plotting
Solution_type=questdlg('Which method you want to solve the time derivative with ?','Question','Euler','2nd order Runge-Kutte','Euler');                     % solve with 2nd order Runge Kutte in time or 2 to solve with Euler

if dt<= 1/(2*alpha*((1/dx^2)+(1/dy^2)))                             % test the stability condition 
else 
    fprintf('Error, the stability condition is not met\nPlease return to "Inputs Section" and choose a "dt" smaller than %f \n',1/(2*alpha*((1/dx^2)+(1/dy^2))))
    return
end
message=msgbox('Your computer is now solving the problem, Please wait..... ');    % Busy message 
% ----------------- Initial Conditions for finite difference section ---------------
T=zeros(Nx+2,Ny+2,75000);                       % set max iterations 75,000 due to memory limitations (T variable takes maximum 1GB in memory)
T(:,1,:)=T_south;
T(:,Ny+1,:)=T_north;
T(:,Ny+2,:)=T_north;                            % Redundant, it has no effect in calculations but is required in plotting section
T(Nx+1,:,:)=T_east;
T(Nx+2,:,:)=T_east;                             % Redundant, it has no effect in calculations but is required in plotting section
T(1,:,:)=T_west;
T(:,:,1)=T_initial;
% ------------------- Initial Conditions for steady state section -------------------
Tss=zeros(Nx+2,Ny+2);        Tss2=zeros(Nx+2,Ny+2);
Tss(:,1)=T_south;            Tss2(:,1)=T_south;
Tss(:,Ny+1)=T_north;         Tss2(:,Ny+1)=T_north;
Tss(:,Ny+2)=T_north;         Tss2(:,Ny+2)=T_north;             % Redundant, it has no effect in calculations but is required in plotting section
Tss(Nx+1,:)=T_east;          Tss2(Nx+1,:)=T_east;
Tss(Nx+2,:)=T_east;          Tss2(Nx+2,:)=T_east;              % Redundant, it has no effect in calculations but is required in plotting section
Tss(1,:)=T_west;             Tss2(1,:)=T_west;


%% 3- Steady-State section


   while err_SS_max(k)>=tolerence_ss || err_SS_min(k)>=tolerence_ss 
    
    for i=2:Nx                                                    % looping
        for j=2:Ny
            Tss2(i,j)=0.25*(Tss(i+1,j)+Tss(i,j+1)+Tss(i-1,j)+Tss(i,j-1));
        end
    end
    k=k+1;                                                        % update k
    err_SS_max(k)=abs(max(max(Tss2-Tss)));                        % calculate error
    err_SS_min(k)=abs(min(min(Tss2-Tss)));                        % calculate error
    Tss=Tss2;                                                     % update T
    end
   

%% 4- Finite difference section (Using 2nd order Runge Kutte or Euler in time)

k=1;
switch Solution_type
    case '2nd order Runge-Kutte'
    err_R_k_max(k)=100;                            % initial error
    err_R_k_min(k)=100;                            % initial error
    while err_R_k_max(k)>=tolerence || err_R_k_min(k)>=tolerence
      for i=2:Nx
        for j=2:Ny
            k1=alpha*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
            Tk=T(:,:,k)+k1*dt;
            k2=alpha*(((Tk(i-1,j)-2*Tk(i,j)+Tk(i+1,j))/dx^2)+((Tk(i,j-1)-2*Tk(i,j)+Tk(i,j+1))/dy^2));
            T(i,j,k+1) =T(i,j,k)+(dt/2)*(k1+k2);
        end
      end
      k=k+1;
      err_R_k_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
      err_R_k_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
      if round(err_R_k_max(k),5)==round(err_R_k_max(k-1),5) && err_R_k_max(k)~= 0      % Test solution convergence
       errordlg('The solution is not converging, Please choose a larger tolerence','Tolerence Error');
       close(message)
       return
      end
      if round(err_R_k_min(k),5)==round(err_R_k_min(k-1),5) && err_R_k_min(k)~= 0      % Test solution convergence
       errordlg('The solution is not converging, Please choose a larger tolerence','Tolerence Error');
       close(message)
       return
      end
     end

    case'Euler'
    err_E_max(k)=100;                            % initial error
    err_E_min(k)=100;                            % initial error
    while err_E_max(k)>=tolerence || err_E_min(k)>=tolerence
      for i=2:Nx
        for j=2:Ny
            T(i,j,k+1) =T(i,j,k)+dt*alpha*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
        end
      end
      k=k+1;
      err_E_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
      err_E_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
      if round(err_E_max(k),5)==round(err_E_max(k-1),5) && err_E_max(k)~= 0      % Test solution convergence
       errordlg('The solution is not converging, Please choose a larger tolerence','Tolerence Error');
       close(message)
       return
      end
      if round(err_E_min(k),5)==round(err_E_min(k-1),5) && err_E_min(k)~= 0      % Test solution convergence
       errordlg('The solution is not converging, Please choose a larger tolerence','Tolerence Error');
       close(message)
       return
      end
     end

    case []
    close(message)
    msgbox('Error, Please re-run the code and choose Euler or 2nd order Runge-Kutte to continue the solution')
    return
end
T=T(:,:,1:k);                                            % delete the unused assigned zero layers
SStime=k*dt;                                             % steady state time
close(message)                                               % close the busy message


%% 5- Printed results section

fprintf('This is the solution of the heat equation through out a plate of diamensions %i X %i of material %s \n',Lx,Ly,name);
fprintf('The solution is  based on "Dirichlet Boundry Conditions" with initial values \n')
fprintf('T(x,0,t)=%i , T(x,%i,t)=%i , T(0,y,t)=%i , T(%i,y,t)=%i , T(x,y,0)=%i \n',T_south,Ly,T_north,T_west,Lx,T_east,T_initial)
fprintf('The plate takes %i seconds to reach steady-state temperature with tolerence %0.2f \n',round(SStime),tolerence);
fprintf('Now, Simulation is running with final time %i seconds and step %0.2f second \n',t_end,dt)


%% 6- Plotting section

x=zeros(1,Nx+2);y=zeros(1,Ny+2);            %Generate the plate
for i = 1:Nx+2                 
x(i) =(i-1)*dx; 
end
for i = 1:Ny+2                 
y(i) =(i-1)*dy; 
end

% %%%            -------------- Constant plot ----------------

subplot(2,2,3)                           
hold on
title(sprintf('Temperature at steady state time : %i seconds ',round(SStime)))
surf(x,y,Tss)
plot3(  Lx/4,  Ly/4,T_max,'ko','markerfacecolor','r') % plot red point
plot3(  Lx/2,  Ly/2,T_max,'ko','markerfacecolor','g') % plot green point
plot3(3*Lx/4,3*Ly/4,T_max,'ko','markerfacecolor','b') % plot blue point
plot3(  Lx/4,  Ly/4,T_min,'ko','markerfacecolor','r') % plot red point
plot3(  Lx/2,  Ly/2,T_min,'ko','markerfacecolor','g') % plot green point
plot3(3*Lx/4,3*Ly/4,T_min,'ko','markerfacecolor','b') % plot blue point
cb=colorbar;
caxis([T_min T_max]);
view(90,-90);
xlim([0 Lx+dx]); xlabel('Length');
ylim([0 Ly+dy]); ylabel('Width');
zlim([T_min T_max]); zlabel('Temprature');
drawnow
hold off

 subplot(2,2,4)
hold on
title(sprintf('Temperature at steady state time : %i seconds ',round(SStime)))
scatter(k,Tss(floor(Nx/4),floor(Ny/4)),'ko','markerfacecolor','r'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(Nx/4),floor(Ny/4))));
text(k,Tss(floor(Nx/4),floor(Ny/4)),val,'HorizontalAlignment','Left');
scatter(k,Tss(floor(Nx/2),floor(Ny/2)),'ko','markerfacecolor','g'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(Nx/2),floor(Ny/2))));
text(k,Tss(floor(Nx/2),floor(Ny/2)),val,'HorizontalAlignment','right');
scatter(k,Tss(floor(3*Nx/4),floor(3*Ny/4)),'ko','markerfacecolor','b'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(3*Nx/4),floor(3*Ny/4))));
text(k,Tss(floor(3*Nx/4),floor(3*Ny/4)),val,'HorizontalAlignment','Left');
axis tight; xlabel('Time Iterations');
ylim([T_min T_max]); ylabel('Temperature');
legend('Red Point','Green Point ','Blue Point ','Location','northwest')
drawnow
hold off

%%%             ------------ Animated plot ----------

for j=1:n_time                               
    
subplot(2,2,1)
surf(x,y,T(:,:,j))
hold on
title(sprintf('Temperature at time : %i seconds ',round(j*dt)))
plot3(  Lx/4,  Ly/4,T_max,'ko','markerfacecolor','r') % plot red point
plot3(  Lx/2,  Ly/2,T_max,'ko','markerfacecolor','g') % plot green point
plot3(3*Lx/4,3*Ly/4,T_max,'ko','markerfacecolor','b') % plot blue point
plot3(  Lx/4,  Ly/4,T_min,'ko','markerfacecolor','r') % plot red point
plot3(  Lx/2,  Ly/2,T_min,'ko','markerfacecolor','g') % plot green point
plot3(3*Lx/4,3*Ly/4,T_min,'ko','markerfacecolor','b') % plot blue point
cb=colorbar;
caxis([T_min T_max]);
view(90,-90);
xlim([0 Lx+dx]); xlabel('Length');
ylim([0 Ly+dy]); ylabel('Width');
zlim([T_min T_max]); zlabel('Temprature');
drawnow
hold off

subplot(2,2,2)
hold on
title(sprintf('Temperature at time : %i seconds ',round(j*dt)))
scatter(j,T(floor((Nx+2)/4),floor((Ny+2)/4),j),'r.'); 
scatter(j,T(ceil((Nx+2)/2),ceil((Ny+2)/2),j),'g.'); 
scatter(j,T(ceil(3*(Nx+2)/4),ceil(3*(Ny+2)/4),j),'b.'); 
axis tight; xlabel('Time Iterations');
axis tight; ylabel('Temperature');
legend('Red Point','Green Point ','Blue Point ','Location','northwest')
drawnow
hold off

end